#!/usr/bin/env python3
"""Coordinate Pi bias settings, PicoScope captures and automatic Vbr fits."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE.parent / "analysis"


def run(command: list[str], *, capture: bool = False) -> subprocess.CompletedProcess[str]:
    print("+ " + " ".join(shlex.quote(part) for part in command), flush=True)
    return subprocess.run(command, check=True, text=True, capture_output=capture)


def bias_tag(voltage: float) -> str:
    return f"bias_{voltage:06.3f}".replace(".", "p")


def ssh_command(host: str, remote_command: str) -> list[str]:
    return ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host, remote_command]


def command_pi_agent(config: dict, arguments: list[str], dry_run: bool) -> dict:
    common = [
        "--channels",
        ",".join(str(signal["detector_channel"]) for signal in config["signals"]),
        "--minimum-bias-v",
        str(config["minimum_bias_v"]),
        "--maximum-bias-v",
        str(config["maximum_bias_v"]),
    ]
    if dry_run:
        completed = run(
            [sys.executable, str(HERE / "pi_bias_control.py"), *arguments, *common, "--dry-run"],
            capture=True,
        )
    else:
        remote = " ".join(
            shlex.quote(part)
            for part in ["python3", config["pi_agent"], *arguments, *common]
        )
        completed = run(ssh_command(config["pi_host"], remote), capture=True)
    if completed.stdout:
        print(completed.stdout, end="")
    return json.loads(completed.stdout)


def stop_temperature_compensation(config: dict) -> None:
    command = "pkill -TERM -f '[b]iasAdj.py' || true; sleep 2; pgrep -af '[b]iasAdj.py' && exit 1 || true"
    run(ssh_command(config["pi_host"], command))


def start_temperature_compensation(config: dict) -> None:
    command = (
        "cd /home/cosmic && "
        "nohup python3 /home/cosmic/biasAdj.py "
        ">/home/cosmic/biasAdj_scan_restore.log 2>&1 </dev/null &"
    )
    run(ssh_command(config["pi_host"], command))


def probe_scope(config: dict) -> None:
    scope = config["scope"]
    run(
        [
            sys.executable,
            str(HERE / "pico_capture.py"),
            "--probe",
            "--driver",
            scope["driver"],
            "--range-mv",
            str(scope["range_mv"]),
            "--probe-attenuation",
            str(scope.get("probe_attenuation", 1.0)),
        ]
    )


def capture_point(
    config: dict,
    signal: dict,
    raw_dir: Path,
    voltage: float,
    events: int,
    synthetic: bool,
) -> None:
    scope = config["scope"]
    command = [
        sys.executable,
        str(HERE / "pico_capture.py"),
        "--output",
        str(raw_dir),
        "--events",
        str(events),
        "--bias-v",
        str(voltage),
        "--driver",
        scope["driver"],
        "--range-mv",
        str(scope["range_mv"]),
        "--probe-attenuation",
        str(scope.get("probe_attenuation", 1.0)),
        "--sample-interval-ns",
        str(scope["sample_interval_ns"]),
        "--pre-trigger-ns",
        str(scope["pre_trigger_ns"]),
        "--post-trigger-ns",
        str(scope["post_trigger_ns"]),
        "--trigger-channel",
        signal["scope_channel"],
        "--trigger-mv",
        str(signal.get("trigger_mv", scope["trigger_mv"])),
        "--auto-trigger-ms",
        str(scope["auto_trigger_ms"]),
    ]
    if synthetic:
        by_scope = {signal["scope_channel"]: signal for signal in config["signals"]}
        signal_a = by_scope["A"]
        signal_b = by_scope["B"]
        command.extend(
            [
                "--synthetic",
                "--synthetic-vbr-a",
                str(signal_a["synthetic_vbr_v"]),
                "--synthetic-vbr-b",
                str(signal_b["synthetic_vbr_v"]),
                "--synthetic-gain-a",
                str(signal_a["synthetic_gain_mv_per_v"]),
                "--synthetic-gain-b",
                str(signal_b["synthetic_gain_mv_per_v"]),
                "--seed",
                str(round(voltage * 1000)),
            ]
        )
    run(command)


def analyze_signal(signal: dict, raw_dir: Path, output_dir: Path, voltage: float) -> None:
    run(
        [
            sys.executable,
            str(ANALYSIS_DIR / "analyze_point.py"),
            "--run-dir",
            str(raw_dir),
            "--output",
            str(output_dir),
            "--bias-v",
            str(voltage),
            "--sipm",
            signal["name"],
            "--detector-channel",
            str(signal["detector_channel"]),
            "--scope-channel",
            signal["scope_channel"],
            "--color",
            signal["color"],
        ]
    )


def fit_results(config: dict, run_root: Path) -> None:
    for signal in config["signals"]:
        for metric in ("height", "area"):
            run(
                [
                    sys.executable,
                    str(ANALYSIS_DIR / "fit_vbr.py"),
                    "--run-root",
                    str(run_root),
                    "--signal",
                    signal["name"],
                    "--color",
                    signal["color"],
                    "--metric",
                    metric,
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path.home() / "brDownVstudy" / "automated_scans")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="synthetic waveforms and no hardware access")
    mode.add_argument("--live", action="store_true", help="control real Pi and PicoScope hardware")
    parser.add_argument("--confirm", default="")
    parser.add_argument("--events", type=int, help="override events_per_point")
    parser.add_argument(
        "--skip-fit",
        action="store_true",
        help="keep all acquired points without requiring enough resolved peaks for an immediate Vbr fit",
    )
    parser.add_argument(
        "--skip-analysis",
        action="store_true",
        help="acquire every configured point before running any waveform analysis",
    )
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.live and args.confirm != "ENABLE_HV_SCAN":
        parser.error("live operation requires --confirm ENABLE_HV_SCAN")
    scope_channels = [item["scope_channel"] for item in config["signals"]]
    if not 1 <= len(scope_channels) <= 2 or len(set(scope_channels)) != len(scope_channels):
        parser.error("configuration must define one or two signals on distinct scope channels")
    if args.dry_run and set(scope_channels) != {"A", "B"}:
        parser.error("synthetic mode requires exactly one signal on scope channel A and one on B")

    events = args.events or int(config["events_per_point"])
    start_tag = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = args.output_root / f"vbr_scan_{start_tag}"
    run_root.mkdir(parents=True, exist_ok=False)
    (run_root / "scan_config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    manifest: dict[str, object] = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "dry-run" if args.dry_run else "live",
        "events_per_point": events,
        "points": [],
        "completed": False,
    }
    manifest_path = run_root / "scan_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    live_success = False

    if args.live:
        try:
            probe_scope(config)
            manifest["scope_preflight"] = "passed"
        except Exception as exc:
            manifest["scope_preflight"] = f"failed before hardware control: {exc}"
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            raise

    try:
        if args.live:
            stop_temperature_compensation(config)

        for voltage in config["bias_points_v"]:
            tag = bias_tag(float(voltage))
            point: dict[str, object] = {"requested_bias_v": voltage, "tag": tag, "status": "started"}
            manifest["points"].append(point)
            point["bias_control"] = command_pi_agent(
                config,
                ["set", "--voltage", str(voltage)],
                args.dry_run,
            )
            effective_bias_v = float(point["bias_control"]["plan"]["effective_bias_v"])
            point["effective_bias_v"] = effective_bias_v
            if args.live:
                time.sleep(float(config["settle_seconds"]))
            point["environment"] = command_pi_agent(config, ["environment"], args.dry_run)

            for signal in config["signals"]:
                raw_dir = run_root / "raw" / tag / signal["name"]
                analysis_dir = run_root / "analysis" / tag / signal["name"]
                capture_point(
                    config,
                    signal,
                    raw_dir,
                    effective_bias_v,
                    events,
                    args.dry_run,
                )
                if not args.skip_analysis:
                    analyze_signal(signal, raw_dir, analysis_dir, effective_bias_v)
            point["status"] = "complete"
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

        if not args.skip_fit and not args.skip_analysis:
            fit_results(config, run_root)
        manifest["completed"] = True
        manifest["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        live_success = True
    finally:
        active_scan_exception = sys.exc_info()[0] is not None
        finalization_error: Exception | None = None
        if args.live:
            try:
                if live_success and config.get("restore_normal_after_scan", False):
                    command_pi_agent(config, ["normal"], False)
                    start_temperature_compensation(config)
                    manifest["final_hardware_state"] = "normal bias restored; temperature compensation restarted"
                else:
                    command_pi_agent(config, ["off"], False)
                    manifest["final_hardware_state"] = "HV off"
            except Exception as exc:
                finalization_error = exc
                try:
                    command_pi_agent(config, ["off"], False)
                    manifest["final_hardware_state"] = f"HV off after finalization failure: {exc}"
                except Exception as shutdown_exc:
                    manifest["final_hardware_state"] = (
                        f"FINAL HARDWARE STATE UNKNOWN; finalization failed: {exc}; "
                        f"HV-off command also failed: {shutdown_exc}"
                    )
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        if finalization_error is not None and not active_scan_exception:
            raise RuntimeError(manifest["final_hardware_state"]) from finalization_error

    print(f"Scan complete: {run_root}")


if __name__ == "__main__":
    main()
