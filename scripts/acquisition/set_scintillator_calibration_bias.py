#!/usr/bin/env python3
"""Set all four SiPM biases for the small-tile-trigger calibration setup."""

from __future__ import annotations

import argparse
import fcntl
import json
from datetime import datetime, timezone

try:
    from pi_bias_control import (
        DEFAULT_DAC_SCRIPT,
        DEFAULT_LOCK_FILE,
        DEFAULT_LOG_FILE,
        DEFAULT_MAX_BINARY,
        DEFAULT_STATE_FILE,
        read_environment,
        switch_off,
    )
except ModuleNotFoundError:
    from vbr_bias_control import (
        DEFAULT_DAC_SCRIPT,
        DEFAULT_LOCK_FILE,
        DEFAULT_LOG_FILE,
        DEFAULT_MAX_BINARY,
        DEFAULT_STATE_FILE,
        read_environment,
        switch_off,
    )

from set_dual_operating_point import apply_plan, choose_plan


def corrected_bias(vbr: float, reference_c: float, temperature_c: float, coefficient: float,
                   overvoltage: float) -> float:
    return vbr + coefficient * (temperature_c - reference_c) + overvoltage


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small-vbr-v", type=float, default=50.543)
    parser.add_argument("--small-vbr-reference-temperature-c", type=float, default=20.4)
    parser.add_argument("--triangle-vbr-20p4c", type=float, default=50.512)
    parser.add_argument("--star-vbr-20p4c", type=float, default=50.574)
    parser.add_argument("--overvoltage-v", type=float, default=3.0)
    parser.add_argument("--temperature-coefficient-v-per-c", type=float, default=0.054)
    parser.add_argument("--ramp-step-codes", type=int, default=32)
    parser.add_argument("--ramp-delay-s", type=float, default=0.08)
    parser.add_argument("--confirm", required=True)
    args = parser.parse_args()
    if args.confirm != "SET_FOUR_BIASES":
        parser.error("live operation requires --confirm SET_FOUR_BIASES")

    environment = read_environment(False)
    temperature = environment.get("temperature_C")
    if temperature is None:
        raise RuntimeError(f"BME280 temperature unavailable: {environment}")
    temperature = float(temperature)

    # CH0 and CH1 have no device-specific Vbr measurement yet. Their values
    # use the mean measured Vbr of the Triangle and Star devices on this same
    # readout. This is a conservative board-referenced proxy, not a substitute
    # for measuring the two reference SiPMs.
    small_target = corrected_bias(
        args.small_vbr_v,
        args.small_vbr_reference_temperature_c,
        temperature,
        args.temperature_coefficient_v_per_c,
        args.overvoltage_v,
    )
    targets = {
        0: small_target,
        1: small_target,
        2: corrected_bias(
            args.triangle_vbr_20p4c,
            20.4,
            temperature,
            args.temperature_coefficient_v_per_c,
            args.overvoltage_v,
        ),
        3: corrected_bias(
            args.star_vbr_20p4c,
            20.4,
            temperature,
            args.temperature_coefficient_v_per_c,
            args.overvoltage_v,
        ),
    }
    plan = choose_plan(targets)
    record: dict[str, object] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "action": "set_scintillator_calibration_biases",
        "environment": environment,
        "temperature_coefficient_V_per_C": args.temperature_coefficient_v_per_c,
        "overvoltage_V": args.overvoltage_v,
        "channel_mapping": {
            "0": {
                "label": "small_top_reference",
                "vbr_V": args.small_vbr_v,
                "vbr_reference_temperature_C": args.small_vbr_reference_temperature_c,
                "vbr_status": "proxy: mean measured Vbr of Triangle and Star; individual Vbr not measured",
            },
            "1": {
                "label": "small_bottom_reference",
                "vbr_V": args.small_vbr_v,
                "vbr_reference_temperature_C": args.small_vbr_reference_temperature_c,
                "vbr_status": "proxy: mean measured Vbr of Triangle and Star; individual Vbr not measured",
            },
            "2": {
                "label": "Triangle_SIPM1_top_large",
                "vbr_V": args.triangle_vbr_20p4c,
                "vbr_reference_temperature_C": 20.4,
                "vbr_status": "measured",
            },
            "3": {
                "label": "Star_SIPM2_bottom_large",
                "vbr_V": args.star_vbr_20p4c,
                "vbr_reference_temperature_C": 20.4,
                "vbr_status": "measured",
            },
        },
        "plan": plan,
    }

    DEFAULT_LOCK_FILE.parent.mkdir(parents=True, exist_ok=True)
    with DEFAULT_LOCK_FILE.open("w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        try:
            record["commands"] = apply_plan(
                plan,
                ramp_step=args.ramp_step_codes,
                ramp_delay=args.ramp_delay_s,
            )
        except Exception:
            switch_off(
                sorted(targets),
                max_binary=DEFAULT_MAX_BINARY,
                dac_script=DEFAULT_DAC_SCRIPT,
                dry_run=False,
            )
            raise
        DEFAULT_STATE_FILE.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        with DEFAULT_LOG_FILE.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
