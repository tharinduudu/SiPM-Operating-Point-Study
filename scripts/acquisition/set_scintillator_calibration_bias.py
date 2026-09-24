#!/usr/bin/env python3
"""Set all four SiPM biases for the EPIC-tile-trigger calibration setup."""

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
    parser.add_argument("--epic-vbr-v", type=float, default=50.543)
    parser.add_argument("--epic-vbr-reference-temperature-c", type=float, default=20.4)
    parser.add_argument("--sipm1-vbr-20p4c", type=float, default=50.512)
    parser.add_argument("--sipm2-vbr-20p4c", type=float, default=50.574)
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

    # The EPIC-tile SiPMs on CH0 and CH1 have no device-specific Vbr measurement yet. Their values
    # use the mean measured Vbr of the SiPM 1 and SiPM 2 devices on this same
    # readout. This is a conservative board-referenced proxy, not a substitute
    # for measuring the two EPIC-tile SiPMs.
    epic_target = corrected_bias(
        args.epic_vbr_v,
        args.epic_vbr_reference_temperature_c,
        temperature,
        args.temperature_coefficient_v_per_c,
        args.overvoltage_v,
    )
    targets = {
        0: epic_target,
        1: epic_target,
        2: corrected_bias(
            args.sipm1_vbr_20p4c,
            20.4,
            temperature,
            args.temperature_coefficient_v_per_c,
            args.overvoltage_v,
        ),
        3: corrected_bias(
            args.sipm2_vbr_20p4c,
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
                "label": "top_EPIC_tile_reference",
                "vbr_V": args.epic_vbr_v,
                "vbr_reference_temperature_C": args.epic_vbr_reference_temperature_c,
                "vbr_status": "proxy: mean measured Vbr of SiPM 1 and SiPM 2; individual Vbr not measured",
            },
            "1": {
                "label": "bottom_EPIC_tile_reference",
                "vbr_V": args.epic_vbr_v,
                "vbr_reference_temperature_C": args.epic_vbr_reference_temperature_c,
                "vbr_status": "proxy: mean measured Vbr of SiPM 1 and SiPM 2; individual Vbr not measured",
            },
            "2": {
                "label": "SiPM_1_(△)_top_GSU_gLOWCOST_tile",
                "vbr_V": args.sipm1_vbr_20p4c,
                "vbr_reference_temperature_C": 20.4,
                "vbr_status": "measured",
            },
            "3": {
                "label": "SiPM_2_(★)_bottom_GSU_gLOWCOST_tile",
                "vbr_V": args.sipm2_vbr_20p4c,
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
