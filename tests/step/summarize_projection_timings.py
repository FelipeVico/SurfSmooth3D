#!/usr/bin/env python3
"""Summarize seven-repeat, four-process ABBA full-mesh timing captures.

The first repetition of every fixture in each process is excluded. Exact
payload comparisons remain required, including all excluded warm-up runs.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import statistics

from summarize_projection_benchmark import SMALL_CASES, Pending, file_record, read_json, rows


PROCESS_ORDER = ("baseline-a", "candidate-a", "candidate-b", "baseline-b")


def distribution(values):
    ordered = sorted(values)

    def quantile(fraction):
        position = fraction * (len(ordered)-1)
        lower = math.floor(position)
        upper = math.ceil(position)
        return ordered[lower] + (position-lower) * (ordered[upper]-ordered[lower])

    return {"count": len(values), "median_seconds": statistics.median(values),
            "minimum_seconds": min(values), "maximum_seconds": max(values),
            "p25_seconds": quantile(0.25), "p75_seconds": quantile(0.75),
            "samples_seconds": values}


def build_report(directory, root, background_load, expected_cases=SMALL_CASES):
    expected_cases = set(expected_cases)
    if not expected_cases or not expected_cases.issubset(SMALL_CASES):
        raise ValueError("Selected timing fixtures must be a nonempty subset of the small matrix")
    captures = {}
    evidence = []
    for label in PROCESS_ORDER:
        path = directory / ("timing-" + label) / "capture.json"
        capture = read_json(path)
        expected_package = root / "build" / ("fast-baseline" if label.startswith("baseline") else "fast-candidate")
        if Path(capture.get("package_root", "")).resolve() != expected_package.resolve():
            raise ValueError(f"Wrong native package for timing-{label}")
        expected_mex = capture.get("resolved_functions", {}).get("expected_mex")
        mex_path = Path(expected_mex or "")
        if (mex_path.parent != expected_package / "matlab/+surfsmooth3d/+stepmesher/private" or
                not mex_path.name.startswith("step_mesher_mex.mex")):
            raise ValueError(f"Wrong expected private MEX for timing-{label}")
        cases = rows(capture.get("cases"))
        by_name = {case["name"]: case for case in cases}
        if set(by_name) != expected_cases or len(cases) != len(expected_cases):
            raise Pending(f"timing-{label} has not captured every expected fixture")
        captures[label] = by_name
        for case in cases:
            if case.get("passed") is not True:
                raise ValueError(f"timing-{label} fixture failed: {case['name']}")
            runs = rows(case.get("runs"))
            if len(runs) != 7:
                raise Pending(f"timing-{label} {case['name']} has {len(runs)} of seven repetitions")
            if case.get("profile_choice") != "rigid" or case.get("options", {}).get("order") != 4:
                raise ValueError(f"Timing fixture configuration changed: {case['name']}")
            for run in runs:
                duration = run.get("mesh_seconds")
                if (not isinstance(duration, (int, float)) or not math.isfinite(duration) or duration <= 0 or
                        run.get("loaded_mex") != expected_mex):
                    raise ValueError(f"Invalid timing or loaded MEX in {label} {case['name']}")
        evidence.append(file_record(path, root))
        mat_file = path.with_suffix(".mat")
        if not mat_file.is_file():
            raise Pending(f"Missing full timing capture {mat_file}")
        evidence.append(file_record(mat_file, root))
    comparisons = []
    for pair in ("a", "b"):
        path = directory / f"timing-comparison-{pair}.json"
        comparison = read_json(path)
        cases = rows(comparison.get("cases"))
        exact = (comparison.get("passed") is True and comparison.get("bitwise_passed") is True and
                 {case["id"].removesuffix("_p4_rigid") for case in cases} == expected_cases and
                 len(cases) == len(expected_cases))
        exact = exact and Path(comparison.get("baseline_file", "")).resolve() == (
            directory / f"timing-baseline-{pair}/capture.mat").resolve()
        exact = exact and Path(comparison.get("candidate_file", "")).resolve() == (
            directory / f"timing-candidate-{pair}/capture.mat").resolve()
        for case in cases:
            fields = rows(case.get("comparisons"))
            exact = exact and case.get("passed") is True and case.get("bitwise_passed") is True
            exact = exact and bool(fields) and all(field.get("isequaln") is True and field.get("bitwise") is True
                                                  for field in fields)
            exact = exact and len(rows(case.get("baseline_seconds"))) == 7 and len(rows(case.get("candidate_seconds"))) == 7
            for variant in ("baseline", "candidate"):
                compared_repeats = {int(match[1]) for field in fields
                                    if (match := re.match(variant + r"_repeat_(\d+)\.", field["label"]))}
                exact = exact and compared_repeats == set(range(1, 8))
        if not exact:
            raise ValueError(f"timing-comparison-{pair} has a non-exact payload comparison")
        comparisons.append(file_record(path, root))
    result = []
    for name in sorted(expected_cases):
        values = {}
        process_medians = {}
        for label in PROCESS_ORDER:
            runs = rows(captures[label][name]["runs"])
            values[label] = [run["mesh_seconds"] for run in runs[1:]]
            process_medians[label] = statistics.median(values[label])
        baseline = distribution(values["baseline-a"] + values["baseline-b"])
        candidate = distribution(values["candidate-a"] + values["candidate-b"])
        delta = candidate["median_seconds"] - baseline["median_seconds"]
        fraction = delta / baseline["median_seconds"]
        flag = candidate["median_seconds"] > 1.05 * baseline["median_seconds"] and delta > 0.001
        result.append({"name": name, "baseline": baseline, "candidate": candidate,
                       "per_process_median_seconds": process_medians,
                       "median_speedup": baseline["median_seconds"] / candidate["median_seconds"],
                       "candidate_median_change_seconds": delta,
                       "candidate_median_change_percent": fraction * 100,
                       "flag_for_quiet_rerun": flag})
    flags = [case["name"] for case in result if case["flag_for_quiet_rerun"]]
    return {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
            "status": "pending" if flags else "passed", "completed": not flags, "passed": not flags,
            "measurements_completed": True, "exact_payloads_passed": True,
            "expected_cases": sorted(expected_cases),
            "process_order": list(PROCESS_ORDER), "repetitions_per_fixture_per_process": 7,
            "excluded_repetitions_per_fixture_per_process": [1], "warm_samples_per_variant_per_fixture": 12,
            "capture_evidence": evidence, "comparison_evidence": comparisons, "cases": result,
            "quiet_rerun_required_cases": flags,
            "regression_threshold": {"relative_median_increase": 0.05, "absolute_median_increase_seconds": 0.001,
                                     "condition": "Both thresholds must be exceeded."},
            "background_load": background_load,
            "timing_caveat": "ABBA order and warm repetitions reduce startup/order bias. Background activity can still affect these measurements; any threshold violation needs a quiet rerun before the performance gate passes."}


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-directory", type=Path, default=root / "build/validation/fast")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--selected-cases", nargs="+", default=sorted(SMALL_CASES))
    parser.add_argument("--background-load", default="Six projection replay workers and an existing user MATLAB process were active; only one timing MATLAB worker ran at a time.")
    args = parser.parse_args()
    directory = args.evidence_directory.resolve()
    output = args.output or directory / "repeated-timing-summary.json"
    try:
        report = build_report(directory, root, args.background_load, args.selected_cases)
    except (Pending, ValueError, KeyError, TypeError, OSError) as error:
        report = {"status": "pending" if isinstance(error, Pending) else "failed",
                  "completed": False, "passed": False, "reason": str(error)}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)
    print(json.dumps({"output": str(output), "status": report["status"],
                      "quiet_rerun_required_cases": report.get("quiet_rerun_required_cases", [])}))


if __name__ == "__main__":
    main()
