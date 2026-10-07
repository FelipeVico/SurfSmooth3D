#!/usr/bin/env python3
"""Summarize projector-reuse validation without running any benchmarks.

Missing or unfinished evidence stays pending. Overall completed/passed become
true only when every required check passes. Concurrent replay timings are
reported as diagnostics and never called a serial full-pipeline speedup.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys


BASELINE_COMMIT = "562f968"
HAIRY_STEM = "wobbly_hairy_torus_10_v2_flat_AP214"
HAIRY_PREFIX = HAIRY_STEM + "_p4"
SMALL_CASES = {
    "test_geometry_1", "test_curved_uv", "spheres_intersect", "Multi_res",
    "UV_curved_bezier", "piecewise_bezier_cap", "torus_ellipse",
}
NATIVE_CASES = SMALL_CASES | {
    "Multiscale_step", "cylinder_bezier_v2", "intersection_two_balls",
    "test_geometry_2", "test_geometry_manas", HAIRY_STEM,
}
EXTRACTION_CASES = {
    "test_geometry_1", "test_curved_uv", "UV_curved_bezier", "torus_ellipse",
}


class Pending(Exception):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path, root):
    path = Path(path)
    try:
        label = str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        label = str(path.resolve())
    result = {"path": label, "exists": path.is_file()}
    if result["exists"]:
        result.update(bytes=path.stat().st_size, sha256=sha256(path))
    return result


def read_json(path):
    if not path.is_file():
        raise Pending(f"Missing {path}")
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise Pending(f"Unfinished or invalid JSON in {path}: {error}") from error


def rows(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def section(status, **values):
    return {"status": status, "completed": status == "passed", **values}


def guarded(operation):
    try:
        return operation()
    except Pending as error:
        return section("pending", reason=str(error))
    except (ValueError, KeyError, TypeError, OSError) as error:
        return section("failed", reason=str(error))


def native_tests(directory, root):
    path = directory / "ctest-candidate.log"
    if not path.is_file():
        raise Pending("Candidate CTest log is missing")
    text = path.read_text()
    summaries = re.findall(r"(\d+)% tests passed, (\d+) tests failed out of (\d+)", text)
    if not summaries:
        raise Pending("Candidate CTest final summary has not been written")
    percent, failed, total = map(int, summaries[-1])
    records = re.findall(r"\d+/\d+\s+Test\s+#\d+:\s+(.+?)\s+\.{2,}\s+(Passed|\*\*\*Failed)", text)
    passed_names = [name.strip() for name, status in records if status == "Passed"]
    required_groups = {
        "projection_reuse": "projection_reuse" in passed_names,
        "frozen_project_nodes": "frozen_project-nodes" in passed_names,
        "retention_project_nodes": "retention_project-nodes" in passed_names,
        "context_lifetime": "occ_context" in passed_names,
    }
    passed = percent == 100 and failed == 0 and total > 0 and all(required_groups.values())
    return section("passed" if passed else "failed", evidence=file_record(path, root),
                   passed_tests=total-failed, failed_tests=failed, total_tests=total,
                   required_groups=required_groups, passed_test_names=passed_names)


def native_projection_matrix(directory, root):
    path = directory / "native-matrix/exit-statuses.json"
    exits = rows(read_json(path))
    passed = ({row["fixture"] for row in exits} == NATIVE_CASES and len(exits) == len(NATIVE_CASES))
    fixtures = []
    required_cases = {"contracts", "all-faces A", "all-faces shuffled B", "all-faces repeated A"} | {
        f"primary batch {count} repetition 1" for count in (1, 2, 16)}
    for exit_row in sorted(exits, key=lambda value: value["fixture"]):
        name = exit_row["fixture"]
        stdout = Path(exit_row["stdout"])
        if not stdout.is_file():
            raise Pending(f"Native projection matrix stdout is missing for {name}")
        records = [json.loads(line) for line in stdout.read_text().splitlines() if line.startswith("{")]
        timed = [row for row in records if "oracle_digest" in row]
        names = [row.get("case") for row in records]
        contracts = [row for row in records if row.get("case") == "contracts"]
        source = root / "step_mesher/examples/step_files" / (name + ".step")
        fixture_passed = (source.is_file() and exit_row.get("exit_status") == 0 and required_cases.issubset(names) and
                          len(names) == len(set(names)) and bool(records) and len(timed) >= 6 and
                          all(row.get("exact") is True and Path(row.get("fixture", "")).resolve() == source.resolve()
                              for row in records) and
                          all(bool(row.get("oracle_digest")) and row.get("oracle_digest") == row.get("candidate_digest")
                              and row.get("nodes", 0) > 0 for row in timed) and
                          len(contracts) == 1 and all(contracts[0].get(key) is True for key in
                              ("empty_batch", "unknown_face_rejected", "independent_contexts")))
        for count in (1, 2, 16):
            fixture_passed = fixture_passed and any(row.get("case") == f"primary batch {count} repetition 1" and
                                                   row.get("nodes") == count for row in timed)
        fixtures.append({"fixture": name, "passed": fixture_passed, "exit_status": exit_row.get("exit_status"),
                         "projection_comparisons": len(timed), "compared_nodes_including_repeated_batches":
                         sum(row.get("nodes", 0) for row in timed), "records": records,
                         "stdout_evidence": file_record(stdout, root), "fixture_evidence": file_record(source, root),
                         "stderr_evidence": file_record(Path(exit_row["stderr"]), root)})
        passed = passed and fixture_passed
    return section("passed" if passed else "failed", evidence=file_record(path, root),
                   expected_fixtures=sorted(NATIVE_CASES), tested_fixtures=len(fixtures),
                   projection_comparisons=sum(row["projection_comparisons"] for row in fixtures),
                   fixtures=fixtures,
                   scope="Native projection comparisons against the preserved original projector, beyond the separate CTest import fixture coverage. Synthetic UV/seam/narrow-span, grouped/shuffled/repeated and 1/2/16-node batches; full actual target nodes are checked by replay.")


def small_comparison(directory, root):
    path = directory / "small-comparison.json"
    data = read_json(path)
    cases = rows(data.get("cases"))
    names = {entry["id"].removesuffix("_p4_rigid") for entry in cases}
    summary = []
    checks = [data.get("passed") is True, data.get("bitwise_passed") is True,
              names == SMALL_CASES, len(cases) == len(SMALL_CASES)]
    for entry in cases:
        comparisons = rows(entry.get("comparisons"))
        exact = bool(comparisons) and all(row.get("isequaln") is True and
                                          row.get("bitwise") is True for row in comparisons)
        checks.append(exact and entry.get("passed") is True and entry.get("bitwise_passed") is True)
        summary.append({key: entry.get(key) for key in (
            "id", "passed", "bitwise_passed", "baseline_seconds", "candidate_seconds",
            "baseline_median_seconds", "candidate_median_seconds", "speedup")})
        summary[-1]["comparisons"] = len(comparisons)
        summary[-1]["changed_fields"] = [row["label"] for row in comparisons
                                           if not (row.get("isequaln") and row.get("bitwise"))]
    return section("passed" if all(checks) else "failed", evidence=file_record(path, root),
                   expected_cases=sorted(SMALL_CASES), cases=summary,
                   timing_caveat="The stored full-mesh measurements are single runs, not a statistical no-regression proof.")


def extraction_parity(directory, root):
    path = directory / "extraction-parity" / "parity.json"
    data = read_json(path)
    cases = rows(data.get("cases"))
    passed = data.get("passed") is True and {row["name"] for row in cases} == EXTRACTION_CASES
    passed = passed and len(cases) == len(EXTRACTION_CASES) and all(
        row.get("numeric_comparisons", 0) > 0 and
        row.get("exact_comparisons") == row.get("numeric_comparisons") and
        row.get("maximum_scaled_difference") == 0 for row in cases)
    candidate = (root / "build/fast-candidate").resolve()
    matlab_root = candidate / "matlab"
    resolved = data.get("resolved_functions", {})
    names = rows(resolved.get("names"))
    paths = rows(resolved.get("paths"))
    required_functions = {
        "surfsmooth3d.stepmesher.mesh_step", "surfsmooth3d.stepmesher.to_srcvals",
        "surfsmooth3d.stepmesher.validate_area_flux", "surfsmooth3d.stepmesher.export_geom_package",
        "surfsmooth3d.surfer",
    }
    expected_mex = Path(resolved.get("expected_mex", ""))
    route_checks = {
        "route_verified": data.get("route_verified") is True,
        "candidate_package": Path(data.get("package_root", "")).resolve() == candidate,
        "live_fixtures": Path(data.get("fixture_root", "")).resolve() == root.resolve(),
        "all_requested_functions": set(names) == required_functions and len(paths) == len(names),
        "functions_in_candidate_package": bool(paths) and all(
            Path(value).resolve().is_relative_to(matlab_root) for value in paths),
        "private_candidate_mex": expected_mex.parent == matlab_root / "+surfsmooth3d/+stepmesher/private"
            and expected_mex.name.startswith("step_mesher_mex.mex"),
        "loaded_mex_matches_candidate": data.get("loaded_mex") == str(expected_mex),
        "frozen_baseline_exists": Path(data.get("baseline_file", "")).is_file(),
    }
    passed = passed and all(route_checks.values())
    return section("passed" if passed else "failed", evidence=file_record(path, root),
                   route_checks=route_checks, package_root=data.get("package_root"),
                   loaded_mex=data.get("loaded_mex"), baseline=file_record(Path(data.get("baseline_file", "")), root),
                   native_module=file_record(expected_mex, root),
                   cases=cases, comparison="Exact MATLAB numerical equality, including the existing GO3 round-trip checks; byte equality is separately checked by the differential captures.")


def candidate_smoke(directory, root):
    path = directory / "candidate-smoke.json"
    data = read_json(path)
    candidate = (root / "build/fast-candidate").resolve()
    loaded_mex = Path(data.get("loaded_mex", ""))
    passed = (data.get("passed") is True and data.get("route_verified") is True and
              Path(data.get("package_root", "")).resolve() == candidate and
              loaded_mex.parent == candidate / "matlab/+surfsmooth3d/+stepmesher/private" and
              loaded_mex.name.startswith("step_mesher_mex.mex"))
    return section("passed" if passed else "failed", evidence=file_record(path, root),
                   package_root=data.get("package_root"), loaded_mex=data.get("loaded_mex"),
                   fixture_source=data.get("fixture_source"), route_verified=data.get("route_verified"))


def request_metrics(request, mapping):
    canonical = {}
    for line in mapping.read_text().splitlines():
        tag, index = map(int, line.split())
        if tag in canonical:
            raise ValueError("Duplicate face tag in captured map")
        canonical[tag] = index
    counts = Counter()
    with request.open() as handle:
        if handle.readline().strip() != "STEP_MESHER_OCC_PROJECT_V1":
            raise ValueError("Invalid captured request marker")
        step_file = handle.readline().strip()
        face_count = int(handle.readline())
        descriptor_tags = []
        for _ in range(face_count):
            descriptor = handle.readline().split()
            if len(descriptor) != 8:
                raise ValueError("Invalid face descriptor in captured request")
            descriptor_tags.append(int(descriptor[0]))
        source_nodes = int(handle.readline())
        read_nodes = 0
        for line in handle:
            point = line.split()
            if len(point) != 4 or not all(math.isfinite(float(x)) for x in point[1:]):
                raise ValueError("Invalid or nonfinite captured projection input")
            tag = int(point[0])
            if tag not in canonical:
                raise ValueError("Captured projection point lacks a face association")
            counts[tag] += 1
            read_nodes += 1
    if read_nodes != source_nodes or source_nodes <= 0:
        raise ValueError("Captured request node count differs from its payload")
    if (set(descriptor_tags) != set(canonical) or len(canonical) != face_count or
            sorted(canonical.values()) != list(range(face_count))):
        raise ValueError("Captured face descriptors and canonical map differ")
    result = {"step_file": step_file, "source_nodes": source_nodes,
              "source_faces": face_count, "represented_faces": len(counts),
              "per_source_face": [{"tag": tag, "canonical_index": canonical[tag], "nodes": counts[tag]}
                                  for tag in sorted(canonical, key=canonical.get)]}
    if Path(step_file).stem == HAIRY_STEM:
        result["main_spline_canonical_index"] = 0
        result["main_spline_nodes"] = sum(counts[tag] for tag in canonical if canonical[tag] == 0)
    return result


def native_timings(path):
    if not path.is_file():
        raise Pending("Candidate hairy-torus stage CSV is missing")
    stages = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            seconds = float(row["seconds"])
            if row["stage"] in stages or not math.isfinite(seconds) or seconds < 0:
                raise ValueError("Invalid native stage timings")
            stages[row["stage"]] = seconds
    return stages


def candidate_hairy(directory, root):
    capture_path = directory / "capture-hairy" / "capture.json"
    data = read_json(capture_path)
    entries = rows(data.get("cases"))
    matching = [entry for entry in entries if entry.get("name") == HAIRY_STEM]
    if len(matching) != 1:
        raise ValueError("Hairy-torus full-mesh capture does not contain exactly one target case")
    entry = matching[0]
    if entry.get("passed") is not True:
        return section("failed", evidence=file_record(capture_path, root), reason=entry.get("error_message", "Hairy-torus capture failed"))
    runs = rows(entry.get("runs"))
    if not runs:
        raise Pending("Hairy-torus full-mesh capture has no completed run")
    run = runs[-1]
    prefix = directory / "capture-hairy" / HAIRY_PREFIX
    stages_path = Path(str(prefix) + ".timings.csv")
    stages = native_timings(stages_path)
    log = directory / "capture-hairy.log"
    if (not log.is_file() or "PROJECTION_BENCHMARK_COMPLETE" not in log.read_text() or
            "projection" not in stages or "result_packaging" not in stages):
        raise Pending("Hairy-torus full-mesh run has not reached its completion marker and final native stages")
    metrics = request_metrics(Path(str(prefix) + ".request.txt"), Path(str(prefix) + ".map.txt"))
    required_options = {"order": 4, "edge_gll_order": 4, "mesh_fraction": 0.10,
                        "refinement_level": 0, "occt_precision": 1e-6,
                        "occt_maxprecision": 1e-6, "sameparameter": True,
                        "surfacecurve_mode": "3d_preferred"}
    options = entry.get("options", {})
    passed = all(options.get(key) == value for key, value in required_options.items())
    passed = passed and entry.get("profile_choice") == "curvature_balanced_local"
    passed = passed and run.get("nodes") == metrics["source_nodes"] and run.get("projection_failures") == 0
    passed = passed and run.get("loaded_mex") == data["resolved_functions"]["expected_mex"]
    flux = run.get("normal_flux")
    return section("passed" if passed else "failed", evidence=file_record(capture_path, root),
                   stage_evidence=file_record(stages_path, root), options=options,
                   profile_choice=entry.get("profile_choice"), capture_metrics=metrics,
                   diagnostics={key: run.get(key) for key in (
                       "mesh_seconds", "srcvals_seconds", "validation_seconds", "patches", "nodes",
                       "projection_failures", "area", "relative_area_error", "normal_flux",
                       "max_projection_distance", "mean_projection_distance", "loaded_mex")},
                   normal_flux_norm=math.sqrt(sum(x*x for x in flux)) if flux is not None else None,
                   native_stages_seconds=stages,
                   measured_native_seconds_excluding_capture=sum(value for key, value in stages.items() if key != "diagnostic_capture"),
                   timing_caveat="The native call wall time includes diagnostic request/BRep writes. Stage timings isolate this cost. No serial original full-mesh completion was timed for this target.")


def preprojection_identity(directory, root):
    baseline = directory / "baseline-capture-hairy" / HAIRY_PREFIX
    candidate = directory / "capture-hairy" / HAIRY_PREFIX
    log = directory / "baseline-capture-hairy.log"
    stopped_as_expected = log.is_file() and "BASELINE_PREPROJECTION_CAPTURE_COMPLETE" in log.read_text()
    artifacts = []
    for extension in (".request.txt", ".map.txt", ".brep"):
        a = file_record(Path(str(baseline) + extension), root)
        b = file_record(Path(str(candidate) + extension), root)
        artifacts.append({"kind": extension, "baseline": a, "candidate": b,
                          "byte_equal": a.get("sha256") == b.get("sha256") if a["exists"] and b["exists"] else None})
    if not stopped_as_expected:
        return section("pending", artifacts=artifacts, reason="Original pipeline has not reached its intentional stop after preprojection capture")
    passed = all(row["byte_equal"] is True for row in artifacts)
    metrics = request_metrics(Path(str(candidate) + ".request.txt"), Path(str(candidate) + ".map.txt"))
    return section("passed" if passed else "failed", artifacts=artifacts, capture_metrics=metrics,
                   baseline_intentionally_stopped_after_capture=True,
                   baseline_is_not_a_completed_full_mesh=True,
                   baseline_stop_evidence=file_record(log, root))


def replay_output(path, expected_nodes, expected_digest):
    if not path.is_file():
        raise ValueError(f"Completed replay output is missing: {path}")
    with path.open("rb") as handle:
        header = handle.readline().split()
        actual_rows = sum(1 for _ in handle)
    if len(header) != 3 or header[0] != b"OK" or int(header[1]) != expected_nodes or actual_rows != expected_nodes:
        raise ValueError(f"Completed replay output row count or status is invalid: {path}")
    digest = sha256(path)
    if digest != expected_digest:
        raise ValueError(f"Completed replay output no longer matches its recorded hash: {path}")
    return digest, int(header[2])


def full_replay(directory, root):
    path = directory / "full-replay" / "report.json"
    data = read_json(path)
    batches = sorted(rows(data.get("batches")), key=lambda row: row["index"])
    errors = rows(data.get("errors"))
    source_nodes = int(data["source_nodes"])
    coverage_end = 0
    verified = 0
    total_fallbacks = 0
    for batch in batches:
        count = int(batch["nodes"])
        if batch["end"] - batch["begin"] != count or count <= 0:
            raise ValueError("Replay batch range and node count differ")
        if batch.get("exact") is not True or batch.get("baseline_sha256") != batch.get("candidate_sha256"):
            errors.append({"index": batch["index"], "error": "Baseline/candidate replay bytes differ"})
            continue
        folder = directory / "full-replay" / f"batch-{batch['index']:04d}"
        a, fallback = replay_output(folder / "baseline.txt", count, batch["baseline_sha256"])
        b, candidate_fallback = replay_output(folder / "candidate.txt", count, batch["candidate_sha256"])
        if a != b or fallback != candidate_fallback or fallback != batch["fallbacks"]:
            raise ValueError("Replay output bytes or fallback counts differ")
        verified += count
        total_fallbacks += fallback
    if data.get("completed") is True:
        for index, batch in enumerate(batches):
            if batch["index"] != index or batch["begin"] != coverage_end:
                raise ValueError("Completed replay does not cover consecutive, nonoverlapping source ranges")
            coverage_end = batch["end"]
        if coverage_end != source_nodes:
            raise ValueError("Completed replay coverage does not include every source node")
    request = Path(data["request"])
    mapping = directory / "capture-hairy" / (HAIRY_PREFIX + ".map.txt")
    source_matches = (request.is_file() and sha256(request) == data["request_sha256"] and
                      mapping.is_file() and sha256(mapping) == data["mapping_sha256"])
    if not source_matches:
        errors.append({"error": "Replay request/map no longer match their recorded source hashes"})
    completed = data.get("completed") is True
    passed = (completed and data.get("passed") is True and not errors and verified == source_nodes and
              data.get("verified_nodes") == verified)
    status = "failed" if errors or (completed and not passed) else "passed" if passed else "pending"
    timings = {}
    for variant in ("baseline", "candidate"):
        values = [row.get(variant + "_wall_seconds") for row in batches]
        known = [value for value in values if value is not None]
        if any(not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0
               for value in known):
            raise ValueError("Replay contains an invalid process duration")
        timings[variant] = {"timed_batches": len(known), "untimed_batches": len(values)-len(known),
                           "sum_timed_process_wall_seconds": sum(known) if known else None}
    return section(status, evidence=file_record(path, root), source_nodes=source_nodes,
                   verified_nodes=verified, reported_verified_nodes=data.get("verified_nodes", 0),
                   completed_batches=len(batches), expected_batches=math.ceil(source_nodes/data["batch_size"]),
                   fallback_count=total_fallbacks, jobs=data["jobs"], batch_size=data["batch_size"],
                   outputs_rehashed=True, errors=errors, request_sha256=data["request_sha256"],
                   mapping_sha256=data["mapping_sha256"], elapsed_wall_seconds=data.get("elapsed_wall_seconds"),
                   process_timings=timings,
                   timing_caveat="Batches run in independent concurrent processes, each reimporting geometry and rebuilding grids. These process or aggregate wall times are not serial full-mesh baseline/candidate timings and do not establish a full-pipeline speedup.")


def full_mesh_replay_comparison(directory, root):
    path = directory / "full-replay-matlab-comparison.json"
    data = read_json(path)
    replay = read_json(directory / "full-replay/report.json")
    request = directory / "capture-hairy" / (HAIRY_PREFIX + ".request.txt")
    mapping = directory / "capture-hairy" / (HAIRY_PREFIX + ".map.txt")
    metrics = request_metrics(request, mapping)
    runs = rows(data.get("runs"))
    candidate_file = directory / "capture-hairy/capture.mat"
    replay_file = directory / "full-replay/report.json"
    flags = ("xyz_bitwise_exact", "distance_bitwise_exact", "face_order_bitwise_exact",
             "fallback_count_exact", "dimensions_exact", "passed")
    passed = (data.get("passed") is True and bool(runs) and
              data.get("replay_xyz_distance_uv_byte_exact") is True and
              data.get("nodes") == metrics["source_nodes"] and
              data.get("batches") == len(rows(replay.get("batches"))) and
              data.get("fallbacks") == 0 and
              data.get("request_sha256") == sha256(request) and
              Path(data.get("candidate_file", "")).resolve() == candidate_file.resolve() and
              Path(data.get("replay_report_file", "")).resolve() == replay_file.resolve() and
              data.get("candidate_sha256") == sha256(candidate_file) and
              data.get("replay_report_sha256") == sha256(replay_file) and
              all(all(run.get(flag) is True for flag in flags) for run in runs))
    return section("passed" if passed else "failed", evidence=file_record(path, root),
                   candidate_mesh_evidence=file_record(candidate_file, root),
                   replay_evidence=file_record(replay_file, root),
                   nodes=data.get("nodes"), batches=data.get("batches"), fallbacks=data.get("fallbacks"),
                   request_sha256=data.get("request_sha256"), runs=runs,
                   replay_xyz_distance_uv_byte_exact=data.get("replay_xyz_distance_uv_byte_exact"),
                   projected_uv_in_matlab_mesh=data.get("projected_uv_in_matlab_mesh"),
                   scope="The completed candidate mesh XYZ, projection distances, face order and fallback count are checked against all independent original projector outputs; replay UV is byte exact, while the mesh stores scaffold UV rather than projected UV.")


def microbenchmark(directory, root):
    path = directory / "micro-hairy.jsonl"
    if not path.is_file():
        raise Pending("Hairy-torus microbenchmark rows are missing")
    data = [json.loads(line) for line in path.read_text().splitlines() if line.startswith("{")]
    timed = [row for row in data if "oracle_digest" in row]
    if not timed:
        raise Pending("No completed hairy-torus microbenchmark comparisons")
    passed = all(row.get("exact") is True and row.get("oracle_digest") == row.get("candidate_digest") for row in timed)
    return section("passed" if passed else "failed", evidence=file_record(path, root), rows=timed,
                   scope="Completed printed comparisons only; synthetic seam/narrow-span and mixed/shuffled/repeated batches, not the full mesh.")


def actual_rv_timing(directory, root):
    path = directory / "rv-sample/benchmark.json"
    data = read_json(path)
    runs = rows(data.get("runs"))
    order = ["baseline", "candidate", "candidate", "baseline", "candidate", "baseline", "baseline", "candidate"]
    if len(runs) < 8 or "summary" not in data:
        raise Pending(f"Actual RV timing has completed {len(runs)} of eight serial helper processes")
    sample = data["sample"]
    request = Path(sample["request_file"])
    mapping = Path(sample["map_file"])
    metrics = request_metrics(request, mapping)
    passed = (data.get("all_response_bytes_equal") is True and data.get("order") == order and
              [run.get("variant") for run in runs] == order and metrics["source_nodes"] == 512 and
              sample.get("sample_nodes") == 512 and sample.get("canonical_face_index") == 0 and
              sha256(request) == sample.get("request_sha256") and sha256(mapping) == sample.get("map_sha256"))
    original = Path(sample["original_request"])
    passed = passed and original.is_file() and sha256(original) == sample.get("original_request_sha256")
    for group in ("helpers", "cores"):
        for record in sample[group].values():
            source = Path(record["path"])
            passed = passed and source.is_file() and sha256(source) == record["sha256"]
    digests = []
    compact_runs = []
    for run in runs:
        digest, fallback = replay_output(Path(run["response_file"]), 512, run["response_sha256"])
        digests.append(digest)
        stdout = Path(run["stdout_file"])
        passed = passed and (run.get("strict_response_bytes_equal") is True and run.get("exit_status") == 0 and
                             fallback == 0 and stdout.is_file() and sha256(stdout) == run["stdout_sha256"])
        compact_runs.append({key: run.get(key) for key in (
            "index", "variant", "wall_seconds", "user_seconds", "system_seconds", "peak_rss_bytes",
            "response_sha256", "strict_response_bytes_equal", "background")})
    passed = passed and len(set(digests)) == 1
    return section("passed" if passed else "failed", evidence=file_record(path, root),
                   request_evidence=file_record(request, root), map_evidence=file_record(mapping, root),
                   sample_nodes=512, original_source_nodes=sample.get("source_total_nodes"),
                   original_main_spline_nodes=sample.get("source_main_face_nodes"), selection=sample.get("selection"),
                   runs=compact_runs, summary=data.get("summary"), paired_speedups=data.get("paired_speedups"),
                   sample_median_wall_speedup=data.get("sample_median_wall_speedup"),
                   sample_median_user_speedup=data.get("sample_median_user_speedup"),
                   candidate_to_baseline_peak_rss_ratio=data.get("candidate_to_baseline_peak_rss_ratio"),
                   background_note=data.get("background_note"), timing_caveat=data.get("interpretation"))


def validate_timing_report(data, root, expected_cases):
    cases = rows(data.get("cases"))
    valid = (data.get("exact_payloads_passed") is True and data.get("measurements_completed") is True and
             data.get("process_order") == ["baseline-a", "candidate-a", "candidate-b", "baseline-b"] and
             data.get("warm_samples_per_variant_per_fixture") == 12 and
             set(data.get("expected_cases", [])) == expected_cases and
             {case["name"] for case in cases} == expected_cases and len(cases) == len(expected_cases))
    valid = valid and all(case["baseline"]["count"] == 12 and case["candidate"]["count"] == 12 for case in cases)
    actual_flags = {case["name"] for case in cases if case.get("flag_for_quiet_rerun") is True}
    valid = valid and actual_flags == set(data.get("quiet_rerun_required_cases", []))
    for record in rows(data.get("capture_evidence")) + rows(data.get("comparison_evidence")):
        source = root / record["path"]
        valid = valid and source.is_file() and sha256(source) == record.get("sha256")
    valid = valid and len(rows(data.get("capture_evidence"))) == 8 and len(rows(data.get("comparison_evidence"))) == 2
    return valid


def repeated_timing(directory, root):
    path = directory / "repeated-timing-summary.json"
    data = read_json(path)
    status = data.get("status")
    if status not in ("passed", "pending", "failed"):
        raise ValueError("Repeated timing report has an invalid status")
    cases = rows(data.get("cases"))
    flags = set(data.get("quiet_rerun_required_cases", []))
    quiet = None
    quiet_path = directory / "quiet-timing/repeated-timing-summary.json"
    if data.get("measurements_completed") is True:
        if not validate_timing_report(data, root, SMALL_CASES):
            status = "failed"
        elif flags:
            if quiet_path.is_file():
                quiet = read_json(quiet_path)
                if quiet.get("measurements_completed") is not True:
                    status = "failed" if quiet.get("status") == "failed" else "pending"
                else:
                    quiet_valid = validate_timing_report(quiet, root, flags)
                    quiet_passed = (quiet_valid and quiet.get("status") == "passed" and quiet.get("completed") is True and
                                    quiet.get("passed") is True and not quiet.get("quiet_rerun_required_cases"))
                    status = "passed" if quiet_passed else "failed"
            else:
                status = "pending"
        elif not (status == "passed" and data.get("completed") is True and data.get("passed") is True):
            status = "failed"
    return section(status, evidence=file_record(path, root),
                   warm_samples_per_variant_per_fixture=data.get("warm_samples_per_variant_per_fixture"),
                   exact_payloads_passed=data.get("exact_payloads_passed"), cases=cases,
                   initial_quiet_rerun_required_cases=sorted(flags),
                   quiet_rerun_evidence=file_record(quiet_path, root), quiet_rerun=quiet,
                   quiet_rerun_resolved=bool(quiet is not None and status == "passed"),
                   background_load=data.get("background_load"), timing_caveat=data.get("timing_caveat"),
                   reason=data.get("reason"))


def build_report(directory, root):
    operations = {
        "native_ctest": native_tests,
        "native_projection_matrix": native_projection_matrix,
        "small_matlab_differential": small_comparison,
        "legacy_extraction_parity": extraction_parity,
        "candidate_smoke": candidate_smoke,
        "hairy_full_projection_replay": full_replay,
        "hairy_actual_mesh_projection_comparison": full_mesh_replay_comparison,
        "hairy_candidate_full_mesh": candidate_hairy,
        "hairy_preprojection_identity": preprojection_identity,
        "repeated_small_mesh_performance": repeated_timing,
    }
    checks = {name: guarded(lambda operation=operation: operation(directory, root))
              for name, operation in operations.items()}
    passed = all(check["status"] == "passed" for check in checks.values())
    status = "passed" if passed else "failed" if any(check["status"] == "failed" for check in checks.values()) else "pending"
    return {
        "schema_version": 1, "component": "SurfSmooth3D STEP node projection",
        "baseline_commit": BASELINE_COMMIT,
        "implementation_files": [file_record(root / name, root) for name in
            ("step_mesher/cpp/src/occ_context.cpp", "step_mesher/CMakeLists.txt")],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status, "completed": passed, "passed": passed,
        "required_checks": checks,
        "optional_evidence": {
            "hairy_microbenchmark": guarded(lambda: microbenchmark(directory, root)),
            "actual_rv_serial_timing": guarded(lambda: actual_rv_timing(directory, root)),
        },
        "scope": "Reuse the same initialized bounded OCC spline/Bezier projector within consecutive source-face runs; preserve original global nearest-point algorithm, tolerance, order, UV, distances, and fallbacks.",
        "limits": [
            "The isolated baseline benchmark stops after preprojection capture; its serialized geometry, face map and all input nodes must match the candidate before replay is accepted.",
            "Full target projection correctness is checked against an independent baseline executable for every captured node; candidate end-to-end mesh correctness also includes the recorded RV-derived and area/flux diagnostics.",
            "Parallel replay timing is not a serial original full-mesh runtime measurement. Fixed-input projection benchmarks and repeated small-case timings provide the direct speed evidence.",
            "The initial small-case captures contain single measurements. The separate ABBA gate uses repeated warm samples and a 5% plus 1ms median threshold for these fixtures; it does not prove a universal absence of performance regressions.",
        ],
    }


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-directory", type=Path, default=root / "build/validation/fast")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-complete", action="store_true", help="Exit nonzero unless every required check passes")
    args = parser.parse_args()
    directory = args.evidence_directory.resolve()
    output = args.output or directory / "summary.json"
    report = build_report(directory, root)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)
    print(json.dumps({"output": str(output), "status": report["status"], "completed": report["completed"],
                      "checks": {name: check["status"] for name, check in report["required_checks"].items()}}))
    if args.require_complete and not report["completed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
