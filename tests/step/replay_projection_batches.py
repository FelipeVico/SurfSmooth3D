#!/usr/bin/env python3
"""Compare every captured projection against a separate baseline executable.

Splits an existing STEP_MESHER_OCC_PROJECT_V1 request into consecutive batches.
Each process imports its own geometry, so parallel validation shares no OCC
state. Byte comparison includes XYZ, UV, distances, and fallback counts.
Concurrent replay wall times are diagnostic, not serial end-to-end benchmarks.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def valid_response(path, count):
    if not path.is_file():
        return False
    rows = path.read_text().splitlines()
    if len(rows) != count + 1:
        return False
    try:
        status = rows[0].split()
        return (len(status) == 3 and status[0] == "OK" and
                int(status[1]) == count and 0 <= int(status[2]) <= count and
                all(len(row.split()) == 6 and
                    all(math.isfinite(float(value)) for value in row.split())
                    for row in rows[1:]))
    except ValueError:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "mapping", "baseline", "candidate", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=5000)
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--resume", action="store_true",
                        help="Reuse complete responses with identical request and executable provenance")
    args = parser.parse_args()
    if min(args.jobs, args.batch_size, args.timeout) <= 0:
        parser.error("jobs, batch-size and timeout must be positive")
    lines = args.request.read_text().splitlines(keepends=True)
    if lines[0].strip() != "STEP_MESHER_OCC_PROJECT_V1":
        parser.error("unrecognized request format")
    face_count = int(lines[2])
    header = lines[:face_count + 3]
    count = int(lines[face_count + 3])
    points = lines[face_count + 4:]
    if len(points) != count:
        parser.error("request point count does not match rows")
    args.output.mkdir(parents=True, exist_ok=True)
    mapping = args.mapping.resolve()
    reference = args.baseline.resolve()
    candidate = args.candidate.resolve()
    report = dict(request=str(args.request.resolve()), request_sha256=sha256(args.request),
                  mapping_sha256=sha256(mapping), baseline=str(reference),
                  candidate=str(candidate), source_nodes=count, jobs=args.jobs,
                  baseline_sha256=sha256(reference), candidate_sha256=sha256(candidate),
                  batch_size=args.batch_size, completed=False, passed=False, batches=[])
    report_path = args.output / "report.json"
    if args.resume and report_path.exists():
        previous = json.loads(report_path.read_text())
        for key in ("request_sha256", "mapping_sha256", "baseline", "candidate", "batch_size",
                    "baseline_sha256", "candidate_sha256"):
            if previous.get(key) != report[key]:
                parser.error(f"Cannot resume: {key} changed or provenance is missing")

    def checkpoint():
        temporary = report_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(report_path)

    def run_batch(index, begin, end):
        directory = args.output / f"batch-{index:04d}"
        directory.mkdir(exist_ok=True)
        request = directory / "request.txt"
        request_text = "".join(header) + f"{end-begin}\n" + "".join(points[begin:end])
        if args.resume and request.exists() and request.read_text() != request_text:
            raise RuntimeError(f"Cannot resume: request differs in batch {index}")
        request.write_text(request_text)
        result = dict(index=index, begin=begin, end=end, nodes=end-begin)
        for label, executable in (("baseline", reference), ("candidate", candidate)):
            output = directory / f"{label}.txt"
            if args.resume and valid_response(output, end-begin):
                result[label + "_reused"] = True
                result[label + "_wall_seconds"] = None
                result[label + "_sha256"] = sha256(output)
                continue
            command = [str(executable), "project-nodes", str(request.resolve()),
                       str(mapping), str(output.resolve())]
            started = time.perf_counter()
            with (directory / f"{label}.stdout").open("w") as stdout, \
                 (directory / f"{label}.stderr").open("w") as stderr:
                subprocess.run(command, stdout=stdout, stderr=stderr,
                               check=True, timeout=args.timeout)
            result[label + "_wall_seconds"] = time.perf_counter() - started
            if not valid_response(output, end-begin):
                raise RuntimeError(f"Incomplete or invalid {label} response in batch {index}")
            result[label + "_sha256"] = sha256(output)
        baseline_bytes = (directory / "baseline.txt").read_bytes()
        candidate_bytes = (directory / "candidate.txt").read_bytes()
        result["exact"] = baseline_bytes == candidate_bytes
        status = baseline_bytes.splitlines()[0].split()
        if status[:1] != [b"OK"] or int(status[1]) != end-begin:
            raise RuntimeError(f"Unexpected response count in batch {index}")
        result["fallbacks"] = int(status[2])
        return result

    started = time.perf_counter()
    checkpoint()
    errors = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_batch, i, begin, min(begin + args.batch_size, count)): i
                   for i, begin in enumerate(range(0, count, args.batch_size))}
        for future in as_completed(futures):
            try:
                result = future.result()
                report["batches"].append(result)
                report["batches"].sort(key=lambda value: value["index"])
                print(json.dumps(result), flush=True)
            except Exception as error:
                errors.append(dict(index=futures[future], error=str(error)))
                print(json.dumps(errors[-1]), flush=True)
            report["errors"] = errors
            report["verified_nodes"] = sum(row["nodes"] for row in report["batches"] if row["exact"])
            checkpoint()
    report["elapsed_wall_seconds"] = time.perf_counter() - started
    report["completed"] = True
    report["passed"] = not errors and report["verified_nodes"] == count
    checkpoint()
    if not report["passed"]:
        raise SystemExit("Projection replay differs or did not complete; inspect report.json")


if __name__ == "__main__":
    main()
