#!/usr/bin/env python3
"""
The Cortex-M7 gate: build the benchmark, run it under QEMU, and check that

  1. NMPC's own code links no double-precision support routine (__aeabi_d*,
     __aeabi_*2d) in this single-precision build, with the default solver or
     with per-stage Jacobians, and
  2. each instruction count is within a tolerance of the accepted baseline
     for the same compiler version.

Counts move by up to 15 % between compiler versions, so a baseline is kept
per compiler. A compiler with no baseline passes with a notice. After a
deliberate change of cost, run with --accept to record the new counts.

Usage: check.py <eigen include directory> [--tolerance 0.03] [--accept]
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
BASELINE = HERE / "accepted.json"
DOUBLE = re.compile(r"__aeabi_(d[a-z0-9]+|[a-z0-9]+2d)\b")
COUNT = re.compile(r"^(.+?): (\d+) instructions, stack (\d+) bytes$")


def compiler_version():
    first = subprocess.run(["arm-none-eabi-gcc", "--version"], check=True, capture_output=True,
                           text=True).stdout.splitlines()[0]
    return first.split(")")[-1].strip().split()[0]


def parse(text):
    counts, stacks = {}, {}
    for line in text.splitlines():
        match = COUNT.match(line.strip())
        if match:
            counts[match.group(1)] = int(match.group(2))
            stacks[match.group(1)] = int(match.group(3))
    return counts, stacks


def callers_of_double_routines(elf):
    """Functions in the image that call a double-precision routine, for information."""
    listing = subprocess.run(["arm-none-eabi-objdump", "-d", str(elf)], check=True, capture_output=True,
                             text=True).stdout
    callers, current = set(), None
    for line in listing.splitlines():
        header = re.match(r"^[0-9a-f]+ <(.+)>:$", line)
        if header:
            current = header.group(1)
        elif current and re.search(r"\bbl?\b.*<" + DOUBLE.pattern, line):
            callers.add(current)
    return sorted(callers)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("eigen")
    parser.add_argument("--tolerance", type=float, default=0.03)
    parser.add_argument("--accept", action="store_true")
    parser.add_argument("--output", default=str(HERE / "out"))
    args = parser.parse_args()

    run = subprocess.run([str(HERE / "build.sh"), args.eigen, args.output], capture_output=True, text=True)
    sys.stdout.write(run.stdout)
    if run.returncode != 0:
        sys.stderr.write(run.stderr)
        print("FAIL: the benchmark did not build or run")
        return 1
    if run.stdout.count("all solves usable") != 2:
        print("FAIL: a solve was not usable")
        return 1

    output = pathlib.Path(args.output)
    undefined = "".join(subprocess.run(["arm-none-eabi-nm", "-u", str(output / name)], check=True,
                                       capture_output=True, text=True).stdout
                        for name in ("main.o", "main_per_stage.o"))
    own = sorted(set(DOUBLE.findall(undefined)))
    if own:
        print("FAIL: NMPC's code calls double-precision routines: " + ", ".join("__aeabi_" + s for s in own))
        return 1
    print("No double-precision routine is called from NMPC's code.")
    others = callers_of_double_routines(output / "benchmark.elf")
    if others:
        print("For information, the C library's " + ", ".join(others) + " call double-precision routines.")

    counts, _ = parse(run.stdout)
    compiler = compiler_version()
    baselines = json.loads(BASELINE.read_text()) if BASELINE.exists() else {}
    if args.accept:
        baselines[compiler] = counts
        BASELINE.write_text(json.dumps(baselines, indent=2, sort_keys=True) + "\n")
        print(f"Accepted the counts for arm-none-eabi-gcc {compiler}.")
        return 0
    baseline = baselines.get(compiler)
    if baseline is None:
        print(f"No baseline for arm-none-eabi-gcc {compiler}: counts not gated. Run with --accept to record them.")
        return 0
    failed = False
    for name, accepted in sorted(baseline.items()):
        if name not in counts:
            print(f"FAIL: {name} was not measured")
            failed = True
            continue
        change = counts[name] / accepted - 1.0
        verdict = "ok" if abs(change) <= args.tolerance else "FAIL"
        failed = failed or verdict == "FAIL"
        print(f"{verdict:4} {name}: {counts[name]} against {accepted} ({change:+.1%})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
