#!/bin/sh
# Cross-compile the NMPC Cortex-M7 benchmark and run it under QEMU.
# Usage: build.sh <eigen include directory> [output directory]
# Environment: OPT sets the optimisation level (default -O2).
set -e
EIGEN=${1:?Give the Eigen include directory, for example build/src/eigen3}
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${2:-$HERE/out}
OPT=${OPT:--O2}
ROOT=$HERE/../..
mkdir -p "$OUT"
FLAGS="-mcpu=cortex-m7 -mthumb -mfloat-abi=hard -mfpu=fpv5-d16 $OPT -ffunction-sections -fdata-sections -fno-math-errno"
arm-none-eabi-gcc $FLAGS -std=c11 -Wall -c "$HERE/startup.c" -o "$OUT/startup.o"
# The default solver, then the solver with per-stage Jacobians, each a program of its own.
for VARIANT in default per-stage; do
    SUFFIX=; DEFINES=
    if [ $VARIANT = per-stage ]; then SUFFIX=_per_stage; DEFINES=-DNMPC_BENCHMARK_PER_STAGE; fi
    arm-none-eabi-g++ $FLAGS -std=c++17 -fno-exceptions -fno-rtti -fno-threadsafe-statics -DNDEBUG -DEIGEN_NO_DEBUG \
        -DEIGEN_NO_MALLOC $DEFINES -isystem "$EIGEN" -I"$ROOT/include" -I"$ROOT/benchmark" \
        -Wall -Wextra -Wshadow -Wconversion -Wdouble-promotion -Werror -fstack-usage \
        -c "$HERE/main.cpp" -o "$OUT/main$SUFFIX.o"
    arm-none-eabi-g++ $FLAGS -nostartfiles -Wl,--no-warn-rwx-segments -T "$HERE/link.ld" -Wl,--gc-sections -Wl,-Map="$OUT/benchmark$SUFFIX.map" \
        --specs=nosys.specs "$OUT/startup.o" "$OUT/main$SUFFIX.o" -lm -o "$OUT/benchmark$SUFFIX.elf"
    arm-none-eabi-size "$OUT/benchmark$SUFFIX.elf"
    qemu-system-arm -M mps2-an500 -icount shift=0 -semihosting -nographic -monitor none -serial none \
        -kernel "$OUT/benchmark$SUFFIX.elf" 2>&1 | tee "$OUT/result$SUFFIX.txt"
done
