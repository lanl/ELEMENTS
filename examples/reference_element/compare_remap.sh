#!/bin/bash
#
# Compare remap_dg_lumped_test (baseline) with remap_dg_lumped_optimized.
#
# Reports median wall time for each and diffs their outputs.  Works for any
# Kokkos backend; on a CUDA build just point it at the CUDA build directory:
#
#     ./compare_remap.sh -b ../../build-cuda
#
# Notes
#   * The baseline hardcodes its problem (8x8x1, max_time 0.5, graphics_dt 0.1)
#     and parses no arguments, so that is the only case both binaries can run.
#     The optimized binary is handed exactly those values.  To compare a larger,
#     more GPU-meaningful mesh you have to edit the constants near the top of
#     main() in src/remap_dg_lumped_test.cpp.
#   * The baseline applies limit_corner_field every Runge-Kutta stage and the
#     optimized version does not, so part of any speedup is work the optimized
#     version simply does not do, and the field values legitimately differ.
#   * The baseline has no internal timer, so both runs are timed externally
#     around the whole process (this includes Kokkos/MPI startup).
#   * If remap_dg_lumped_test_cuda6 is present it is also run, with REMAP_OPT=6.
#     That binary is the optimization reference the clean version was extracted
#     from, so its output should match remap_dg_lumped_optimized *exactly*.
#
set -u

BUILD_DIR=""
REPS=3
KEEP=0

usage() {
    cat <<EOF
usage: $0 [-b build_dir] [-r reps] [-k]

  -b  build directory (default: auto-detect ../../build*, newest first)
  -r  repetitions per binary, median reported (default: $REPS)
  -k  keep the scratch run directories instead of deleting them
EOF
    exit 1
}

while getopts "b:r:kh" opt; do
    case $opt in
        b) BUILD_DIR=$OPTARG ;;
        r) REPS=$OPTARG ;;
        k) KEEP=1 ;;
        *) usage ;;
    esac
done

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

# ---- locate the binaries -------------------------------------------------
if [ -z "$BUILD_DIR" ]; then
    for d in "$REPO_ROOT"/build*/; do
        [ -x "$d/examples/reference_element/remap_dg_lumped_optimized" ] && BUILD_DIR=$d
    done
fi
if [ -z "$BUILD_DIR" ]; then
    echo "ERROR: no build directory found; pass one with -b" >&2
    exit 1
fi
BIN="$(cd "$BUILD_DIR" && pwd)/examples/reference_element"

BASE="$BIN/remap_dg_lumped_test"
OPT="$BIN/remap_dg_lumped_optimized"
REF="$BIN/remap_dg_lumped_test_cuda6"

for exe in "$BASE" "$OPT"; do
    if [ ! -x "$exe" ]; then
        echo "ERROR: missing $exe" >&2
        echo "       build it first, e.g. cmake --build $BUILD_DIR -j" >&2
        exit 1
    fi
done

# the baseline's built-in problem; the optimized binary is matched to it
CASE_ARGS="8 8 1 0.5 0.1"

echo "build     : $BIN"
echo "problem   : $CASE_ARGS  (baseline's built-in case)"
echo "reps      : $REPS (median reported)"
[ -n "${OMP_NUM_THREADS:-}" ] && echo "OMP_THREADS: $OMP_NUM_THREADS"
echo ""

# ---- helpers -------------------------------------------------------------
median() {
    printf '%s\n' "$@" | sort -g |
      awk '{a[NR]=$1} END{print (NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2}'
}

# run a binary in its own directory, echo elapsed seconds
timeit() {
    local exe=$1 dir=$2; shift 2
    local t0 t1
    mkdir -p "$dir" && cd "$dir" || exit 1
    t0=$(date +%s.%N)
    "$exe" "$@" > stdout.txt 2>&1
    local rc=$?
    t1=$(date +%s.%N)
    cd - > /dev/null || exit 1
    if [ $rc -ne 0 ]; then
        echo "ERROR: $exe exited $rc; see $dir/stdout.txt" >&2
        echo "NaN"
        return
    fi
    echo "$t1 - $t0" | bc
}

WORK=$(mktemp -d -t remap_cmp_XXXXXX)
cleanup() { [ "$KEEP" -eq 0 ] && rm -rf "$WORK" || echo "kept: $WORK"; }
trap cleanup EXIT

# ---- timed runs ----------------------------------------------------------
bt=(); ot=(); rt=()
for r in $(seq "$REPS"); do
    bt+=( "$(timeit "$BASE" "$WORK/base")" )
    ot+=( "$(timeit "$OPT"  "$WORK/opt" $CASE_ARGS)" )
    [ -x "$REF" ] && rt+=( "$(timeit "$REF" "$WORK/ref" $CASE_ARGS)" )
done

b=$(median "${bt[@]}")
o=$(median "${ot[@]}")

printf "%-34s %s\n" "remap_dg_lumped_test (baseline)" "$(printf '%.3f s' "$b")"
printf "%-34s %s\n" "remap_dg_lumped_optimized"       "$(printf '%.3f s' "$o")"
if [ ${#rt[@]} -gt 0 ]; then
    rmed=$(median "${rt[@]}")
    printf "%-34s %s\n" "remap_dg_lumped_test_cuda6 (ref)" "$(printf '%.3f s' "$rmed")"
fi
echo ""
printf "speedup (baseline / optimized): %sx\n" "$(echo "scale=2; $b/$o" | bc)"
echo "  (includes the limiter the baseline runs and the optimized version omits)"
echo ""

# ---- correctness ---------------------------------------------------------
echo "---- output comparison ----"

# baseline vs optimized: expected to differ, the baseline limits the field
if [ -f "$WORK/base/ErrorNorms.txt" ] && [ -f "$WORK/opt/ErrorNorms.txt" ]; then
    echo "baseline vs optimized error norms (differ by design - limiter):"
    paste <(tail -n +2 "$WORK/base/ErrorNorms.txt") \
          <(tail -n +2 "$WORK/opt/ErrorNorms.txt") |
      awk 'BEGIN{printf "  %-10s %-12s %-12s %-12s %-12s\n","time","L1_base","L1_opt","L2_base","L2_opt"}
           {printf "  %-10s %-12s %-12s %-12s %-12s\n",$1,$2,$5,$3,$6}'
fi
echo ""

# optimized vs the reference it was extracted from: must be identical
if [ -x "$REF" ]; then
    same=1
    diff -q "$WORK/ref/ErrorNorms.txt" "$WORK/opt/ErrorNorms.txt" > /dev/null 2>&1 || same=0
    for f in "$WORK/ref"/output_time_*.vtu; do
        [ -e "$f" ] || continue
        diff -q "$f" "$WORK/opt/$(basename "$f")" > /dev/null 2>&1 || same=0
    done
    if [ "$same" -eq 1 ]; then
        echo "optimized vs cuda6 reference: IDENTICAL (norms + all VTU files)"
    else
        echo "optimized vs cuda6 reference: *** DIFFERENCES FOUND ***"
        echo "  rerun with -k and diff the directories under the kept path"
    fi
else
    echo "optimized vs cuda6 reference: skipped (remap_dg_lumped_test_cuda6 not built)"
fi

# mass conservation is the strongest self-check either binary provides
echo ""
echo "---- final mass ----"
grep -h "Domain mass at time" "$WORK/base/stdout.txt" 2>/dev/null | sed 's/^/  baseline : /'
grep -h "Domain mass at time" "$WORK/opt/stdout.txt"  2>/dev/null | sed 's/^/  optimized: /'
