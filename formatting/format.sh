#!/usr/bin/env bash
# Formats ELEMENTS' own source (src/, examples/, tests/) using MATAR's
# clang-format style and macro-reflow post-processor, both sourced live from
# the matar/ submodule -- ELEMENTS keeps no copy of either. See
# formatting/README.md and matar/formatting/README.md for details.
#
# Usage:
#   formatting/format.sh            # format in place
#   formatting/format.sh --check    # report files that would change; exits 1 if any would
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
MATAR_DIR="${REPO_ROOT}/matar"
CLANG_FORMAT_CONFIG="${MATAR_DIR}/.clang-format"
REFLOW_SCRIPT="${MATAR_DIR}/formatting/matar-format.py"

CHECK=0
case "${1:-}" in
    --check) CHECK=1 ;;
    "") ;;
    *)
        echo "usage: $(basename "$0") [--check]" >&2
        exit 1
        ;;
esac

if [[ ! -f "${CLANG_FORMAT_CONFIG}" || ! -f "${REFLOW_SCRIPT}" ]]; then
    echo "error: matar/.clang-format or matar/formatting/matar-format.py not found." >&2
    echo "  the matar submodule may not be checked out -- run:" >&2
    echo "    git submodule update --init --recursive" >&2
    exit 1
fi
command -v clang-format >/dev/null 2>&1 || { echo "error: clang-format not found on PATH" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "error: python3 not found on PATH" >&2; exit 1; }

mapfile -t FILES < <(find "${REPO_ROOT}/src" "${REPO_ROOT}/examples" "${REPO_ROOT}/tests" \
    -type f \( -name "*.h" -o -name "*.cpp" \) | sort)

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "no files to format"
    exit 0
fi

if [[ "${CHECK}" -eq 0 ]]; then
    clang-format -i --style="file:${CLANG_FORMAT_CONFIG}" "${FILES[@]}"
    python3 "${REFLOW_SCRIPT}" "${FILES[@]}"
    echo "formatted ${#FILES[@]} files"
    exit 0
fi

# --check: run the identical two-stage pipeline against scratch copies so the
# working tree is never touched, then diff each copy against the real file.
SCRATCH="$(mktemp -d)"
trap 'rm -rf "${SCRATCH}"' EXIT

for f in "${FILES[@]}"; do
    rel="${f#"${REPO_ROOT}"/}"
    copy="${SCRATCH}/${rel}"
    mkdir -p "$(dirname "${copy}")"
    cp "${f}" "${copy}"
done

mapfile -t COPIES < <(find "${SCRATCH}" -type f | sort)
clang-format -i --style="file:${CLANG_FORMAT_CONFIG}" "${COPIES[@]}"
python3 "${REFLOW_SCRIPT}" "${COPIES[@]}"

UNFORMATTED=()
for f in "${FILES[@]}"; do
    rel="${f#"${REPO_ROOT}"/}"
    if ! diff -q "${f}" "${SCRATCH}/${rel}" >/dev/null 2>&1; then
        UNFORMATTED+=("${rel}")
    fi
done

if [[ ${#UNFORMATTED[@]} -gt 0 ]]; then
    echo "the following files are not formatted (run formatting/format.sh to fix):"
    printf '  %s\n' "${UNFORMATTED[@]}"
    exit 1
fi

echo "all ${#FILES[@]} files are correctly formatted"
