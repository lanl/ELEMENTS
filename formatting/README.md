# ELEMENTS code formatting

ELEMENTS uses MATAR's own formatting pipeline directly — the style config
and macro-reflow post-processor both live in the `matar/` submodule
(`matar/.clang-format`, `matar/formatting/matar-format.py`) and are not
duplicated here. `formatting/format.sh` is a thin wrapper that points
clang-format at MATAR's config and scopes the two-stage pipeline to
ELEMENTS' own source.

## Usage

```bash
formatting/format.sh            # format src/, examples/, tests/ in place
formatting/format.sh --check    # report non-conforming files; exits 1 if any
```

`--check` never modifies the working tree — it runs the pipeline against
scratch copies and diffs the result.

## What gets formatted

`src/`, `examples/`, and `tests/` (the only trees under `.h`/`.cpp`
ELEMENTS owns). Excluded: `matar/` (the submodule — MATAR formats its own
code), `legacy/` (unreferenced by the active build; reformatting it is pure
diff noise), and any `build*/` directory.

## The pipeline

See `matar/formatting/README.md` for the full explanation. In short:

1. **clang-format**, styled via `matar/.clang-format` (Google base, 4-space
   indents, 150-column limit, and macro-awareness for `MATAR_*` statement
   macros and `KOKKOS_*` attribute macros).
2. **`matar/formatting/matar-format.py`**, a post-processor that rewrites
   MATAR parallel-macro calls (`FOR_ALL`, `FOR_ALL_CLASS`, `FOR_REDUCE_*`,
   etc.) into the canonical multi-line layout clang-format cannot produce on
   its own — it has no notion of a macro's numeric index-triple arguments
   and would otherwise pack or split them arbitrarily.

Always run both stages as a pair (which `format.sh` does automatically):
clang-format alone will re-pack the macro headers, and the reflow pass
restores them. Applied together the result is stable — repeated runs are
byte-identical.
