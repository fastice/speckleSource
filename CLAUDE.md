# CLAUDE.md — speckleSource

Supplements the root `/Users/ian/progs/GIT64/CLAUDE.md`. SAR speckle-tracking and offset-culling
programs. Each subdirectory builds one binary; there is no shared `common/` here (unlike
`mosaicSource/`) — each program is largely self-contained.

## Programs

Each has its own `make <target>` rule in `speckleSource/Makefile` (run from `speckleSource/`).

| Binary | Source dir | Purpose | Doc |
|---|---|---|---|
| `strack` | `Strack/` | Main speckle tracker (amplitude + complex tracking) | `Documents/strack.md` |
| `strackw` | `Strackw/` | Fast-offsets speckle tracker, optional pre-smoothing for decorrelated speckle | `Documents/strackw.md` |
| `cullst` | `Cullst/` | Culls offsets based on neighborhood statistics | — |
| `cullls` | `Cullls/` | Culls Landsat offsets | — |

`make all` builds `strack strackw cullst cullls`.

## Notes

- **OpenMP** — `strack`, `strackw`, and `cullst` all support `-ompThreads N` (strack/cullst
  default 4, strackw default 2). Inner j-loop (range direction) is parallelized
  (`schedule(dynamic, 4)`); outer i-loop (azimuth) stays serial. Per-thread work arrays are
  module-level globals marked `#pragma omp threadprivate(...)`, malloc'd/freed via
  `mallocPerThreadArrays[W]()`/`freePerThreadArrays[W]()` inside `#pragma omp parallel` blocks.
  Image line buffers (`imageBuf1`/`imageBuf2`) are shared singletons — safe because all j for a
  given i share the same azimuth position.
  - **a2 pre-load pattern**: before each j-loop, `imageBuf2` is anchored at the minimum valid a2
    across the row by probing j=0/center/last via `findImage2Pos`, discarding out-of-range probes
    and clamping to `[0, nSlpA2-1]`. This avoids buffer ping-pong when threads span a wide a2
    range. `updateSLCBuffer`/`updateSLCBufferVRT` have an early-return guard for `a1` out of
    `[0, nSlpA)`. Any reload still occurring inside the parallel section is protected by
    `#pragma omp critical(imagebuf_reload[_w])`.
  - **cullst**: validated — 100% mask agreement vs serial on a 2844x114 truth grid; smoothed-value
    differences up to ~0.03-0.05 px are expected FP-summation-order noise from `cullSmooth.c`.
  - **strack/strackw**: not yet validated against serial output. Known benign race: `mrqminMod`
    has static locals that race when `gaussFlag==TRUE`; `gaussFlag` is always FALSE for glacier
    use, so low priority (do not "fix" without checking `testStrack.py` results first).
- **`Cullst/` vs `Cullls/`** — despite identically-named files (`cullSTData.c`, `loadCullData.c`,
  `cullSmooth`/`cullStats` variants), the two directories' copies have **diverged** — they are not
  interchangeable. Edit each in place; don't assume a fix in one applies to the other.
- **`BUGFIXES_speckleSource_2026-05-24.md`** — dated record of fixes from a past automated review
  covering strack, strackw, cullst, cullls. Historical, not living docs.
- Validate strack/strackw/cullst changes with `tests/sTrackTest/testStrack.py` (compares
  raw/culled/interpolated offsets against truth files); run with `-ompThreads 1` vs `N` to check
  parallel correctness.
