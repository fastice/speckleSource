# Bug fixes applied to speckleSource 2026-05-24

Applied findings from automated code review of speckleSource (strack, strackw, cullst, cullls).
All changes compile clean with no warnings.

---

## Strack/speckleTrack.c

### 1. aZDefocus written to wrong indices in ampTrackFast()
**Problem:** `trackPar->aZDefocus[i][j] = -LARGEINT;` used `i` and `j` (local loop variables
from the FFT/peak-finding scope) instead of `iOut` and `jOut` (the output-pixel parameters).
The companion valid-peak assignment on the next lines correctly used `iOut`/`jOut`, so the
no-peak sentinel was written to wrong (potentially out-of-bounds) indices.

**Fix:** Changed `trackPar->aZDefocus[i][j] = -LARGEINT;` to
`trackPar->aZDefocus[iOut][jOut] = -LARGEINT;`

**To revert:** Change `[iOut][jOut]` back to `[i][j]` in the ampTrackFast sentinel assignment.

---

### 2. Hanning window uses wA/2 offset for range (j) dimension
**Problem:** Inside `makeHanning()`, the range column index was offset by `trackPar->wA / 2`:
```c
j1 = trackPar->wA / 2 + j;
j2 = trackPar->wA / 2 - j - 1;
```
The inner loop runs `j = 0 .. wR/2`, so the offset should use `wR/2` (range half-width).
When `wA != wR` (common in speckle tracking), this writes column indices relative to the
wrong centre, and can write past the end of the `hanning` array if `wR > wA`.

**Fix:** Changed both to `trackPar->wR / 2`.

**To revert:** Change both `wR / 2` back to `wA / 2`.

---

### 3. sqrt of amplitude-squared values has no negative guard
**Problem:** `im1[i][j].re = sqrt(im1[i][j].re);` — `re` holds amplitude-squared (power).
SAR power values should always be non-negative, but no guard existed. A negative value
(e.g., from floating-point underflow or a sentinel leaking through) would produce NaN silently.

**Fix:** Changed to:
```c
im1[i][j].re = (im1[i][j].re > 0.0) ? sqrt(im1[i][j].re) : 0.0;
im2[i][j].re = (im2[i][j].re > 0.0) ? sqrt(im2[i][j].re) : 0.0;
```

**To revert:** Remove the ternary guard and restore bare `sqrt(im1[i][j].re)` and `sqrt(im2[i][j].re)`.

---

### 4. Unchecked fread return value (LSB path)
**Problem:** `size_t rv = fread(...);` — return value declared but never used. Compiler
warning; also means short reads go undetected on the LSB code path.

**Fix:** Changed to `(void)fread(...)` to make the deliberate discard explicit.

**To revert:** Change `(void)fread(...)` back to `size_t rv = fread(...)`.

---

## Cullst/cullst.c

### Inverted condition in addSubtractSimOffsets()
**Problem:** The condition to apply simulated offsets was:
```c
if (cullPar->offA[j][i] < (-LARGEINT + 1) && cullPar->offSimA[j][i] < (-LARGEINT + 1))
```
This fires when **both** measured and simulated offsets are the no-data sentinel (invalid
pixels). The intent is the opposite: apply the sim offset when both are valid (greater than
the sentinel). The inverted logic meant sim offsets were never applied to valid pixels.

**Fix:** Changed both `<` to `>`:
```c
if (cullPar->offA[j][i] > (-LARGEINT + 1) && cullPar->offSimA[j][i] > (-LARGEINT + 1))
```

**To revert:** Change both `>` back to `<`.

---

## Cullst/loadCullData.c

### Uninitialised lineCount in three functions
**Problem:** `int32_t lineCount, eod, ...;` — `lineCount` was uninitialised in:
- `loadCullData()` (used immediately as the line counter passed to `getDataString`)
- `loadSimData()` (same)
- `loadCullMask()` (same)

All other callers of `getDataString` in this codebase initialise `lineCount = 0`. The
uninitialized value meant error messages would print garbage line numbers (no runtime crash
because the return value is reassigned).

**Fix:** Changed all three declarations to `int32_t lineCount = 0, eod, ...;`

**To revert:** Remove the `= 0` initialisers in the three function declarations.

---

## Cullls/loadCullData.c

### Uninitialised lineCount
**Problem:** Same as Cullst/loadCullData.c — `int32_t lineCount, eod;` passed uninitialised
to `getDataString`.

**Fix:** Changed to `int32_t lineCount = 0, eod;`

**To revert:** Remove the `= 0` initialiser.

---

## Strackw/corrTrackFast.c

### 1. Non-atomic timeFFT update causes data race and double-counting
**Problem:** After `ampMatchEdge(...)`:
```c
timeFFT += (double)((clock()-time1))/(CLOCKS_PER_SEC);   // non-atomic, line ~563
```
This line was neither protected by `#pragma omp atomic` nor by the compound-statement
pattern used for all other timing stats in the file. Immediately afterwards, the same
variable was updated again (atomically) including time for both `ampMatchEdge` and
`getPeakCorr`. The double update was both a data race and a double-counting bug.

**Fix:** Removed the non-atomic `timeFFT +=` line. The subsequent atomic block already
accounts for the full FFT+peak time.

**To revert:** Restore `timeFFT += (double)((clock()-time1))/(CLOCKS_PER_SEC);` between
`ampMatchEdge(...)` and the `// Step 4` comment.

---

### 2. Unchecked fread return value (LSB path)
**Problem:** Same as Strack/speckleTrack.c — `size_t rv = fread(...)` declared but unused.

**Fix:** Changed to `(void)fread(...)`.

**To revert:** Change `(void)fread(...)` back to `size_t rv = fread(...)`.

---

## Cullst/cullStats.c

### Division by zero when ngood == 0
**Problem:**
```c
meanA /= (double)ngood;   // crashes/NaN if ngood == 0
meanR /= (double)ngood;
if (ngood > 8) { ... }
```
The division occurred before the `ngood > 8` guard. A window with no valid pixels
(ngood == 0) produces division by zero (IEEE 754 yields ±Inf, not a trap), and the
results are meaningless since they're only used inside `if (ngood > 8)`.

**Fix:** Moved both divisions inside the `if (ngood > 8)` block so they only execute
when `ngood` is known to be positive.

**To revert:** Move `meanA /= (double)ngood;` and `meanR /= (double)ngood;` back to
before the `if (ngood > 8)` line.

---

## Cullls/cullLSStats.c

### Division by zero when ngood == 0
**Problem:** Identical pattern to Cullst/cullStats.c — `meanY /= (double)ngood;` and
`meanX /= (double)ngood;` executed before the `if (ngood > 8)` guard.

**Fix:** Moved both divisions inside the `if (ngood > 8)` block.

**To revert:** Move `meanY /= (double)ngood;` and `meanX /= (double)ngood;` back to
before the `if (ngood > 8)` line.

---

## Cullst/cullSmooth.c

### malloc size off by one (operator precedence)
**Problem:** `wA = (float *)malloc(sizeof(float) * sA + 1);`
Due to C operator precedence, this evaluates as `(sizeof(float) * sA) + 1` = `4*sA + 1`
bytes — enough for `sA` floats plus 1 byte. The subsequent loop writes `sA + 1` floats
(indices `-sA/2` through `sA/2` inclusive), overflowing the allocation by 3 bytes.

**Fix:** `wA = (float *)malloc(sizeof(float) * (sA + 1));`

**To revert:** Remove the parentheses: `malloc(sizeof(float) * sA + 1)`.

---

## Strack/parseTrack.c

### 1. parFile1/parFile2 uninitialised when VRT path is taken
**Problem:** When VRT files are found for image1/image2, `trackPar->parFile1` and
`trackPar->parFile2` are never assigned. `printTrackPar()` then prints them
unconditionally with `%s`, which is undefined behaviour (garbage pointer → likely segfault).

**Fix:** Added `trackPar->parFile1 = NULL; trackPar->parFile2 = NULL;` at the start of
`parseTrack()` before the VRT/par branching logic.

**To revert:** Remove the two NULL initialisation lines.

---

### 2. edgePadR assigned instead of edgePadA in fallback
**Problem:** When `sscanf` fails to parse the `edgeA` parameter:
```c
trackPar->edgePadR = trackPar->edgePad;   // typo: should set edgePadA
```
This left `edgePadA` uninitialised and set `edgePadR` a second time.

**Fix:** Changed to `trackPar->edgePadA = trackPar->edgePad;`

**To revert:** Change `edgePadA` back to `edgePadR` in the sscanf fallback for edgeA.

---

## Strack/strack.c

### Wrong error message for image 2 parse failure
**Problem:** `error("Could not read par or vrt file for image 1\n")` was called in the
image-2 parsing block (where neither `parFile2` nor `vrtFile2` was provided).

**Fix:** Changed the message to `"Could not read par or vrt file for image 2\n"`.

**To revert:** Change "image 2" back to "image 1" in the second `error(...)` call.

---

## Strackw/strackw.c

### Wrong error message for image 2 parse failure
**Problem:** Identical to Strack/strack.c — same wrong "image 1" message in the
image-2 block.

**Fix:** Changed to `"Could not read par or vrt file for image 2\n"`.

**To revert:** Change "image 2" back to "image 1".

---

## Cullst/writeCullData.c

### fwrite arguments in wrong order
**Problem:**
```c
return fwrite(ptr, nitems, size, fp);
```
`fwrite` signature is `fwrite(ptr, element_size, count, fp)`. The arguments `nitems`
(count) and `size` (element size) were swapped. Since `nitems` is the total byte count
and `size` is the element size, this caused `fwrite` to attempt writing `nitems` elements
of `size` bytes each — likely writing far more data than intended.

**Fix:** Changed to `return fwrite(ptr, size, nitems, fp);`

**To revert:** Swap `size` and `nitems` back to `nitems, size`.

---

## Cullst/cullIslands.c

### Dead zero-byte malloc for nB arrays before nB is computed
**Problem:**
```c
holes[i].xb = (int32_t *)malloc(holes[i].nB * sizeof(int));
holes[i].yb = (int32_t *)malloc(holes[i].nB * sizeof(int));
```
`nB` (border pixel count) is always 0 at this point in Pass 2 — it is only populated
in Pass 3. These are zero-byte allocations whose returned pointers are never written to
before being overwritten (or the arrays are unused). `malloc(0)` behaviour is
implementation-defined (may return NULL or a unique non-dereferenceable pointer).

**Fix:** Removed both lines. The `xh`/`yh` allocations (which use `nH`, correctly computed
in Pass 2) are unchanged.

**To revert:** Restore the two `holes[i].xb` and `holes[i].yb` malloc lines before the
`xh`/`yh` allocations.
