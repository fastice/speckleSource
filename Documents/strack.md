# strack / strackw — SAR Speckle / Amplitude Offset Trackers

## Purpose

**strack** estimates sub-pixel range and azimuth offsets between two co-registered SAR
images using complex cross-correlation (phase coherence), falling back to amplitude
correlation when coherence is insufficient.

**strackw** is a stripped-down variant that uses amplitude-only cross-correlation
(`corrTrackFast`) and is optimised for incoherent targets such as wide-area ice-sheet
tracking where the complex signal decorrelates.

Both programs output the same file set and share the same RDF parameter-file format.

---

## Usage

```
strack   [options] paramFile
strackw  [options] paramFile
```

### strack Options

| Option | Description |
|--------|-------------|
| `-noComplex` | Skip complex matching; use amplitude only |
| `-singleAmp` | Single amplitude-match attempt (no large-window fallback) |
| `-gauss` | Use Gaussian fit for peak location instead of oversampling |
| `-integerComplex` | Input images are int16 complex (default: float32 complex) |
| `-noHanning` | Suppress Hanning window on patches before FFT |
| `-checkAzFocus` | Enable azimuth-defocus quality filter; write `.azd` file |
| `-LSB` | Write output in little-endian byte order (default: big-endian MSB) |

### strackw Options

| Option | Description |
|--------|-------------|
| `-integerComplex` | Input images are int16 complex (default: float32 complex) |
| `-LSB` | Write output in little-endian byte order |

---

## Parameter File Format

Both programs use an **RDF** (keyword = value) ASCII parameter file. Lines beginning
with `#` or `;` are comments. Keyword matching is case-insensitive.

### Image and Geometry

| Keyword | Description |
|---------|-------------|
| `image1` | Path to primary SLC image |
| `image2` | Path to secondary SLC image |
| `image1par` | Geodat file for image 1 (CW `.par` format, or use `image1vrt`) |
| `image2par` | Geodat file for image 2 (CW `.par` format, or use `image2vrt`) |
| `image1vrt` | VRT file for image 1 (alternate to `image1par`) |
| `image2vrt` | VRT file for image 2 (alternate to `image2par`) |
| `intfile` | Complex interferogram file (int16 pairs, big-endian) — strack only |
| `intgeodat` | Geodat file for the interferogram — strack only |

### Output

| Keyword | Description |
|---------|-------------|
| `outputfile` | Base name for all output files (suffixes added automatically) |

### Search Grid

| Keyword | Default | Description |
|---------|---------|-------------|
| `rstart` | — | Starting range pixel (single-look) |
| `astart` | — | Starting azimuth line (single-look) |
| `deltar` | — | Range step between output pixels (single-look pixels) |
| `deltaa` | — | Azimuth step between output lines (single-look pixels) |
| `nR` | — | Number of output columns |
| `nA` | — | Number of output rows |

### Correlation Window

| Keyword | Description |
|---------|-------------|
| `wr` | Complex-match window half-width in range (pixels) — full width = `wr` |
| `wa` | Complex-match window half-height in azimuth (pixels) |
| `wra` | Amplitude-match window width in range |
| `waa` | Amplitude-match window height in azimuth |
| `edgepad` | Symmetric edge pad (search radius limit); overridden by `edgeR`/`edgeA` |
| `edgeR` | Edge pad in range (pixels) |
| `edgeA` | Edge pad in azimuth (pixels) |
| `navgR` | Range averaging factor |
| `navgA` | Azimuth averaging factor |
| `scalefactor` | Scale factor applied to output offsets (usually number of looks) |

### Baseline / Phase Correction (strack only)

| Keyword | Description |
|---------|-------------|
| `baseparams` | ASCII file with baseline: single line `Bn Bp dBn dBp` (metres) |
| `initShift` | Flag (`1`/`0`) to enable polynomial initial-offset read |
| `initialoffsetfile` | ASCII file with initial offset polynomials (see below) |

### Masks

| Keyword | Description |
|---------|-------------|
| `maskfile` | Binary byte mask (0 = skip, ≥1 = process) co-registered to images |
| `maskgeodat` | Geodat or `.dat` header for the mask file |
| `offsetmaskfile` | Alternative mask on the offset grid |
| `offsetmaskdat` | Dat/VRT header for the offset-grid mask |
| `offsetmaskvrt` | VRT file for the offset-grid mask (preferred over `offsetmaskdat`) |

---

## Input File Formats

### SLC Images

Two modes controlled by `-integerComplex`:

- **float32 complex** (default): interleaved float32 pairs `(re, im)`, big-endian (or
  VRT-described via `image1vrt`)
- **int16 complex** (`-integerComplex`): interleaved int16 pairs, big-endian

### Interferogram (`intfile`) — strack only

Big-endian int16 complex, `nr × na` samples, 4 bytes per pixel (2 bytes real + 2 bytes
imaginary). Geometry described by `intgeodat`. Used for phase flattening of the complex
cross-correlation before sub-pixel peak location.

### Baseline Parameter File (`baseparams`) — strack only

Single ASCII line:

```
Bn   Bp   dBn   dBp
```

where $B_n$ and $B_p$ are the normal and parallel baselines at scene start (metres), and
$\delta B_n$, $\delta B_p$ are their linear rates of change per scene line.
(Same format as `computeBaseline` output.)

### Initial Offset Polynomial File (`initialoffsetfile`)

ASCII file with two labelled lines (order not required):

```
range_offset_polynomial:   c0  c1  c2
azimuth_offset_polynomial: c0  c1  c2
```

The polynomial is evaluated at each output pixel $(r, a)$ as:

$$
\delta r = c_0 + c_1 \cdot r + c_2 \cdot a
$$

(and similarly for azimuth). Provides the predicted image-to-image shift used to
position patch 2 before matching.

### Mask File (`maskfile`)

Binary byte array, same dimensions as the images (or the offset grid if an offset mask).
Values: `0` = skip this pixel; `≥1` = attempt matching. Read big-endian (byte, so no
swap needed). Geometry described by `maskgeodat` (geodat or `.dat` header with
`r0 a0 nr na deltaR deltaA`).

---

## Output Files

All output file paths are derived from `outputfile` by appending suffixes:

| Suffix | Type | Contents |
|--------|------|----------|
| `.dr` | float32, nA×nR | Range offsets (pixels, scaled by `scaleFactor`) |
| `.da` | float32, nA×nR | Azimuth offsets (pixels, scaled by `scaleFactor`) |
| `.cc` | float32, nA×nR | Peak correlation value |
| `.mt` | byte, nA×nR | Match type (see codes below) |
| `.azd` | float32, nA×nR | Azimuth defocus ratio σ_az/σ_rg (strack `-checkAzFocus` only) |
| `.dat` | ASCII | Offset grid header: `r0 a0 nr na deltaR deltaA` |
| `.vrt` | GDAL VRT | Three-band VRT for `.dr`, `.da`, `.cc` (bands: RangeOffsets, AzimuthOffsets, Correlation) |
| `.mt.vrt` | GDAL VRT | Single-band VRT for `.mt` (band: MatchType) |

Failed pixels are written as `−LARGEINT` (≈ −2×10⁹) in `.dr` and `.da`.

### Match Type Codes (`.mt`)

| Value | Meaning |
|-------|---------|
| `0` | BAD — no match |
| `1` | CMATCH — complex cross-correlation |
| `2` | AMPMATCH — amplitude correlation (standard window) |
| `3` | AMPMATCHLARGE — amplitude correlation (2× larger window) |

### VRT Metadata Fields

The `.vrt` file carries GDAL dataset-level metadata for downstream use:

| Key | Description |
|-----|-------------|
| `r0`, `a0` | First range/azimuth pixel (single-look, after scaleFactor) |
| `deltaR`, `deltaA` | Step size in range/azimuth (after scaleFactor) |
| `geo1`, `geo2` | Geodat paths for image 1 and image 2 |
| `Image1`, `Image2` | Paths to the SLC image files |
| `mask` | Path to the mask file |
| `wR`, `wA` | Complex-match window dimensions |
| `wRa`, `wAa` | Amplitude-match window dimensions |
| `scaleFactor` | Scale factor applied to offsets |
| `ByteOrder` | `MSB` or `LSB` |
| `sigmaStreaks`, `sigmaRange` | Reserved (written as 0) |

---

## Algorithm

### strack — `speckleTrack`

Constants: NOVER = 12 (oversampling factor), OS = 2 (amplitude oversample).

For each output pixel $(i, j)$:

1. **Find patch-2 position** using the initial offset polynomial evaluated at $(r_1, a_1)$.

2. **Mask check**: skip if `maskValue` returns 0.

3. **Complex matching** (unless `-noComplex`):
   - Read SLC patches (size `wR × wA`) from both images into buffers.
   - Apply Hanning window (unless `-noHanning`).
   - Phase-flatten with baseline model and/or interferogram (if available).
   - 2-D FFT cross-correlation: compute power spectrum, oversample to `wR·NOVER × wA·NOVER`.
   - Locate peak with `cmpTrackFast` (or Gaussian fit with `-gauss`).
   - Run Gaussian 1-D fits in range and azimuth to compute $\sigma_{rg}$, $\sigma_{az}$.
   - If `-checkAzFocus`: reject if $\sigma_{az}/\sigma_{rg}$ exceeds threshold; write ratio to `.azd`.
   - Sub-pixel shift:
     $$r_{\text{shift}} = \frac{j_{\text{peak}} - w_R/2 \cdot \text{NOVER} \cdot \text{OSA}}{\text{NOVER}}
     \qquad
     a_{\text{shift}} = \frac{i_{\text{peak}} - w_A/2 \cdot \text{NOVER} \cdot \text{OSA}}{\text{NOVER}}$$
   - Accepted if peak is within ±15 % of window and correlation exceeds threshold.

4. **Amplitude matching** (fallback, or primary with `-noComplex`):
   - Try standard window (`wRa × wAa`); correlation threshold 0.07.
   - If still no match (and not `-singleAmp`): try 2× larger window; threshold 0.028.
   - Accept if peak is within ±15 % of window and correlation exceeds threshold.

5. **Save** range and azimuth offsets, correlation, and match type.
   Failed pixels: offset = `−LARGEINT`.
   All offsets multiplied by `scaleFactor` before writing.

6. Write line to output files after each azimuth row; write `.dat` and `.vrt` at end.

---

### strackw — `corrTrackFast`

Constants: NOVER = 16, OS = 2.

Amplitude-only; no complex patches, no interferogram, no baseline phase correction.

1. Detect amplitude patches (square-law detect complex SLC) for both images.
2. Normalised cross-correlation of detected patches:
   $$C(i,j) = \frac{\sum (s_1 - \bar{s}_1)(s_2 - \bar{s}_2)}{\sqrt{\sum(s_1-\bar{s}_1)^2 \sum(s_2-\bar{s}_2)^2}}$$
3. Oversample correlation surface (FFT zero-padding × NOVER) to locate peak.
4. Peak must lie within `edgePadR` / `edgePadA` of window centre; otherwise rejected.
5. Same output file set as strack; `.mt` values are only 0 or 2 (BAD / AMPMATCH).

---

## Supporting Functions

### `parseTrack`(paramFile, trackPar) → void
Reads RDF parameter file into `TrackParams`; sets all keywords listed above.
Derives output filenames by appending `.dr`, `.da`, `.cc`, `.mt`, `.dat`, `.azd`,
`.vrt`, `.mt.vrt` to `outputfile`.  
*Calls:* `getKeyWord`, `openInputFile` (RDF routines)

---

### `parseBase`(trackPar) → void
Reads `baseparams` file: single ASCII line → `Bn Bp dBn dBp` (float64).
Sets `trackPar->dBnQ = dBpQ = 0`.  
*Calls:* `getDataString`, `openInputFile`

---

### `parseInitialOffsets`(trackPar) → void
Reads `initialoffsetfile`; looks for lines starting with
`range_offset_polynomial:` and `azimuth_offset_polynomial:`, each followed by 3
float64 coefficients. Initialises shift polynomials to zero if `polyShift = FALSE`.

---

### `getInt`(trackPar) → void
Reads complex interferogram (`intfile`) into memory:
- Geometry from `intgeodat` via `parseInputFile`.
- Data: int16 complex pairs, big-endian. Converts to float32 fftw_complex.
- If file is shorter than geodat dimensions, missing lines are filled with (1, 1).

*Calls:* `parseInputFile`, `freadBS`

---

### `getMask`(trackPar) → void
Two modes:
- **VRT** (`offsetmaskvrt`): reads via GDAL, extracts `r0`, `a0`, `deltaR`, `deltaA` from VRT metadata.
- **File** (`maskfile` + `maskgeodat` or `.dat` header): reads byte array with `freadBS`.

*Calls:* `parseInputFile`, `GDALOpen`, `readDataSetMetaData`, `freadBS`

---

### `writeVrtFile`(trackPar) → void
Writes `.mt.vrt` (1 band, byte) and `.vrt` (3 bands: RangeOffsets float32,
AzimuthOffsets float32, Correlation float32) with full metadata dictionary.
`geo2` is constructed from the directory of `imageFile2` and the basename of `intGeodat`.  
*Calls:* `writeSingleVRT`, `insert_node`

---

### `writeOffsets`(i, trackPar, …) → void
Writes one azimuth row of `.dr`, `.da`, `.cc`, `.azd`, `.mt` to open file handles.
Called after each output row so results stream to disk during processing.

---

## Correlation Quality Metrics (`metrics.c`)

`compute_corr_metrics(corr, i_peak, j_peak, n_i, n_j, ex_i, ex_j)` returns a
`CorrMetrics` struct:

| Field | Description |
|-------|-------------|
| `p1` | Peak correlation value |
| `p2` | Second-highest peak outside exclusion box |
| `p2_over_p1` | Ambiguity ratio $p_2 / p_1$ |
| `mean_sl` | Mean sidelobe level |
| `std_sl` | Standard deviation of sidelobes |
| `psr` | Peak-to-sidelobe ratio $(p_1 - \bar{s}) / \sigma_s$ |
| `fwhm_i`, `fwhm_j` | FWHM width of correlation peak in azimuth / range (pixels) |

---

## Dependencies

| Function | Source | Purpose |
|----------|--------|---------|
| `parseTrack` | `Strack/parseTrack.c` | RDF parameter file parser |
| `parseBase` | `Strack/parseBase.c` | Baseline parameter reader |
| `parseInitialOffsets` | `Strack/parseInitialOffsets.c` | Initial offset polynomial reader |
| `getInt` | `Strack/getInt.c` | Complex interferogram reader |
| `getMask` | `Strack/getMask.c` | Byte mask reader |
| `speckleTrack` | `Strack/speckleTrack.c` | Complex+amplitude matching loop (strack) |
| `corrTrackFast` | `Strackw/corrTrackFast.c` | Amplitude-only matching loop (strackw) |
| `writeVrtFile` | `Strack/writeVrt.c` | VRT output writer |
| `writeOffsets` | `Strack/sTrackOut.c` | Per-row binary output writer |
| `compute_corr_metrics` | `Strack/metrics.c` | Correlation quality metrics |
| `parseInputFile` | `common/parseInputFile.c` | SAR geodat reader |
| `readOldPar` | `common/readOldPar.c` | Legacy CW `.par` reader |
| `parseSLCVrtNew` | `common/parseSLCVrtNew.c` | SLC VRT metadata reader |
| `writeSingleVRT` | `gdalIO/gdalIO/gdalIO.c` | GDAL VRT writer |
| `readDataSetMetaData` | `gdalIO/gdalIO/gdalIO.c` | GDAL metadata reader |
| FFTW | — | FFT-based cross-correlation |
| GDAL | — | Raster I/O for VRT-described inputs |
