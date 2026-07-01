# strackw — SAR Amplitude Correlation Tracker

## Purpose

**strackw** estimates sub-pixel range and azimuth offsets between two co-registered SAR
images using normalised amplitude cross-correlation only. Unlike `strack`, it never
attempts complex cross-correlation and does not require a coherent interferometric signal.
It is designed for applications such as ice-sheet velocity mapping where speckle
decorrelates between acquisitions.

---

## Usage

```
strackw [-integerComplex] [-LSB] parFile
```

### Options

| Option | Description |
|--------|-------------|
| `-integerComplex` | Input images are int16 complex (default: float32 complex) |
| `-LSB` | Read/write in little-endian byte order (default: big-endian MSB) |

---

## Parameter File Format

Same RDF (keyword = value) ASCII format as `strack`. Comments begin with `;` or `#`.
Only the amplitude-tracking keywords are used; complex-match, baseline, and
interferogram keywords are ignored.

### Image and Geometry

| Keyword | Description |
|---------|-------------|
| `image1` | Path to primary SLC image |
| `image2` | Path to secondary SLC image |
| `image1par` | Legacy CW `.par` file for image 1 |
| `image2par` | Legacy CW `.par` file for image 2 |
| `image1vrt` | VRT file for image 1 (alternative to `image1par`) |
| `image2vrt` | VRT file for image 2 (alternative to `image2par`) |

### Output

| Keyword | Description |
|---------|-------------|
| `outputfile` | Base name for all output files |

### Search Grid

| Keyword | Description |
|---------|-------------|
| `rstart` | Starting range pixel (single-look) |
| `astart` | Starting azimuth line (single-look) |
| `deltar` | Range step between output pixels |
| `deltaa` | Azimuth step between output lines |
| `nR` | Number of output columns |
| `nA` | Number of output rows |

### Correlation Window

| Keyword | Description |
|---------|-------------|
| `wra` | Total search window width in range (pixels) — see geometry diagram below |
| `waa` | Total search window height in azimuth (pixels) |
| `edgepad` | Symmetric edge pad (sets both `edgeR` and `edgeA`) |
| `edgeR` | Edge pad in range (pixels) — defines search radius (see diagram) |
| `edgeA` | Edge pad in azimuth (pixels) |
| `navgR` | Range averaging factor applied before correlation |
| `navgA` | Azimuth averaging factor applied before correlation |
| `scalefactor` | Scale factor applied to output offsets (number of looks) |

### Initial Offset Polynomial

| Keyword | Description |
|---------|-------------|
| `initShift` | `1` = read polynomial from file; `0` = zero shift |
| `initialoffsetfile` | ASCII file with `range_offset_polynomial:` and `azimuth_offset_polynomial:` coefficients |

### Mask

| Keyword | Description |
|---------|-------------|
| `maskfile` | Binary byte mask (0 = skip, 1 = process) |
| `maskgeodat` | Geodat or `.dat` header for the mask |
| `offsetmaskvrt` | VRT file for an offset-grid mask |

---

## Search Window and Chip Geometry

The core concept is the relationship between the **search window** (large patch from
image 2) and the **search chip** (small patch from image 1 that slides across the search
window to find the best match).

### Parameter Definitions

```
  wRa       total search window width  (range,   pixels, from parameter file)
  wAa       total search window height (azimuth, pixels, from parameter file)
  edgePadR  edge pad in range          (from edgeR or edgepad keyword)
  edgePadA  edge pad in azimuth        (from edgeA or edgepad keyword)

  Derived:
    chipR  =  wRa - 2 * edgePadR        (search chip width in range)
    chipA  =  wAa - 2 * edgePadA        (search chip height in azimuth)

  Maximum detectable displacement:
    maxShiftR  =  edgePadR  (range pixels)
    maxShiftA  =  edgePadA  (azimuth pixels)
```

### Window Layout Diagram

Image 2 is centred on the position predicted by the initial-offset polynomial. The
search chip from image 1 is slid across the interior of the image 2 search window.

```
  Image 2 search window  (wRa x wAa pixels, centred on predicted position)
  +---------------------------------------------------------------------+
  |                                                                     |
  |  <--edgePadR--><------------- chipR = wRa - 2*edgePadR ----------->|<--edgePadR-->
  |                                                                     |
  |    ^           +---------------------------------------+            |
  |    |           |                                       |            |
  | edgePadA       |   Image 1 search chip  (chipR x chipA)|           |
  |    |           |   Correlated at every position within |            |
  |    v           |   the +-edgePadR / +-edgePadA region  |            |
  |    ^           |                                       |            |
  |    |           +---------------------------------------+            |
  |  chipA                                                              |
  |    |              chip slides +/-edgePadR  in range                |
  |    v              chip slides +/-edgePadA  in azimuth              |
  |    ^                                                                |
  | edgePadA                                                            |
  |    v                                                                |
  +---------------------------------------------------------------------+
  <-------------------------------- wRa -------------------------------->
```

### Extra Internal Padding (NFAST/4 = 5 pixels)

At the start of `corrTrackFast`, the code adds `NFAST/4 = 5` pixels to both
`edgePadR` and `edgePadA` **beyond** what the parameter file specifies:

```c
trackPar->edgePadR += NFAST / 4;   /* +5 pixels */
trackPar->edgePadA += NFAST / 4;   /* +5 pixels */
```

This extra pad is required so that the NFAST-wide oversampling window used to locate
the correlation peak cannot run off the edge of the correlation surface. The
consequence is that the **effective chip size** is smaller than naive subtraction
suggests:

```
  Internal edgePadR  =  edgePadR (parfile) + 5
  Internal edgePadA  =  edgePadA (parfile) + 5

  Actual chip size correlated:
    internalChipR  =  wRa - 2 * (edgePadR + 5)
    internalChipA  =  wAa - 2 * (edgePadA + 5)

  Effective detectable displacement (unchanged by internal pad):
    maxShiftR  =  edgePadR  (parfile value)
    maxShiftA  =  edgePadA  (parfile value)
```

**Example:** `wra = 64`, `edgeR = 20`
- Naive chip width: 64 - 2x20 = 24 pixels
- Actual correlated chip width: 64 - 2x(20+5) = **14 pixels**

The NFAST/4 overhead must be budgeted when choosing `wra`/`waa` and `edgeR`/`edgeA`.
A safe rule is to add at least 5 pixels to each edge pad beyond the minimum needed for
the desired search radius.

---

## Algorithm — `corrTrackFast`

Constants: NOVER = 16, NFAST = 20, OS = 2 (amplitude oversampling factor).

For each output pixel $(i, j)$:

1. **Centre position** in image 1:
   $r_1 = r_\text{start} + j \cdot \Delta r$,
   $a_1 = a_\text{start} + i \cdot \Delta a$.

2. **Predict image-2 position** $(r_2, a_2)$ from the initial-offset polynomial (or
   interpolated offset table), rounded to the nearest integer pixel.

3. **Mask check** — skip pixel and write `−LARGEINT` if mask value is 0.

4. **Read patches** from SLC line buffers (`NBUFFERLINES = 2000` lines):
   - Image 1: `wRa × wAa` centred on $(r_1, a_1)$.
   - Image 2: `wRa × wAa` centred on $(r_2, a_2)$.

5. **Remove Doppler carrier** (`estDopCarrier1`) to whiten the azimuth spectrum before
   detection.

6. **Forward FFT** both patches, **zero-pad** to `wRa·OS × wAa·OS` (factor OS = 2),
   **inverse FFT** → oversampled detected (amplitude squared) patches:
   - `dataS` (`wRa·OS × wAa·OS`) — image 2 search window (large patch).
   - `dataR` (`chipR·OS × chipA·OS`) — image 1 search chip, extracted from the
     interior after removing the `edgePadR`/`edgePadA` border.

7. **Normalised cross-correlation** (`correlateFast`):

   a. Compute mean $\bar{d}_R$ and variance $\sigma_R^2$ of `dataR` (chip, single values).

   b. For each possible chip position $(s, l)$ within `dataS`, compute the local
      running mean $\bar{d}_S(s,l)$ and variance $\sigma_S^2(s,l)$ using 2-D
      sliding-window box filters (separable column then row passes).

   c. Compute the normalised cross-correlation via FFT convolution:

$$C(s,l) = \frac{\displaystyle\sum_{p,q}\bigl(d_R(p,q)-\bar{d}_R\bigr)\bigl(d_S(s+p,\,l+q)-\bar{d}_S(s,l)\bigr)}{\sqrt{\sigma_R^2 \cdot \sigma_S^2(s,l)}}$$

   d. Find the integer-pixel peak $(i_\text{max}, j_\text{max})$ in the correlation
      surface (`corrResult`).

8. **Sub-pixel oversampling** — extract an NFAST×NFAST region around the integer peak
   into `cFast`, zero-pad to NFAST·NOVER × NFAST·NOVER, locate sub-pixel peak.

9. **Sub-pixel shift** from the oversampled peak:

$$r_\text{shift} = \frac{j_\text{max} - \text{edgePadR} \cdot \text{OS} \cdot \text{NOVER}}{\text{NOVER} \cdot \text{OS}}$$

$$a_\text{shift} = \frac{i_\text{max} - \text{edgePadA} \cdot \text{OS} \cdot \text{NOVER}}{\text{NOVER} \cdot \text{OS}}$$

10. **Save** final offset = initial predicted shift − measured sub-pixel shift,
    multiplied by `scaleFactor`. Failed pixels: `−LARGEINT`.

---

## Output Files

All paths derived from `outputfile`:

| Suffix | Type | Contents |
|--------|------|----------|
| `.dr` | float32, nA×nR | Range offsets (pixels × scaleFactor, MSB) |
| `.da` | float32, nA×nR | Azimuth offsets (pixels × scaleFactor, MSB) |
| `.cc` | float32, nA×nR | Peak normalised cross-correlation value |
| `.mt` | byte, nA×nR | Match type (0 = BAD, 2 = AMPMATCH) |
| `.dat` | ASCII | Offset grid header: `r0 a0 nr na deltaR deltaA` |
| `.vrt` | GDAL VRT | Three-band VRT: RangeOffsets, AzimuthOffsets, Correlation |
| `.mt.vrt` | GDAL VRT | Single-band VRT for match type |

Failed pixels are written as `−LARGEINT` (≈ −2×10⁹) in `.dr` and `.da`.

---

## Dependencies

| Function | Source | Purpose |
|----------|--------|---------|
| `corrTrackFast` | `Strackw/corrTrackFast.c` | Main amplitude correlation loop |
| `parseTrack` | `Strack/parseTrack.c` | RDF parameter file parser |
| `parseInitialOffsets` | `Strack/parseInitialOffsets.c` | Initial offset polynomial reader |
| `getMask` | `Strack/getMask.c` | Byte mask reader |
| `writeVrtFile` | `Strack/writeVrt.c` | VRT output writer |
| `writeOffsets` | `Strack/sTrackOut.c` | Per-row binary output writer |
| `readOldPar` | `common/readOldPar.c` | Legacy CW `.par` reader |
| `parseSLCVrtNew` | `common/parseSLCVrtNew.c` | SLC VRT metadata reader |
| `writeSingleVRT` | `gdalIO/gdalIO/gdalIO.c` | GDAL VRT writer |
| FFTW | — | FFT-based amplitude cross-correlation and peak oversampling |
| GDAL | — | VRT-described SLC image input |
