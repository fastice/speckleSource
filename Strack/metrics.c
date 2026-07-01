typedef struct {
    int i_peak;          /* row index of peak */
    int j_peak;          /* col index of peak */
    double p1;           /* peak value */
    double p2;           /* 2nd peak value (outside exclusion box) */
    double p2_over_p1;   /* ambiguity ratio */
    double mean_sl;      /* sidelobe mean */
    double std_sl;       /* sidelobe std */
    double psr;          /* (p1 - mean_sl) / std_sl */
    int fwhm_i;          /* FWHM width in rows (approx, pixels) */
    int fwhm_j;          /* FWHM width in cols (approx, pixels) */
    int n_sidelobe;      /* number of sidelobe samples used */
} CorrMetrics;


/*
  Compute FWHM width along a 1D cut (array values assumed >=0).
  center index is c, length is n. Returns width in samples (>=1).
  This is a simple integer FWHM: finds first crossing below half on each side.
*/

static int fwhm_1d_from_2d_row(float **corr, int n_i, int j_peak, int i_center)
{
    double half = 0.5 * corr[i_center][j_peak];
    int left = i_center;
    int right = i_center;
    if (half <= 0.0) return 0;
    while (left > 0 && corr[left][j_peak] >= half) left--;
    while (right < n_i - 1 && corr[right][j_peak] >= half) right++;
    return right - left + 1;
}

static int fwhm_1d_from_2d_col(float **corr, int n_j, int i_peak, int j_center)
{
    double half = 0.5 * corr[i_peak][j_center];
    int left = j_center;
    int right = j_center;

    if (half <= 0.0) return 0;
    while (left > 0 && corr[i_peak][left] >= half) left--;
    while (right < n_j - 1 && corr[i_peak][right] >= half) right++;
    return right - left + 1;
}

/*
  Compute correlation quality metrics.

  corr:  correlation surface (e.g., abs(crosscorr)) flattened row-major
  n_i:   number of rows (az dimension of corr surface)
  n_j:   number of cols (rg dimension of corr surface)
  ex_i, ex_j: half-size of exclusion box around main peak when computing sidelobes and p2
             (e.g., 2..4 depending on expected peak width)

  Returns CorrMetrics (psr can be huge if std is tiny; handle thresholds in caller).
*/
CorrMetrics compute_corr_metrics(float **corr, int32_t i_peak, int32_t j_peak,
    							int32_t n_i, int32_t n_j,int32_t ex_i,	int32_t ex_j)
{
    CorrMetrics m;
    int i, j;
    int i0, i1, j0, j1;
    double p1 = -DBL_MAX;
    double p2 = -DBL_MAX;
    double sum = 0.0L;
    double sum2 = 0.0L;
    int nsl = 0;
    /* Initialize output */
    m.i_peak = i_peak; m.j_peak = j_peak;
    m.p1 = 0.0; m.p2 = 0.0; m.p2_over_p1 = 0.0;
    m.mean_sl = 0.0; m.std_sl = 0.0; m.psr = 0.0;
    m.fwhm_i = 0; m.fwhm_j = 0;
    m.n_sidelobe = 0;
 	if (!corr || n_i <= 0 || n_j <= 0) return m;
	/* 1) Get peak value */
	p1 = corr[i_peak][j_peak];
    /* Define exclusion box around peak */
    i0 = MAX(0, i_peak - ex_i);
    i1 = MIN(n_i - 1, i_peak + ex_i);
    j0 = MAX(0, j_peak - ex_j);
    j1 = MIN(n_j - 1, j_peak + ex_j);
    /* 2) Compute sidelobe mean/std (exclude peak box) and find 2nd peak */
	for (i = 0; i < n_i; ++i) {
		for (j = 0; j < n_j; ++j) {
			double v = corr[i][j];
			/* outside exclusion box => sidelobe region */
			if (i < i0 || i > i1 || j < j0 || j > j1) {
				sum += (long double)v;
				sum2 += (long double)v * (long double)v;
				nsl++;
				if (v > p2) p2 = v;
			}
		}
	}
	m.n_sidelobe = nsl;
	
	if (nsl > 0) {
		double mean = sum / (double)nsl;
		double var = (sum2 / (double)nsl) - mean * mean;
		if (var < 0.0L) var = 0.0L;
		m.mean_sl = (double)mean;
		m.std_sl = (double)sqrt((double)var);
	} else {
		m.mean_sl = 0.0;
		m.std_sl = 0.0;
	}
	
    /* 3) PSR and p2/p1 */
    m.i_peak = i_peak;
    m.j_peak = j_peak;
    m.p1 = p1;
    m.p2 = (p2 > -DBL_MAX/2 ? p2 : 0.0);
	m.p2_over_p1 = (p1 > 1e-12 ? m.p2 / p1 : 0.0);
	m.psr = (m.std_sl > 1e-12 ? (p1 - m.mean_sl) / m.std_sl : 0.0);  
    /* 4) FWHM estimates from 1D cuts through the peak */
	m.fwhm_i = fwhm_1d_from_2d_row(corr, n_i, j_peak, i_peak);
    m.fwhm_j = fwhm_1d_from_2d_col(corr, n_j, i_peak, j_peak);
    return m;
}