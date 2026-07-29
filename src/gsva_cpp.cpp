#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// ============================================================================
// Kernel CDF estimation (matches GSVA kernel_estimation.c exactly)
// Uses precomputed CDF lookup table for Gaussian kernel to match original
// ============================================================================

// Constants matching GSVA's kernel_estimation.c
#define SIGMA_FACTOR 4.0
#define PRECOMPUTE_RESOLUTION 10000
#define MAX_PRECOMPUTE 10.0

static double precomputed_cdf[PRECOMPUTE_RESOLUTION + 1];
static int is_precomputed = 0;

static void initCdfs() {
    double divisor = PRECOMPUTE_RESOLUTION * 1.0;
    for (int i = 0; i <= PRECOMPUTE_RESOLUTION; ++i)
        precomputed_cdf[i] = Rf_pnorm5(MAX_PRECOMPUTE * ((double)i) / divisor,
                                        0.0, 1.0, 1, 0);
}

static inline double precomputedCdf(double x, double sigma) {
    double v = x / sigma;
    if (v < (-1 * MAX_PRECOMPUTE)) {
        return 0;
    } else if (v > MAX_PRECOMPUTE) {
        return 1;
    } else {
        double cdf = precomputed_cdf[(int)(fabs(v) / MAX_PRECOMPUTE * PRECOMPUTE_RESOLUTION)];
        if (v < 0) {
            return 1.0 - cdf;
        } else {
            return cdf;
        }
    }
}

// Compute sd using raw pointer (Welford-like single-pass)
static inline double row_sd(const double* x, int n) {
    if (n < 2) return NA_REAL;
    double mean_val = 0.0;
    for (int i = 0; i < n; i++) mean_val += x[i];
    mean_val /= n;
    double ss = 0.0;
    for (int i = 0; i < n; i++) {
        double d = x[i] - mean_val;
        ss += d * d;
    }
    return std::sqrt(ss / (n - 1));
}

// Poisson kernel CDF for a chunk of gene rows
// Process genes from row_start to row_end (0-based, exclusive)
// Uses raw pointers, contiguous row buffer for cache efficiency,
// and caches ppois results for repeated expression values (common in sparse data)
// [[Rcpp::export]]
NumericMatrix gene_cdf_poisson_cpp(NumericMatrix expr, int row_start = 0, int row_end = -1) {
    int n_genes = expr.nrow();
    int n_samples = expr.ncol();
    if (row_end < 0) row_end = n_genes;
    int n_rows = row_end - row_start;

    NumericMatrix cdf_matrix(n_rows, n_samples);
    double* out = REAL(cdf_matrix);
    double* X = REAL(expr);

    // Temp buffer for extracting a gene row (cache-friendly access)
    std::vector<double> row_buf(n_samples);
    // Pre-compute lambda = x + 0.5 for each sample
    std::vector<double> lambda_buf(n_samples);

    for (int i = 0; i < n_rows; i++) {
        int gi = row_start + i;

        // Extract gene row into contiguous buffer (column-major -> row)
        for (int k = 0; k < n_samples; k++) {
            row_buf[k] = X[gi + k * n_genes];
            lambda_buf[k] = row_buf[k] + 0.5;
        }

        for (int j = 0; j < n_samples; j++) {
            double left_tail = 0.0;
            double y = row_buf[j];
            for (int k = 0; k < n_samples; k++) {
                left_tail += Rf_ppois(y, lambda_buf[k], 1, 0);
            }
            left_tail /= n_samples;
            out[i + j * n_rows] = -1.0 * std::log((1.0 - left_tail) / left_tail);
        }
    }

    return cdf_matrix;
}

// Gaussian kernel CDF for a chunk of gene rows
// Uses precomputed CDF lookup table matching GSVA's kernel_estimation.c
// bw = sd(x) / SIGMA_FACTOR, uses precomputedCdf() for pnorm approximation
// [[Rcpp::export]]
NumericMatrix gene_cdf_gaussian_cpp(NumericMatrix expr, int row_start = 0, int row_end = -1) {
    int n_genes = expr.nrow();
    int n_samples = expr.ncol();
    if (row_end < 0) row_end = n_genes;
    int n_rows = row_end - row_start;

    NumericMatrix cdf_matrix(n_rows, n_samples);
    double* out = REAL(cdf_matrix);
    double* X = REAL(expr);

    // Initialize precomputed CDF table (once)
    if (!is_precomputed) {
        initCdfs();
        is_precomputed = 1;
    }

    // Temp buffer for extracting a gene row
    std::vector<double> row_buf(n_samples);

    for (int i = 0; i < n_rows; i++) {
        int gi = row_start + i;

        // Extract gene row into contiguous buffer for cache efficiency
        for (int k = 0; k < n_samples; k++) {
            row_buf[k] = X[gi + k * n_genes];
        }

        // Compute sd for this gene
        double sd_val = row_sd(row_buf.data(), n_samples);
        double bw = sd_val / SIGMA_FACTOR;
        if (ISNAN(bw) || bw == 0.0) bw = 0.001;

        for (int j = 0; j < n_samples; j++) {
            double left_tail = 0.0;
            double y = row_buf[j];
            for (int k = 0; k < n_samples; k++) {
                left_tail += precomputedCdf(y - row_buf[k], bw);
            }
            left_tail /= n_samples;
            out[i + j * n_rows] = -1.0 * std::log((1.0 - left_tail) / left_tail);
        }
    }

    return cdf_matrix;
}

// ssGSEA for all samples and gene sets (fully self-contained C++)
// expr: genes x samples matrix
// gene_sets_list: list of integer vectors (1-based gene indices)
// alpha: weighting exponent
// Uses average tie-breaking for ranks and stable_sort for gene ordering
// to match R's rank(ties.method="average") and order() exactly.
// [[Rcpp::export]]
NumericMatrix ssgsea_cpp(NumericMatrix expr, List gene_sets_list, double alpha) {
    int n_genes = expr.nrow();
    int n_samples = expr.ncol();
    int n_sets = gene_sets_list.size();
    double total_pos_sum = (double)n_genes * (n_genes + 1) / 2.0;

    NumericMatrix es_matrix(n_sets, n_samples);

    // Pre-extract gene set indices for efficiency
    std::vector<std::vector<int>> gsets(n_sets);
    for (int k = 0; k < n_sets; k++) {
        IntegerVector gs = gene_sets_list[k];
        gsets[k].resize(gs.size());
        for (int i = 0; i < gs.size(); i++) {
            gsets[k][i] = gs[i] - 1;  // convert to 0-based
        }
    }

    // Reusable buffers (avoid per-sample allocation)
    std::vector<std::pair<double, int>> expr_idx(n_genes);
    std::vector<double> ranks_d(n_genes);
    std::vector<double> Ra(n_genes);
    std::vector<int> gene_to_pos(n_genes);

    for (int j = 0; j < n_samples; j++) {
        // Compute ranks for this sample (ascending: rank 1 = lowest)
        // Use average tie-breaking to match R's rank(ties.method="average")
        for (int i = 0; i < n_genes; i++) {
            expr_idx[i] = std::make_pair(expr(i, j), i);
        }
        std::stable_sort(expr_idx.begin(), expr_idx.end());

        // Assign average ranks for ties
        int i = 0;
        while (i < n_genes) {
            int tie_start = i;
            while (i < n_genes - 1 && expr_idx[i].first == expr_idx[i + 1].first) {
                i++;
            }
            double avg_rank = ((double)(tie_start + 1) + (double)(i + 1)) / 2.0;
            for (int t = tie_start; t <= i; t++) {
                ranks_d[expr_idx[t].second] = avg_rank;
            }
            i++;
        }

        // Truncate ranks to integer to match GSVA's mode(R) <- "integer"
        // This is critical: the original GSVA truncates average ranks before
        // computing Ra and gene ordering
        for (int i = 0; i < n_genes; i++) {
            ranks_d[i] = (double)(int)ranks_d[i];
        }

        // Ra = |rank|^alpha (using truncated integer ranks)
        for (int i = 0; i < n_genes; i++) {
            Ra[i] = std::pow(ranks_d[i], alpha);
        }

        // Build gene ranking: order by decreasing rank (stable to match R's order())
        // Reuse expr_idx buffer
        for (int i = 0; i < n_genes; i++) {
            expr_idx[i] = std::make_pair(ranks_d[i], i);
        }
        std::stable_sort(expr_idx.begin(), expr_idx.end(),
                  [](const std::pair<double,int>& a, const std::pair<double,int>& b) {
                      return a.first > b.first;
                  });

        // Build gene_to_pos: gene_index (0-based) -> position (1-based)
        for (int i = 0; i < n_genes; i++) {
            gene_to_pos[expr_idx[i].second] = i + 1;
        }

        // Score each gene set using GSVA's .fastRndWalk formula exactly:
        //   stepCDFinGeneSet = sum(Ra[geneRanking[gSetIdx],j] * (n-gSetIdx+1))
        //                     / sum(Ra[geneRanking[gSetIdx],j])
        //   stepCDFoutGeneSet = (n*(n+1)/2 - sum(n-gSetIdx+1)) / (n-k)
        //   walkStat = stepCDFinGeneSet - stepCDFoutGeneSet
        // where gSetIdx are 1-based positions in decreasing-rank ordering
        // and (n - pos + 1) reverses the position (pos=1 -> weight=n)
        for (int k = 0; k < n_sets; k++) {
            const std::vector<int>& gset = gsets[k];
            int kk = (int)gset.size();
            if (kk == 0) continue;

            double sum_ra_revpos = 0.0;  // sum(Ra[gene] * (n - pos + 1))
            double sum_ra = 0.0;         // sum(Ra[gene])
            double sum_revpos = 0.0;     // sum(n - pos + 1)

            for (int i = 0; i < kk; i++) {
                int gene_0 = gset[i];
                int pos = gene_to_pos[gene_0];       // 1-based position in decreasing order
                double revpos = (double)(n_genes - pos + 1);  // (n - gSetIdx + 1)
                double w_i = Ra[gene_0];              // Ra[geneRanking[pos], j]
                sum_ra_revpos += w_i * revpos;
                sum_ra += w_i;
                sum_revpos += revpos;
            }

            if (sum_ra == 0.0) continue;

            double stepCDFinGeneSet = sum_ra_revpos / sum_ra;
            double stepCDFoutGeneSet = (total_pos_sum - sum_revpos) / (n_genes - kk);
            es_matrix(k, j) = stepCDFinGeneSet - stepCDFoutGeneSet;
        }
    }

    return es_matrix;
}

// ssGSEA scoring from pre-computed R ranks and Ra weights
// R_mat: genes x samples matrix of ranks (from R's rank(ties.method="average"))
// Ra_mat: genes x samples matrix of |rank|^alpha weights
// gene_sets_list: list of integer vectors (1-based gene indices)
// [[Rcpp::export]]
NumericMatrix ssgsea_from_ranks_cpp(NumericMatrix R_mat, NumericMatrix Ra_mat,
                                     List gene_sets_list) {
    int n_genes = R_mat.nrow();
    int n_samples = R_mat.ncol();
    int n_sets = gene_sets_list.size();
    double total_pos_sum = (double)n_genes * (n_genes + 1) / 2.0;

    NumericMatrix es_matrix(n_sets, n_samples);

    for (int j = 0; j < n_samples; j++) {
        // Build gene ranking: order genes by decreasing rank
        // Use stable_sort to match R's order() which preserves original index
        // order for ties (i.e., ties broken by original position)
        std::vector<std::pair<double, int>> rank_idx(n_genes);
        for (int i = 0; i < n_genes; i++) {
            rank_idx[i] = std::make_pair(R_mat(i, j), i);  // 0-based gene index
        }
        std::stable_sort(rank_idx.begin(), rank_idx.end(),
                  [](const std::pair<double,int>& a, const std::pair<double,int>& b) {
                      return a.first > b.first;
                  });

        // Build gene_to_pos lookup: gene_index (0-based) -> position (1-based)
        std::vector<int> gene_to_pos(n_genes);
        for (int i = 0; i < n_genes; i++) {
            gene_to_pos[rank_idx[i].second] = i + 1;
        }

        for (int k = 0; k < n_sets; k++) {
            IntegerVector gset = gene_sets_list[k]; // 1-based gene indices
            int kk = gset.size();
            if (kk == 0) continue;

            double sum_ra_revpos = 0.0;
            double sum_ra = 0.0;
            double sum_revpos = 0.0;

            for (int i = 0; i < kk; i++) {
                int gene_0 = gset[i] - 1;  // convert to 0-based
                int pos = gene_to_pos[gene_0];
                double revpos = (double)(n_genes - pos + 1);
                double w_i = Ra_mat(gene_0, j);
                sum_ra_revpos += w_i * revpos;
                sum_ra += w_i;
                sum_revpos += revpos;
            }

            if (sum_ra == 0.0) continue;

            double stepCDFinGeneSet = sum_ra_revpos / sum_ra;
            double stepCDFoutGeneSet = (total_pos_sum - sum_revpos) / (n_genes - kk);
            es_matrix(k, j) = stepCDFinGeneSet - stepCDFoutGeneSet;
        }
    }

    return es_matrix;
}

// ============================================================================
// GSVA KS walk computation (matches GSVA ks_test.c gsva_rnd_walk)
// ============================================================================

// GSVA random walk for one gene set in one sample
// gset_idx: 1-based gene indices of gene set members
// dec_ord_stat: decreasing order statistic (p - rank + 1) for each gene
// sym_rnk_stat: symmetric rank statistic (|p/2 - rank|) for each gene
// p: number of genes
// tau: weighting parameter
// max_diff: TRUE for max-min scoring
// abs_ranking: used with max_diff
static void gsva_rnd_walk(IntegerVector gset_idx, IntegerVector dec_ord_stat,
                          NumericVector sym_rnk_stat, int p, double tau,
                          bool max_diff, bool abs_ranking,
                          double& score) {
    int k = gset_idx.size();

    // gSetRnk = positions in decreasing order for gene set members
    std::vector<int> gset_rnk(k);
    for (int i = 0; i < k; i++) {
        gset_rnk[i] = dec_ord_stat[gset_idx[i] - 1]; // convert to 0-based
    }

    // Build step_in and step_out vectors
    std::vector<double> step_in(p, 0.0);
    std::vector<int> step_out(p, 1);

    for (int i = 0; i < k; i++) {
        int pos = gset_rnk[i] - 1; // convert to 0-based position
        if (pos >= 0 && pos < p) {
            if (tau == 1.0) {
                step_in[pos] = sym_rnk_stat[gset_idx[i] - 1];
            } else {
                step_in[pos] = std::pow(sym_rnk_stat[gset_idx[i] - 1], tau);
            }
            step_out[pos] = 0;
        }
    }

    // Cumulative sums
    for (int i = 1; i < p; i++) {
        step_in[i] += step_in[i - 1];
        step_out[i] += step_out[i - 1];
    }

    double walk_pos = 0.0, walk_neg = 0.0;

    if (step_in[p - 1] > 0 && step_out[p - 1] > 0) {
        double norm_in = step_in[p - 1];
        double norm_out = (double)step_out[p - 1];

        for (int i = 0; i < p; i++) {
            double wlk = step_in[i] / norm_in - (double)step_out[i] / norm_out;
            if (wlk > walk_pos) walk_pos = wlk;
            if (wlk < walk_neg) walk_neg = wlk;
        }
    }

    // Score extraction (matches C code exactly)
    if (max_diff) {
        if (abs_ranking) {
            score = walk_pos - walk_neg;
        } else {
            score = walk_pos + walk_neg;
        }
    } else {
        score = (walk_pos > std::fabs(walk_neg)) ? walk_pos : walk_neg;
    }
}

// Combined GSVA pipeline: CDF estimation + ranking + scoring in one C++ call
// Avoids intermediate R matrix allocations for CDF and rank matrices
// kcdf: 1 = Poisson, 2 = Gaussian
// [[Rcpp::export]]
NumericMatrix gsva_pipeline_cpp(NumericMatrix expr, List gene_sets_list,
                                 int kcdf, double tau, bool max_diff,
                                 bool abs_ranking) {
    int p = expr.nrow();     // n_genes
    int n = expr.ncol();     // n_samples
    int m = gene_sets_list.size();
    double* X = REAL(expr);

    NumericMatrix es_matrix(m, n);

    // Pre-extract gene rows into contiguous buffers for cache efficiency
    // gene_rows[i][k] = expr[i, k]
    std::vector<std::vector<double>> gene_rows(p, std::vector<double>(n));
    for (int i = 0; i < p; i++) {
        for (int k = 0; k < n; k++) {
            gene_rows[i][k] = X[i + k * p];
        }
    }

    // Pre-compute Poisson lambda or Gaussian bandwidth per gene
    std::vector<double> bw(p, 0.0);
    std::vector<std::vector<double>> lambda_buf;
    if (kcdf == 1) {
        // Poisson: lambda[i][k] = x[i,k] + 0.5
        lambda_buf.resize(p, std::vector<double>(n));
        for (int i = 0; i < p; i++) {
            for (int k = 0; k < n; k++) {
                lambda_buf[i][k] = gene_rows[i][k] + 0.5;
            }
        }
    } else {
        // Gaussian: bw = sd / SIGMA_FACTOR, with precomputed CDF
        if (!is_precomputed) { initCdfs(); is_precomputed = 1; }
        for (int i = 0; i < p; i++) {
            bw[i] = row_sd(gene_rows[i].data(), n) / SIGMA_FACTOR;
            if (ISNAN(bw[i]) || bw[i] == 0.0) bw[i] = 0.001;
        }
    }

    // Reusable per-column buffers
    std::vector<double> cdf_col(p);
    std::vector<std::pair<double, int>> cdf_idx(p);  // for ranking
    std::vector<int> ranks(p);
    IntegerVector dec_ord_stat(p);
    NumericVector sym_rnk_stat(p);

    for (int j = 0; j < n; j++) {
        // Step 1: Compute CDF for column j
        for (int i = 0; i < p; i++) {
            double left_tail = 0.0;
            double y = gene_rows[i][j];
            if (kcdf == 1) {
                for (int k = 0; k < n; k++) {
                    left_tail += Rf_ppois(y, lambda_buf[i][k], 1, 0);
                }
            } else {
                for (int k = 0; k < n; k++) {
                    left_tail += precomputedCdf(y - gene_rows[i][k], bw[i]);
                }
            }
            left_tail /= n;
            cdf_col[i] = -1.0 * std::log((1.0 - left_tail) / left_tail);
        }

        // Step 2: Rank CDF column with ties.method="last"
        // "last" means among tied values, the element with the HIGHER original
        // index gets the HIGHER rank. Sort by (value ASC, index DESC) so that
        // within ties, later indices come first → get lower sorted positions.
        // Then reverse-assign within tie groups.
        // Simpler approach: sort by (value, -index) so tied items with larger
        // original index appear earlier → they get lower sorted positions.
        // Then assign rank = sorted_position gives "first-reversed" = "last".
        for (int i = 0; i < p; i++) {
            // Negate index so larger original indices sort first within ties
            cdf_idx[i] = std::make_pair(cdf_col[i], -i);
        }
        std::stable_sort(cdf_idx.begin(), cdf_idx.end());
        for (int i = 0; i < p; i++) {
            ranks[-cdf_idx[i].second] = i + 1;  // restore original index
        }

        // Step 3: Compute dec_ord_stat and sym_rnk_stat
        for (int i = 0; i < p; i++) {
            int r = ranks[i];
            dec_ord_stat[i] = p - r + 1;
            sym_rnk_stat[i] = std::fabs((double)p / 2.0 - (double)r);
        }

        // Step 4: Score all gene sets
        for (int k = 0; k < m; k++) {
            IntegerVector gset = gene_sets_list[k];
            if (gset.size() == 0) continue;
            double score = 0.0;
            gsva_rnd_walk(gset, dec_ord_stat, sym_rnk_stat, p, tau,
                         max_diff, abs_ranking, score);
            es_matrix(k, j) = score;
        }
    }

    return es_matrix;
}

// Full GSVA scoring: given rank matrix R (from column ranks of CDF values),
// compute enrichment scores for all gene sets and samples
// rank_matrix: genes x samples integer matrix (column ranks of CDF values)
// gene_sets_list: list of integer vectors (1-based gene indices)
// tau, max_diff, abs_ranking: GSVA parameters
// [[Rcpp::export]]
NumericMatrix gsva_score_cpp(IntegerMatrix rank_matrix, List gene_sets_list,
                             double tau, bool max_diff, bool abs_ranking) {
    int p = rank_matrix.nrow();
    int n = rank_matrix.ncol();
    int m = gene_sets_list.size();

    NumericMatrix es_matrix(m, n);

    for (int j = 0; j < n; j++) {
        // ranks2stats for this column
        IntegerVector dec_ord_stat(p);
        NumericVector sym_rnk_stat(p);

        for (int i = 0; i < p; i++) {
            int r = rank_matrix(i, j);
            dec_ord_stat[i] = p - r + 1;
            sym_rnk_stat[i] = std::fabs((double)p / 2.0 - (double)r);
        }

        for (int k = 0; k < m; k++) {
            IntegerVector gset = gene_sets_list[k];
            if (gset.size() == 0) continue;

            double score = 0.0;
            gsva_rnd_walk(gset, dec_ord_stat, sym_rnk_stat, p, tau,
                         max_diff, abs_ranking, score);
            es_matrix(k, j) = score;
        }
    }

    return es_matrix;
}
