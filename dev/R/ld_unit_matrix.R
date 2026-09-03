#' One variable per LD-complexity-reduction cluster, for your own scan
#'
#' [ld_outlier_test()] tests clusters, not markers, and a cluster is not itself
#' something you can hand to an association test -- it has to become ONE variable
#' per individual first. This builds that variable, four ways, and stops there:
#' running the test on it (EMMAX, LFMM, a GLM, whatever) is your job, exactly as
#' with the per-marker p-values [ld_scan()] expects. The package supplies LD
#' structure and region assembly; it does not run your scan.
#'
#' The four representations answer different questions and are not
#' interchangeable in general, only observed to be close on some panels:
#'
#' - `"consensus_dosage"`: polarise every member to a common allele, then
#'   row-mean. No size penalty -- a larger cluster is a BETTER-ESTIMATED
#'   consensus, so size helps rather than hurts. Averages away signal that is
#'   genuinely SNP-specific (differing even between strongly linked markers).
#' - `"eMLG"`: [make_eMLGs()]'s own block consensus. Conceptually close to
#'   `consensus_dosage` but not verified identical; included so that question is
#'   answered empirically rather than assumed.
#' - `"representative"`: the cluster's core SNP, chosen by [ld_complexity_reduction()]
#'   for CENTRALITY (highest median r^2 to the rest of the cluster) -- not for
#'   signal. A real, single SNP, at whatever centrality picked it.
#' - `"best_snp"`: [eMLG_best_snp()] -- the member SNP most correlated with the
#'   block consensus, missing calls filled from the consensus. A real SNP, but
#'   chosen by agreement with the block rather than by structural centrality.
#'
#' @param GTs Genotype dosage matrix (individuals x SNPs, column names = markers).
#' @param stage1 An `ld_complexity_reduction` object (or equivalent list with
#'   `clusters$members`, `clusters$core_snp`).
#' @param size_floor Minimum markers per cluster to be included (default 8L).
#'   `"eMLG"`'s own default inside [make_eMLGs()] is `l_min = 10`, LARGER than a
#'   typical `size_floor` -- pass `emlg_args$l_min <= size_floor` or clusters will
#'   be silently dropped there that survive here.
#' @param repr One of `"consensus_dosage"`, `"eMLG"`, `"representative"`,
#'   `"best_snp"`.
#' @param emlg_args Extra arguments to [make_eMLGs()] (used only when
#'   `repr = "eMLG"`).
#' @param best_snp_args Extra arguments to [eMLG_best_snp()] (used only when
#'   `repr = "best_snp"`); note it is a post-processing step over
#'   [ld_prune_and_eMLG()]'s result, not over `stage1` directly, so this arm needs
#'   that result supplied via `prune_result` rather than being derivable from
#'   `stage1` alone.
#' @param prune_result An [ld_prune_and_eMLG()] result, required only when
#'   `repr = "best_snp"`.
#'
#' @return For `repr %in% c("consensus_dosage", "eMLG")`: an individuals x units
#'   numeric matrix. For `repr = "representative"`: a character vector of marker
#'   names (one per unit; index into `GTs` yourself). For `repr = "best_snp"`: an
#'   individuals x units matrix (each column is an observed SNP with gaps filled
#'   from its block's consensus). Every return carries `attr(, "units")`, a
#'   data.table of `unit_id`, `Chr`, `n_markers`, aligned to columns.
#'
#' @seealso [ld_outlier_test()], [ld_complexity_reduction()], [make_eMLGs()],
#'   [eMLG_best_snp()]
#' @export
ld_unit_matrix <- function(GTs, stage1, size_floor = 8L,
                           repr = c("consensus_dosage", "eMLG", "representative", "best_snp"),
                           emlg_args = list(), best_snp_args = list(),
                           prune_result = NULL) {
  repr <- match.arg(repr)
  stop("ld_unit_matrix(): not yet implemented (signature under review)")
}
