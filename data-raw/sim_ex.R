# Regenerates data/sim_ex.rda as a lightweight subset of the original
# simulated dataset, for fast test iteration.
#
# The original (20 chromosomes, ~1,500 markers each, 30,732 markers total)
# made the test suite slow: compute_LD_decay()/ld_complexity_reduction()/
# merge_ld_clusters() all scale with marker and chromosome count, and the
# full test suite took ~162s. This subset (3 chromosomes, 300 markers each,
# 900 total) runs the same suite in ~4.5s with identical pass/fail results,
# since no test asserts an exact count tied to the original size -- they're
# all self-referential (e.g. n_snps == lengths(members)) or computed from
# sim_ex itself (e.g. nrow(decay_sum) == length(unique(map$Chr))).
#
# Markers are taken as a CONTIGUOUS window per chromosome, not a random
# subsample -- LD structure depends on physical proximity, and randomly
# scattered markers would destroy it (huge artificial gaps, no real blocks
# left to cluster).
#
# The original full-size sim_ex.rda is recoverable from git history (last
# committed at 0a274e7, before this subsetting), if a larger/different
# fixture is ever needed:
#   git show 0a274e7:data/sim_ex.rda > data/sim_ex.rda

library(data.table)

## load the original, pre-subset sim_ex -- restore it first if this is
## being re-run from a state where data/sim_ex.rda is already the small
## version, e.g. via the git show command above
sim_ex_full <- local({
  e <- new.env()
  load("data/sim_ex.rda", envir = e)
  e$sim_ex
})

chrs      <- c("Chr1", "Chr2", "Chr3")
n_per_chr <- 300

map_full <- as.data.table(sim_ex_full$map)
idx <- unlist(lapply(chrs, function(ch) {
  which(map_full$Chr == ch)[seq_len(min(n_per_chr, sum(map_full$Chr == ch)))]
}))

sim_ex <- list(
  map = sim_ex_full$map[idx, ],
  GTs = sim_ex_full$GTs[, idx]
)

save(sim_ex, file = "data/sim_ex.rda", compress = "xz")
