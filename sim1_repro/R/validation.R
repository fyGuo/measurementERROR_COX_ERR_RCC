#### Internal validation study sampling ####

# sample n_subj main-study subjects, then n_pair (subject, year) pairs from the
# years they were under follow-up, with at most max_per_id pairs per subject
sample_validation <- function(time, n_subj, n_pair, max_per_id, T) {
  vid <- sample.int(length(time), n_subj)
  nyr <- pmin(ceiling(time[vid]), T)
  pool_id <- rep(vid, nyr)
  pool_l <- sequence(nyr)
  o <- sample.int(length(pool_id))
  pool_id <- pool_id[o]; pool_l <- pool_l[o]
  rank_in_id <- ave(seq_along(pool_id), pool_id, FUN = seq_along)
  keep <- which(rank_in_id <= max_per_id)
  if (length(keep) < n_pair) stop("not enough eligible validation pairs")
  keep <- keep[1:n_pair]
  data.frame(id = pool_id[keep], l = pool_l[keep])
}
