insert_sorted_cpp <- function(sorted_x, value) {
    .Call(`_RepsampleR_insert_sorted_cpp`, sorted_x, value)
}

ks2_d_with_candidate_cpp <- function(sample_sorted, candidate, pop_sorted) {
    .Call(`_RepsampleR_ks2_d_with_candidate_cpp`, sample_sorted, candidate, pop_sorted)
}

ks1_norm_d_with_candidate_cpp <- function(sample_sorted, candidate, mean, sd) {
    .Call(`_RepsampleR_ks1_norm_d_with_candidate_cpp`, sample_sorted, candidate, mean, sd)
}

ks1_cdf_d_with_candidate_cpp <- function(sample_sorted, sample_cdf_sorted, candidate, candidate_cdf) {
    .Call(`_RepsampleR_ks1_cdf_d_with_candidate_cpp`, sample_sorted, sample_cdf_sorted, candidate, candidate_cdf)
}
