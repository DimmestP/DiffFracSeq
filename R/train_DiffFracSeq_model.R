#' Creates the input data for the stan model
#'
#' Creates the input list object from a tibble of transcripts counts
#' required for stan to train the DiffFracSeq model
#'
#' @export
#' @param count_data A data frame of transcript counts with columns: total,
#' sup, pel, gene_name, rep, and condition.
#' @return A list with values NRNA, NREP, NCON, tot_obs, sup_obs, and pel_obs.
#'
#' @importFrom tidyr %>%
#'
#' @examples
create_stan_data <- function(count_data) {
  # fix order of counts so grouped together by replicate
  count_data = count_data %>%
    dplyr::arrange(gene_name, rep, condition)

  # define number of conditions, rep and RNA
  NRNA=length(unique(count_data$gene_name))
  NREP=length(unique(count_data$rep))
  NCON=length(unique(count_data$condition))

  ## make data for stan fit
  list(NRNA=NRNA,
       NREP=NREP,
       NCON=NCON,
       tot_obs=array(as.integer(round(count_data$total)), c(NCON, NREP, NRNA)),
       sup_obs=array(as.integer(round(count_data$sup)), c(NCON, NREP, NRNA)),
       pel_obs=array(as.integer(round(count_data$pel)), c(NCON, NREP, NRNA))
  )

}

#' Trains DiffFracSeq model using the input data
#'
#' First calls the create_stan_data function to create the stan input list object
#' then calls the rstan::sampling function to learn the model parameters using
#' the data set
#'
#' @export
#' @param count_data A data frame of transcript counts with columns: total,
#' sup, pel, gene_name, rep, and condition.
#' @param nchain the number of chains to run when fittting the model
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @examples
train_DiffFracSeq_model <- function(count_data, nchain = 4, ...) {
  # Create stan data
  stan_data <- create_stan_data(count_data)

  # Train DiffFracSeq
  stanfit <- rstan::sampling(stanmodels$DiffFracSeq, data=stan_data, chains = nchain, ...)
  return(stanfit)
}
#' Compares fraction counts between conditions
#'
#' Determines if the ratio of counts between two fractions changes significantlu
#' across all permutations of conditioms
#'
#' @export
#' @param condition A vector of equal size to ratio stating which condition that ratio was taken from.
#' @param ratio A vector of equal size to condition stating the ratio of counts between the two fractions
#'
#' @return A tibble of pvalues comparing rations between all pairs of conditions with three columns:
#' condition_A, condition_B, and p.value.
#'
#' @examples
test_frac_across_condition <- function(condition, ratio){
  unique_permutations <- tibble(condition_A = unique(condition), condition_B = unique(condition)) %>%
    expand(condition_A, condition_B) %>%
    filter(condition_A > condition_B)

  unique_permutations %>%
    dplyr::select(condition_A) %>%
    inner_join(tibble(condition_A = condition, ratio_A = ratio), by = "condition_A") %>%
    bind_cols(unique_permutations %>%
                dplyr::select(condition_B) %>%
                inner_join(tibble(condition_B = condition, ratio_B = ratio), by = "condition_B")) %>%
    group_by(condition_A, condition_B) %>%
    summarise(p.value = sum(ratio_A > ratio_B) / (length(ratio)/ length(unique(condition))), .groups = "drop")
}
