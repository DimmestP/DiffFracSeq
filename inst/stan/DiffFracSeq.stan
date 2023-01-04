//
// This Stan program defines the two fraction DiffFracSeq model
// It expects three integer count arrays: tot_obs, sup_obs, and pel_obs
// and determines latent counts for pellet and supernatant samples together with
// a linear sum of terms deducing changes total expression and in fractionation
// across conditions.

data {
  // Number of RNA species (i.e. genes)
  int<lower=1> NRNA;

  // Number of replicates
  int<lower=1> NREP;

  // Number of conditions
  int<lower=1> NCON;

  // Counts from the two fractions and the total sample
  int<lower=0> tot_obs[NCON, NREP, NRNA];
  int<lower=0> sup_obs[NCON, NREP, NRNA];
  int<lower=0> pel_obs[NCON, NREP, NRNA];
}
parameters {
  // Normalising factors
  real scale_factor_mean;
  real<lower=0> tot_scale_factor[NCON, NREP];
  real sup_scale_factor[NCON, NREP];
  real pel_scale_factor[NCON, NREP];

  // latent count linear terms
  vector[NRNA] gene_count_condition[NCON];
  vector[NRNA] pel_count_condition[NCON];

  // latent count prior parameters
  real norm_alpha;
  real<lower=0> norm_beta;

  // dispersion parameter for counts
  real<lower=0> phi[3];
}

transformed parameters{
  // latent counts
  vector[NRNA] sup_latent[NCON];
  vector[NRNA] pel_latent[NCON];

  for(con in 1:NCON){
    sup_latent[con] = gene_count_condition[con];

    pel_latent[con] = sup_latent[con] + pel_count_condition[con];
  }
}

model{
  // priors on gene_count_condition paramters
  norm_alpha ~ normal(7,2);
  norm_beta ~ normal(2,1);

  for(con in 1:NCON){

    // scale factor priors
    scale_factor_mean ~ normal(0, 0.5);
    tot_scale_factor[con] ~ normal(10, 0.1);
    pel_scale_factor[con] ~ normal(scale_factor_mean, 0.1);
    sup_scale_factor[con] ~ normal(scale_factor_mean, 0.1);

   // neg bin overdispersion prior
    phi ~ normal(100, 10);

    // latent counts linear terms priors
    gene_count_condition[con] ~ normal(norm_alpha, norm_beta);
    pel_count_condition[con] ~ normal(0, 1);


    for(rep in 1:NREP){
      // observed counts from fractions
      sup_obs[con, rep] ~ neg_binomial_2_log(sup_scale_factor[con, rep] + sup_latent[con],
                                              phi[2]);

      pel_obs[con, rep] ~ neg_binomial_2_log(pel_scale_factor[con, rep] + pel_latent[con],
                                              phi[3]);

      // observed counts from total
      tot_obs[con, rep] ~ neg_binomial_2(tot_scale_factor[con, rep] * (exp(pel_latent[con]) + exp(sup_latent[con])),
					 	                              phi[1]);
    }
  }
}
