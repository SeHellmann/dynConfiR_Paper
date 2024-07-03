## Draw proportion of confidence judgments from dirichlet dist
gen_rating_proportions <- function(n, alpha=NULL, min_prop=.02) {
  ## See https://github.com/dkahle/dirichlet
  rdirichlet <- function(n, alpha) {
    normalize <- function(.) . / sum(.)
    samps <- vapply(alpha, function(al) rgamma(n, al, 1), numeric(n))
    if (n == 1) {
      matrix(normalize(samps), nrow = 1, ncol = length(samps))
    } else {
      t(apply(samps, 1, normalize))
    }  
  }
  if (is.null(alpha)) load("helper_fcts/collected_rating_proportions_and_alpha.RData")
  rating_probs <- rdirichlet(n, alpha=as.numeric(alpha))
  while (any(rating_probs< min_prop)) {
    rating_probs <- rating_probs[apply(rating_probs, 1, min)>.02,]
    print(paste("Sample again", n-nrow(rating_probs), "propotion vectors"))
    rating_probs <- rbind(rating_probs,
                rdirichlet(n=n-nrow(rating_probs), alpha=as.numeric(alpha)))
  }
  colnames(rating_probs) <- paste0("prob_rating", 1:length(alpha))
  return(rating_probs)
}
