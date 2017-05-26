# supporting plotting and diagnostic functions.

returnSamplesOfLengthN <- function(sample_list, N) {
  res <- c()
  for (qq in 1:length(sample_list)) {
    if (length(sample_list[[qq]]) == N) {
      res <- rbind(res, sample_list[[qq]])
    }
  }
  return(res)
}

