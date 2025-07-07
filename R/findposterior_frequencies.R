# Function to find posterior frequencies for a given locus index

findposteriorfrequencies = function(locus_index, tempdata, maxMOI, frequencies_RR) {
  data = tempdata[,c(1:maxMOI)+(locus_index-1)*maxMOI];
  nalleles = frequencies_RR[[1]][locus_index]
  freq_prior_alpha = rep(1,nalleles);
  freq_posterior_alpha = freq_prior_alpha + table(factor(c(data),levels=c(1:nalleles)));
  frequencies_RR[[2]][locus_index, 1:nalleles] <- gtools::rdirichlet(1, freq_posterior_alpha);
  
  return(frequencies_RR)
}
