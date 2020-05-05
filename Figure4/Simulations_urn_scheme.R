setwd("~/Downloads/temp_BJ_revision/Final_push/Push2/Figure4")
library(tidyverse)
library(reshape2)
library(magrittr)

mrna_params <- data.frame(v0 = c(0.02, 0.005, 0.02),
                          k0 = c(0.03, 0.03, 0.003),
                          k1 = c(0.01, 0.01, 0.001),
                          d0 = 0.005)
mrna_params$nRed <- mrna_params %$% divide_by(k0, d0)
mrna_params$nTot <- (mrna_params$k0 + mrna_params$k1)/mrna_params$d0
mrna_params$nGrey <- mrna_params$nTot - mrna_params$nRed
mrna_params$mu <- mrna_params %$% divide_by(v0, d0)

prot_params <- data.frame(v0 = 0.02,
                          v1 = 0.01, 
                          k0 = c(0.003, 0.0026),
                          k1 = c(0.001, 0.0014), 
                          d0 = 0.005, 
                          d1 = 0.0005)
prot_params$nRed <- 0.5*(prot_params$v0 + prot_params$k0 + prot_params$k1 -
                           sqrt((prot_params$v0 + prot_params$k0 + prot_params$k1)^2 -
                                  4*prot_params$v0*prot_params$k0))/prot_params$d1
prot_params$nTot <- (prot_params$k0 + prot_params$k1)/prot_params$d1
prot_params$nGrey <- prot_params$nTot - prot_params$nRed
prot_params$alpha <-  0.5*(prot_params$v0 + prot_params$k0 + prot_params$k1 +
                             sqrt((prot_params$v0 + prot_params$k0 + prot_params$k1)^2 -
                                    4*prot_params$v0*prot_params$k0))/prot_params$d1
prot_params$beta <- prot_params %$% divide_by(v1, d0)

n_sims <- 10^5

#-----------------------------
#Function to perform assignment to recipient urns and count the balls assigned to red urns.
assign_to_recepient_urns_and_count <- function(nBlack, nRed, nGrey) {
  map(nBlack, function(x) {
    
    #Assign to red and grey urns.
    permute_balls_and_dividers <- c(rep("R", nRed), rep("G", nGrey), rep(".", x)) %>%
      sample()
    
    #To simplify counting the balls assigned to red urns ...
    #... in cases where the last element in permutation is a ball 
    #... move the balls to the front.
    while (tail(permute_balls_and_dividers , 1) == ".") {
      permute_balls_and_dividers <-
        c(tail(permute_balls_and_dividers , 1),
          permute_balls_and_dividers[-length(permute_balls_and_dividers)])
    }
    
    #Count the number of balls assigned to red urns.
    rle_this <- rle(permute_balls_and_dividers)
    rle_this <- data.frame(cnt = rle_this$lengths, 
                           element = rle_this$values)
    which_ball <- which(rle_this$element == ".")
    which_assigned_to_red <- rle_this$element[which_ball + 1] == "R"
    cntRed <- rle_this$cnt[which_ball[which_assigned_to_red]] %>% sum()
    
  }) %>% unlist() %>% table()
}

#--------------
#Implement urn scheme for Peccoud-Ycart model.
mRNA_obs_counts_in_red <- list()
set.seed(1234)
for (i in 1:nrow(mrna_params)) {
  print(i)
  
  #Sampling from the master urn.
  nBlack <- rpois(n_sims, mrna_params$mu[i])
  nRed <- mrna_params$nRed[i] %>% round()
  nGrey <- mrna_params$nGrey[i] %>% round()
  
  mRNA_obs_counts_in_red[[i]] <- assign_to_recepient_urns_and_count(nBlack, nRed, nGrey)
  
}

mRNA_obs_counts_in_red <- lapply(mRNA_obs_counts_in_red, function(x) {
  y <- numeric(50)
  y[x %>% names() %>% as.numeric() %>% add(., 1)] <- x %>% as.numeric()
  y
}) %>% data.frame() %>% set_colnames(., NULL) %>% 
  cbind(., mRNA_copies = 0:49)

write.table(mRNA_obs_counts_in_red,
            file = paste0("mrna_obs_via_urn_scheme_n_sims_", 
                          n_sims, ".txt"),
            quote = F,
            row.names = F,
            sep = "\t")


#--------------
#Implement urn scheme for Shahrazaei-Swain model.

prot_obs_counts_in_red <- list()
set.seed(1234)
for (i in 1:2) {
  print(i)
  
  #Sampling from the master urn.
  nBlack <- rnbinom(n = n_sims, 
                    size = prot_params$alpha[i],
                    prob = 1/(1+prot_params$beta[i]))
  nRed <- prot_params$nRed[i] %>% round()
  nGrey <- prot_params$nGrey[i] %>% round()
  
  prot_obs_counts_in_red[[i]] <- assign_to_recepient_urns_and_count(nBlack, nRed, nGrey)
  
}

prot_obs_counts_in_red <- lapply(prot_obs_counts_in_red, function(x) {
  y <- numeric(300)
  y[x %>% names() %>% as.numeric() %>% add(., 1)] <- x %>% as.numeric()
  y
}) %>% data.frame() %>% set_colnames(., NULL) %>% 
  cbind(., prot_copies = 0:299)

write.table(prot_obs_counts_in_red,
            file = paste0("prot_obs_via_urn_scheme_n_sims_", 
                          n_sims, ".txt"),
            quote = F,
            row.names = F,
            sep = "\t")

