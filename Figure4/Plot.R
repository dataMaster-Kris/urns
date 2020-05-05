v1 <- mRNA_obs_counts_in_red
v1[, 1:3] <- apply(v1[, 1:3], 2, function(x) x/sum(x))
v2 <- paste0("Exp_Peccoud_Ycart_", 1:3, ".xml") %>% 
  file.path(., "histograms", "hist_2_2.dat") %>% 
  map(., readLines) %>% 
  lapply(., function(x) magrittr::extract(x, i = 3) %>% 
           strsplit(., split = "\t") %>% 
           unlist() %>% 
           as.numeric()) %>% 
  data.frame() %>%
  set_colnames(., paste0("E", 1:3)) %>% 
  head(., n = 50)
v2 <- apply(v2, 2, function(x) x/sum(x)) %>% 
  data.frame()
v2$mRNA_copies <- 0:(nrow(v2)-1)

df1 <- melt(v1, id.vars = "mRNA_copies")
ggplot(df1, aes(x = mRNA_copies, y = value, color = variable)) +
  geom_line() + xlim(0, 20)

df2 <- melt(v2, id.vars = "mRNA_copies")
ggplot(df2, aes(x = mRNA_copies, y = value, color = variable)) +
  geom_line() + xlim(0, 20)
