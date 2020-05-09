v1 <- read.table("Urn_Peccoud_Ycart.txt", header = TRUE, sep = "\t")
v1[, 1:3] <- apply(v1[, 1:3], 2, function(x) x/sum(x))
v2 <- paste0("Kinetic_Peccoud_Ycart_", 1:3, ".xml") %>% 
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
df1$type <- substr(df1$variable, 1, 1) %>% factor()
df1$id <- substr(df1$variable, 2, 2) %>% paste0("PY-", .) %>% factor()
# ggplot(df1, aes(x = mRNA_copies, y = value, color = variable)) +
#   geom_line() + xlim(0, 20)

df2 <- melt(v2, id.vars = "mRNA_copies")
df2$type <- substr(df2$variable, 1, 1) %>% factor()
df2$id <- substr(df2$variable, 2, 2) %>% paste0("PY-", .) %>% factor()
# ggplot(df2, aes(x = mRNA_copies, y = value, color = variable)) +
#   geom_line() + xlim(0, 20)

df <- rbind(df1, df2)

p1 <-ggplot(df2, aes(x = mRNA_copies, y = value)) +
  facet_grid(. ~ id) + 
  geom_point(data= df1, aes(x = mRNA_copies, 
                            y = value), 
             inherit.aes = FALSE, color = "red", 
             shape = 1, size = 0.5) + 
  geom_line() +
  xlim(0, 15) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6)) +
  xlab(expression(italic(m)[italic(1)])) +
  ylab(expression(italic(widehat(p))[italic(m)[italic(1)]])) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5))

#----------------------------
v1 <- read.table("Urn_Shahrezaei_Swain.txt", header = TRUE, sep = "\t")
v1[, 1:3] <- apply(v1[, 1:3], 2, function(x) x/sum(x))
v2 <- paste0("Kinetic_Shahrezaei_Swain_", 1:3, ".xml") %>% 
  file.path(., "histograms", "hist_3_2.dat") %>% 
  map(., readLines) %>% 
  lapply(., function(x) magrittr::extract(x, i = 3) %>% 
           strsplit(., split = "\t") %>% 
           unlist() %>% 
           as.numeric()) %>% 
  data.frame() %>%
  set_colnames(., paste0("E", 1:3)) %>% 
  head(., n = 600)
v2 <- apply(v2, 2, function(x) x/sum(x)) %>% 
  data.frame()
v2$prot_copies <- 0:(nrow(v2)-1)

df1 <- melt(v1, id.vars = "prot_copies")
df1$type <- substr(df1$variable, 1, 1) %>% factor()
df1$id <- substr(df1$variable, 2, 2) %>% paste0("SS-", .) %>% factor()
# ggplot(df1, aes(x = prot_copies, y = value, color = variable)) +
#   geom_line() + xlim(0, 150) + ylim(0, 0.05)

df2 <- melt(v2, id.vars = "prot_copies")
df2$type <- substr(df2$variable, 1, 1) %>% factor()
df2$id <- substr(df2$variable, 2, 2) %>% paste0("SS-", .)  %>% factor()
# ggplot(df2, aes(x = prot_copies, y = value, color = variable)) +
# geom_line() + xlim(0, 150) + ylim(0, 0.05)

df <- rbind(df1, df2)

p2 <- ggplot(df2, aes(x = prot_copies, y = value)) +
  facet_grid(. ~ id) + 
  geom_point(data= df1, aes(x = prot_copies, 
                            y = value), 
             inherit.aes = FALSE, color = "red", 
             size = 0.5, shape = 1) + 
  geom_line() +
  xlim(0, 150) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6)) +
  xlab(expression(italic(m)[italic(2)])) +
  ylab(expression(italic(widehat(p))[italic(m)[italic(2)]])) +
  scale_y_continuous(breaks = c(0, 0.035, 0.07), limits = c(0, 0.06666667))


#----------------------------
v1 <- read.table("Urn_Leaky_2state.txt", header = TRUE, sep = "\t")
v1[, 1:3] <- apply(v1[, 1:3], 2, function(x) x/sum(x))
v2 <- paste0("Kinetic_Leaky_2state_", 1:3, ".xml") %>% 
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
df1$type <- substr(df1$variable, 1, 1) %>% factor()
df1$id <- substr(df1$variable, 2, 2) %>% paste0("L2S-", .) %>% factor()
# ggplot(df1, aes(x = mRNA_copies, y = value, color = variable)) +
# geom_line() + xlim(0, 20) + ylim(0, 0.5)

df2 <- melt(v2, id.vars = "mRNA_copies")
df2$type <- substr(df2$variable, 1, 1) %>% factor()
df2$id <- substr(df2$variable, 2, 2) %>% paste0("L2S-", .) %>% factor()
# ggplot(df2, aes(x = mRNA_copies, y = value, color = variable)) +
# geom_line() + xlim(0, 20) + ylim(0, 0.5)

df <- rbind(df1, df2)

p3 <- ggplot(df2, aes(x = mRNA_copies, y = value)) +
  facet_grid(. ~ id) + 
  geom_point(data= df1, aes(x = mRNA_copies, 
                            y = value), 
             inherit.aes = FALSE, color = "red", 
             size = 0.5, shape = 1) + 
  geom_line() +
  xlim(0, 15) + 
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6)) +
  xlab(expression(italic(m)[italic(1)])) +
  ylab(expression(italic(widehat(p))[italic(m)[italic(1)]])) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50), limits = c(0, 0.4761905))



#----------------------------
v1 <- read.table("Urn_Cao_Grima.txt", header = TRUE, sep = "\t")
v1[, 1:3] <- apply(v1[, 1:3], 2, function(x) x/sum(x))
v2 <- paste0("Kinetic_Cao_Grima_", 1:3, ".xml") %>% 
  file.path(., "histograms", "hist_4_2.dat") %>% 
  map(., readLines) %>% 
  lapply(., function(x) magrittr::extract(x, i = 3) %>% 
           strsplit(., split = "\t") %>% 
           unlist() %>% 
           as.numeric()) %>% 
  data.frame() %>%
  set_colnames(., paste0("E", 1:3)) %>% 
  head(., n = 600)
v2 <- apply(v2, 2, function(x) x/sum(x)) %>% 
  data.frame()
v2$prot_copies <- 0:(nrow(v2)-1)

df1 <- melt(v1, id.vars = "prot_copies")
df1$type <- substr(df1$variable, 1, 1) %>% factor()
df1$id <- substr(df1$variable, 2, 2) %>% paste0("CG-", .) %>% factor()
# ggplot(df1, aes(x = prot_copies, y = value, color = variable)) +
# geom_line() + xlim(0, 600) + ylim(0, 0.05)

df2 <- melt(v2, id.vars = "prot_copies")
df2$type <- substr(df2$variable, 1, 1) %>% factor()
df2$id <- substr(df2$variable, 2, 2) %>% paste0("CG-", .) %>% factor()
# ggplot(df2, aes(x = prot_copies, y = value, color = variable)) +
# geom_line() + xlim(0, 600) + ylim(0, 0.05)

df <- rbind(df1, df2)

p4 <- ggplot(df2, aes(x = prot_copies, y = value)) +
  facet_grid(. ~ id) + 
  geom_point(data= df1, aes(x = prot_copies, y = value), 
             inherit.aes = FALSE, color = "red", size = 0.5, shape = 1) + 
  geom_line() +
  xlim(0, 600) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 6))  +
  xlab(expression(italic(m)[italic(2)])) +
  ylab(expression(italic(widehat(p))[italic(m)[italic(2)]])) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02), limits = c(0, 0.01904762))

library(cowplot)
p <- plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, labels = "auto")

legend <- ggplot() +
  theme_void() +
  annotate(label = "Urn model", x = -8, y = 0,
           geom = "text") +
  annotate(label = "Kinetic model", x = 26, y = 0, 
           geom = "text") +
  annotate(geom = "line", x = c(7, 17), y = c(0, 0)) +
  annotate(geom = "point", x = c(-20.5, -19.25, 
                                 -18, -16.75, -15.5), y = rep(0, 5),
           color = "red", size = 0.5, shape = 1) +
  coord_cartesian(xlim = c(-40, 50)) +
  annotate(geom = "rect", xmin = -26, xmax = 39,
           ymin = -5, ymax = 5, fill = NA, color = "black") #+
  # theme(panel.border = element_rect(size = c(1,1,1,1),
                                    # fill = NA),
        # plot.margin = margin(0, 0, 0, 0, "cm"))

p_final <- plot_grid(legend, p, ncol = 1,
                     rel_heights = c(0.1, 0.9))
p_final

ggsave("Figure4.pdf",
       plot = p_final,
       width = 16.4,
       height = 8.3, units = "cm")
