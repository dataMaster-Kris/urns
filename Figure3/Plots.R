setwd("~/Downloads/temp_BJ_revision/Final_push/Push2/Figure3")
library(tidyverse)
library(reshape2)
library(magrittr)
library(cowplot)

#-------------------
#Coefficients of generating functions using Mathematica

files <- list.files() %>% 
  subset(., endsWith(., ".csv"))

dat <-  map(files, read.csv, header = FALSE) %>% 
  as.data.frame() %>%
  set_colnames(files)

dat$m <- 0:30

df <- melt(dat, id.vars = "m")
df$type <- df$variable %>% 
  as.character() %>%
  strsplit(., split = "_") %>%
  lapply(., function(x) x[1]) %>% 
  unlist()

df$set <- df$variable %>% 
  as.character() %>%
  strsplit(., split = "_") %>%
  lapply(., tail, n = 1) %>% 
  unlist() %>%
  str_remove_all(., ".csv")

labels <- c("italic(k)[0]~'/'~italic(k)[1]~'='~100",
            "italic(k)[0]~'/'~italic(k)[1]~'='~10")
df$set <- labels[(df$set == 10) + 1]

df1 <- subset(df, type != "Hypergeometric")
df1$type <- factor(df1$type, levels = c("Delaporte", "NB"))
df2 <- subset(df, type == "Hypergeometric")



p1 <- ggplot(df1, aes(x = m, y = value, shape = type)) +
  geom_point(data = df1, aes(x = m, y = value, shape = type)) +
  geom_line(data = df2, aes(x=m, y = value)) +
  facet_grid(set ~ ., labeller = label_parsed) +
  coord_cartesian(xlim = c(0, 20)) +#, ylim = c(-6, 0)) +
  labs(x= expression(italic(m)[1]), y= expression(italic(p)[italic(m[1])])) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        strip.text = element_text(size=7.5, lineheight=10.0),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_shape_manual(values = c(0,1, 2))

p_legend <- ggplot() +
  theme_void()+
  annotate(x = -4.7, xend = -4, y = 0, yend = 0, geom = "segment")+
  annotate(x = -2.3, y = 0, 
           label = "Exact solution", geom = "text",
           size = 2.5)+
  annotate(x = 0.2, y = 0, geom = "point", shape = 0) +
  annotate(x = 3, y = 0, label = "Transcriptional lapse", 
           geom = "text", size = 2.5)+
  annotate(x = 5.8, y = 0, geom = "point", shape = 2) +
  annotate(x = 8.5, y = 0, label = "Transcriptional burst", 
           geom = "text", size = 2.5)+
  coord_cartesian(ylim = c(-0.5, 0.5), xlim = c(-5, 10))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


p = plot_grid(p_legend, p1, ncol = 1, labels = c(NA, NA), align = "hv",
              nrow = 2, rel_heights = c(0.1, 0.9))
p

#----------------
#Fig 3a.
s1 <- data.frame(y=c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
                 x=c(0, 10, 10, 11, 11, 33, 33, 33.5, 33.5, 50),
                 type = "Transcriptional\nburst")

s2 <- data.frame(y=c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1), 
                 #c("OFF", "ON")[c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1)+1], 
                 x=c(0, 15, 15, 16.5, 16.5, 24, 24, 25, 25, 50),
                 type = "Transcriptional\nlapse")

df <- rbind(s1, s2)

pa <- ggplot(df, aes(x = x, y = y)) +
  geom_path() +
  facet_grid(type ~ .) +
  scale_y_continuous(breaks = c(0, 1)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  xlab("Time (arbitrary units)") +
  ylab("Transcription rate")

p2 <- plot_grid(pa, p, labels = "auto", scale =0.95)

ggsave("Figure3.pdf", 
       p2,
       width = 17.4,
       height = 8.7,
       units = "cm")
