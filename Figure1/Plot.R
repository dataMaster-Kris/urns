library(tidyverse)
library(magrittr)
library(cowplot)
library(seewave)

# set.seed(123234)
# params <- data.frame(k0 = 0.05,
#                      k1 = 0.2,
#                      d0 = 0.05,
#                      v0 = 0.00005)
# 
# count1 <- c()
# count2 <- c()
# 
# n_runs <- 10000
# 
# pdf("Samples_take1.pdf")
# for (i in 1:100) {
#   if (!(i %% 1000)) print(i)
#   time_intrvls <- params %$%
#     c(k0 = k0, k1 = k1) %>%
#     rep(., n_runs) %>% 
#     map_dbl(~ rexp(1, .x)) 
#   
#   polymrase_arrvl <- params %$%
#     rpois(1, v0*(time_intrvls %>% sum())) %>%
#     seq(to = .) %>%
#     map_dbl(., ~ rexp(1, params$v0)) %>%
#     cumsum() %>%
#     data.frame(t = ., y = 1)
#   
#   df <- data.frame(Stop = time_intrvls %>% cumsum()) 
#   df$Start <- c(0, df$Stop %>% extract(., -length(.)))
#   df$type <- c("OFF", "ON") %>% rep(., n_runs) %>%
#     factor(., levels = c("ON", "OFF"))
#   
#   polymrase_arrvl %<>% subset(., t < time_intrvls %>% sum())
#   polymrase_arrvl$type <- apply(polymrase_arrvl, 1,
#                                 function(x) {
#                                   subset(df, (Start < x[1]) & (x[1] < Stop))$type %>% as.character()
#                                   })
#   
#   time_scale <- params %$%
#     sum(k0, k1) %>%
#     divide_by(1, .)
#   
#   
#   count1 <- c(polymrase_arrvl$type %>% 
#                 equals(., "ON") %>%
#                 sum(), 
#               count1)
#   
#   plot_layers <- list(coord_cartesian(xlim = c(0, time_intrvls %>% sum()),
#                                       ylim = c(0, 1)),
#                       scale_x_continuous(breaks = time_intrvls %>% 
#                                            cumsum() %>% 
#                                            as.numeric() %>%
#                                            c(0, .)),
#                       theme(axis.line.x = element_line(),
#                             panel.grid = element_blank(),
#                             axis.text = element_blank(),
#                             axis.ticks.y = element_blank(),
#                             legend.title = element_blank()
#                             # legend.direction = "horizontal",
#                             # legend.position = "bottom"
#                       ),
#                       # geom_segment(data = df,
#                       #              aes(x = Start, xend = Stop, color = type),
#                       #              y = 0, yend = 0, size = 2),
#                       geom_segment(aes(x = t, y = y, xend =t, color = type), yend = 0),
#                       geom_point(),
#                       ylab(NULL),
#                       xlab("Time"))
#   
#   
#   p1 <- ggplot(polymrase_arrvl, aes(t, y))  + 
#     plot_layers +
#     ggtitle(paste0("Transcriptions : ", count1[1]))
#   
#   nudge_by <- runif(nrow(polymrase_arrvl), min = -100*time_scale, max = 100*time_scale)
#   polymrase_arrvl$t %<>% add(., nudge_by)
#   polymrase_arrvl$type <- apply(polymrase_arrvl, 1,
#                                 function(x) subset(df, (Start < as.numeric(x[1])) & 
#                                                      (as.numeric(x[1]) < Stop))$type)
#   
#   count2 <- c(polymrase_arrvl$type %>% 
#                 equals(., "ON") %>%
#                 sum(), 
#               count2)
#   
#   p2 <- ggplot(polymrase_arrvl, aes(t, y))  +
#     plot_layers +
#     ggtitle(paste0("Transcriptions : ", count2[1]))
#   
#   p <- plot_grid(p1, p2, nrow =2)
#   
#   print(p)
# }
# 
# dev.off()


df <- data.frame(x = seq(0, 2*pi, length.out = 1000) %>% sin(),
                 y = seq(0, 2*pi, length.out = 1000) %>% cos())

df1 <- data.frame(x = (runif(10^2.4, 0, 1)*runif(10^2.4, 0, 4*pi)) %>% sin(),
                  y = (runif(10^2.4, 0, 1)*runif(10^2.4, 0, 4*pi)) %>% cos()) %>%
  subset(., x^2 + y^2 < 0.9) %>%
  subset(., (x+0.5)^2 + y^2 > 0.0625)

df1$group <- c("Birth factor", "Co-solute")[sample(2, nrow(df1), replace = T)] %>% 
  factor(., levels = c("Birth factor", "Co-solute"))

df2 <- data.frame(x = seq(0, 2*pi, length.out = 1000) %>% sin() %>%
                    multiply_by(., 0.25) %>% add(., -0.5),
                  y = seq(0, 2*pi, length.out = 1000) %>% cos() %>%
                    multiply_by(., 0.25))

p1 <- ggplot(df, aes(x, y)) +
  geom_path() +
  theme_minimal() +
  geom_point(data = df1, aes(x, y, shape = group), inherit.aes = F) +
  geom_line(data = df2, aes(x, y)) +
  ggtitle("Cytosol", subtitle = "(The master urn)") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "top") +
  scale_shape_manual(values = c(19, 1)) +
  coord_cartesian(xlim = c(-1.5, 1.5))


#--------------

df <- data.frame(x = NA, y = NA)
p_basic <- ggplot(df, aes(x, y)) +
  theme_minimal() +
  annotate(geom = "segment",
           x = 0,
           xend = 1, 
           y = 0, 
           yend = 0) +
  annotate(geom = "segment",
           x = 0.33,
           xend = 0.33, 
           y = 0, 
           yend = 0.5) +
  annotate(geom = "segment",
           x = 0.33,
           xend = 1, 
           y = 0.5, 
           yend = 0.5,
           arrow = arrow(length = unit(0.07, "inches"))) +
  theme(#panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank())

p2a <- p_basic + annotate(geom = "text",
                          x = 0.57,
                          y = 0.5,
                          label = "X", size = 3) +
  annotate(geom = "text",
           x = 0.67,
           y = 0.67,
           label = "OFF", size = 3) +
  coord_cartesian(ylim = c(-0.1, 0.8))

p2c <- ggplot(df, aes(x, y)) +
  annotate(geom = "segment",
           x = 0,
           xend = 1, 
           y = 0.15, 
           yend = 0.15,
           arrow = arrow(length = unit(0.05, "inches"))) +
  annotate(geom = "segment",
           x = 1,
           xend = 0, 
           y = -0.15, 
           yend = -0.15,
           arrow = arrow(length = unit(0.05, "inches"))) +
  annotate(geom = "text",
           x = 0.5, 
           y = 0.67,
           label = expression(italic("k")[0])) +
  annotate(geom = "text",
           x = 0.5, 
           y = -0.67,
           label = expression(italic("k")[1])) +
  theme_nothing() +
  coord_cartesian(ylim = c(-3, 3))

p2_prefix <- ggplot(data.frame(), aes()) +
  # theme_nothing() +
  ggtitle("Time intervals", 
          subtitle = "(The recipient urns)") +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

p2x <- plot_grid(p2a, 
                 p2c, 
                 p_basic + 
                   annotate(geom = "text",
                            x = 0.67,
                            y = 0.67,
                            label = "ON",
                            size = 3) +
                   coord_cartesian(ylim = c(-0.1, 0.8)), nrow = 1)

p_interval_legend <- ggplot() +
  annotate(geom = "rect", xmin = 0, ymin = -0.325, xmax = 0.6, ymax = -0.175,
           color = "black", fill = "#e10600") +
  coord_cartesian(xlim = c(-1, 4), ylim = c(-1, 1)) +
  annotate(geom = "rect", xmin = 0, ymin = 0.7, xmax = 0.6, ymax = 0.85,
           color = "black", fill = "grey") +
  annotate(geom = "text", x = 2.05, y = -0.25, label = "ON interval", size = 3) +
  annotate(geom = "text", x = 2.15, y = 0.72, label = "OFF interval", size = 3) +
  theme_nothing()+
  theme(axis.title = element_blank(),
        axis.text = element_blank()) +
  scale_y_continuous(expand = c(0.1, 0.1))

p2 <- plot_grid(p2_prefix, p_interval_legend,
                p2x, nrow = 3, rel_heights = c(0.25, 0.26, 0.49))

#--------

df3 <- data.frame(x = -2:9 + 0.5, y = 0)
p3 <- ggplot(df3, aes(x, y)) +
  geom_point() +
  coord_cartesian(xlim = c(-6, 14)) +
  theme_nothing()

#-----------------


set.seed(12345)
params <- data.frame(k0 = 0.3,
                     k1 = 0.2,
                     d0 = 0.05,
                     v0 = 1)#0.007)

for (ix in 1:3) {
  count1 <- c()
  # n_runs <- 100
  # 
  # time_intrvls <- params %$%
  #   c(k0 = k0, k1 = k1) %>%
  #   rep(., n_runs) %>%
  #   map_dbl(~ rexp(1, .x))
  time_intrvls <- rexp(1, params$k0)
  while (time_intrvls %>% max() %>% is_less_than(1/params$d0)) {
    k_this <- ifelse(length(time_intrvls) %% 2, params$k1, params$k0)
    time_intrvls <- c(time_intrvls, tail(time_intrvls, 1) +
                        rexp(1, k_this))
  } 
  
  time_intrvls[length(time_intrvls)] <- 1/params$d0
  
  polymrase_arrvl <-  rpois(1, params$v0/params$d0) %>% #params %$%
    #rpois(1, v0*(time_intrvls %>% sum())) %>%
    seq(to = .) %>%
    map_dbl(., ~ rexp(1, params$v0)) %>%
    cumsum() %>%
    data.frame(t = ., y = 0.9)
  
  # ggplot(polymrase_arrvl)  + 
  #   geom_point(aes(x=t, y =y)) +
  #   # geom_segment(aes(x = Start, xend = Stop, color = type),
  #   #              y = 0, yend = 0, size = 2) +
  #   theme_nothing() +
  #   coord_cartesian(ylim = c(0, 0)) + 
  #   # scale_x_continuous(expand = c(0, 0)) +
  #   scale_y_continuous(expand = c(0, 0)) +
  #   scale_color_manual(values = c("#e10600", "grey")) +
  #   labs(x = NULL, y = NULL)
  
  df <- data.frame(Stop = time_intrvls)# %>% cumsum())
  df$Start <- c(0, df$Stop %>% magrittr::extract(., -length(.)))
  df$Stop[nrow(df)] <- 1/params$d0
  df$type <- c("OFF", "ON") %>% rep(., 10) %>%
    magrittr::extract(., i = 1:nrow(df)) %>%
    factor(., levels = c("ON", "OFF"))
  
  polymrase_arrvl %<>% subset(., t < time_intrvls %>% max())
  polymrase_arrvl$type <- apply(polymrase_arrvl, 1,
                                function(x) {
                                  subset(df, (Start < x[1]) & (x[1] < Stop))$type %>% 
                                    as.character()
                                })
  
  time_scale <- params %$%
    sum(k0, k1) %>%
    divide_by(1, .)
  
  
  count1 <- c(polymrase_arrvl$type %>%
                equals(., "ON") %>%
                sum(),
              count1)
  
  plot_layers <- list(coord_cartesian(xlim = c(0, time_intrvls %>% max())),
                      scale_x_continuous(breaks = time_intrvls %>%
                                           cumsum() %>%
                                           as.numeric() %>%
                                           c(0, .)),
                      theme(axis.line.x = element_line(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks.y = element_blank(),
                            legend.title = element_blank()
                            # legend.direction = "horizontal",
                            # legend.position = "bottom"
                      ),
                      geom_segment(data = df,
                                   aes(x = Start, xend = Stop, color = type),
                                   y = 0, yend = 0, size = 2),
                      # geom_segment(aes(x = t, y = y, xend =t, color = type), yend = 0),
                      # geom_point(),
                      ylab(NULL),
                      xlab(NULL)
                      # xlab("Time")
  )
  
  
  
  
  
  # df <- as.data.frame(lapply(df, unlist))
  # polymrase_arrvl <- as.data.frame(lapply(polymrase_arrvl, unlist))
  
  plot_layers2 <- list(coord_cartesian(xlim = c(0, df$Stop %>% tail(., 1)),
                                       ylim = c(-0.5, 1)),
                       # scale_x_continuous(breaks = time_intrvls %>%
                       #                      cumsum() %>%
                       #                      as.numeric() %>%
                       #                      c(0, .)),
                       theme(#axis.line.x = element_line(),
                         panel.grid = element_blank(),
                         axis.text = element_blank(),
                         axis.ticks.y = element_blank(),
                         legend.title = element_blank()
                         # legend.direction = "horizontal",
                         # legend.position = "bottom"
                       ),
                       geom_segment(data = df,
                                    aes(x = Start, xend = Stop, color = type),
                                    y = 0.05, yend = 0.05, size = 2),
                       geom_segment(aes(x = t, y = y, xend =t, color = type), yend = 0.15),
                       geom_point(),
                       ylab(NULL),
                       xlab("Time"))
  
  
  p_intrctn1 <- ggplot(polymrase_arrvl, aes(t, y))  +
    plot_layers2 +
    ggtitle(paste0("Births : ", count1[1])) +
    scale_color_manual(values = c("#e10600", "grey"), guide = FALSE) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.1))+
    annotate("segment", x=0, xend= df$Stop %>% tail(., 1),
             y=-0.08, yend = -0.08) +
    annotate("segment", x=0, xend= 0,
             y=-0.15, yend = -0.06) +
    annotate("segment", x= df$Stop %>% tail(., 1), xend= df$Stop %>% tail(., 1),
             y=-0.15, yend = -0.06) +
    annotate("text", x=0, y = -0.34, label = "0") +
    annotate("text", x= df$Stop %>% tail(., 1) %>% add(., -0.1), 
             y = -0.34, label="italic(T)[stat]", parse = TRUE)
  
  assign(paste0("p_intr_1", ix), p_intrctn1)
  
  nred <- params$k0/params$d0
  ngrey <- params$k1/params$d0
  w_urns <- 1/(params$k0+params$k1)
  set_of_recip_urns <- c(rep("OFF", ngrey), rep("ON", 6))
  df_urn <- data.frame(type = sample(set_of_recip_urns),
                       Start = c(0, w_urns*(1:(nred+ngrey-1))),
                       Stop = c(w_urns*(1:(nred+ngrey))))
  
  df_urn$type %<>%
    factor(., levels = c("ON", "OFF"))
  # previous <- df_urn$type[1]
  # shade <- c(paste0(previous, 1))
  # for (i in df_urn$type[-1]) {
  #   shade <- ifelse(i == previous, )
  # }
  
  polymrase_arrvl <-  rpois(1, params$v0/params$d0) %>% #params %$%
    #rpois(1, v0*(time_intrvls %>% sum())) %>%
    seq(to = .) %>%
    map_dbl(., ~ rexp(1, params$v0)) %>%
    cumsum() %>%
    data.frame(t = ., y = 0.9)
  
  polymrase_arrvl %<>% subset(., t < df_urn$Stop %>% max())
  polymrase_arrvl$type <- apply(polymrase_arrvl, 1,
                                function(x) {
                                  subset(df_urn, 
                                         (Start < x[1]) & (x[1] < Stop))$type %>% as.character()
                                })
  
  
  plot_layers3 <- list(coord_cartesian(xlim = c(0, df_urn$Stop %>% tail(., 1)),
                                       ylim = c(-0.5, 1)),
                       # scale_x_continuous(breaks = time_intrvls %>%
                       #                      cumsum() %>%
                       #                      as.numeric() %>%
                       #                      c(0, .)),
                       theme(#axis.line.x = element_line(),
                         panel.grid = element_blank(),
                         axis.text = element_blank(),
                         axis.ticks.y = element_blank(),
                         legend.title = element_blank()
                         # legend.direction = "horizontal",
                         # legend.position = "bottom"
                       ),
                       geom_segment(data = df_urn,
                                    aes(x = Start, xend = Stop, color = type),
                                    y = 0.05, yend = 0.05, size = 2),
                       geom_segment(aes(x = t, y = y, xend =t, color = type), yend = 0.15),
                       geom_point(),
                       ylab(NULL),
                       xlab("Time"))
  df_div <- data.frame(x=c(df_urn$Start, df$Stop %>% tail(., 1)),
                       xend=c(df_urn$Start, df$Stop %>% tail(., 1)))
  
  count1 <- polymrase_arrvl$type %>%
    equals(., "ON") %>%
    sum()
  
  p_intrctn2 <- ggplot(polymrase_arrvl, aes(t, y))  +
    plot_layers3 +
    ggtitle(paste0("Births : ", count1[1])) +
    scale_color_manual(values = c("#e10600", "grey"), guide = FALSE) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.1))+
    geom_segment(df_div, mapping = aes(x = x, xend = xend), y = 0.01, yend = 0.09) +
    annotate("segment", x=0, xend= df_urn$Stop %>% tail(., 1),
             y=-0.08, yend = -0.08) +
    annotate("segment", x=0, xend= 0,
             y=-0.15, yend = -0.06) +
    annotate("segment", x= df_urn$Stop %>% tail(., 1), xend= df_urn$Stop %>% tail(., 1),
             y=-0.15, yend = -0.06) +
    annotate("text", x=0, y = -0.34, label = "0") +
    annotate("text", x= df_urn$Stop %>% tail(., 1) %>% add(., -0.1), 
             y = -0.34, label="italic(T)[stat]", parse = TRUE)
  
  assign(paste0("p_intr_2", ix), p_intrctn2)
  
}

p4 <- ggplot()  +
  # geom_segment(df_urn, mapping = aes(x = Start, xend = Stop, color = type),
  #              y = 0, yend = 0, size = 2) +
  geom_rect(df_urn, mapping = aes(xmin = Start, xmax = Stop, fill = type), 
            ymin=0, ymax=0.01) +
  theme_nothing() +
  coord_cartesian(ylim = c(0, 0.01)) + 
  # scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  # scale_color_manual(values = c("#e10600", "grey")) +
  scale_fill_manual(values = c("#e10600", "grey")) +
  labs(x = NULL, y = NULL) +
  #annotate("segment", x=df_div$x, xend = df_div$xend, y = 0.04, yend = 0.07)
  geom_segment(df_div, mapping = aes(x = x, xend = xend), y = 0.003, yend = 0.007)  


p_connector <- ggplot(data.frame(), aes()) +
  theme_nothing() +
  annotate(geom = "segment",
           x = 0,
           xend = 0,
           y = 1,
           yend = 0,
           arrow = arrow(length = unit(0.1, "inches")))


p11 <- plot_grid(p1, p_connector, p3, 
                 nrow = 3,
                 rel_heights = c(0.8, 0.1, 0.1))



p12 <- plot_grid(p2, p_connector, p4, 
                 nrow = 3,
                 rel_heights = c(0.8, 0.1, 0.1))
p <- plot_grid(p11, p12, labels = c("b", "c"))

p_nothing <- ggplot() + theme_nothing()
pp0 <- plot_grid(p, p_nothing, labels = c(NA, NA), nrow = 2, rel_heights = c(1, 0.1))

pp1 <- plot_grid(p_intr_11, p_intr_12, p_intr_13, nrow = 3, labels = c("a", NA, NA))
pp2 <- plot_grid(p_intr_21, p_intr_22, p_intr_23, nrow = 3, labels = c("d", NA, NA))

p21 <- plot_grid(pp1, pp0, pp2, nrow = 1,
                 scale = 0.9, rel_heights = c(1, 1, 1))

# p21 <- plot_grid(pp1, )

# p21 <- plot_grid(p, p_intrctn1, p_intrctn2, nrow = 3, rel_heights = c(0.6, 0.3, 0.3),
#                  labels = c(NA, "c"))

ggsave("The_Urn_Scheme2.pdf",
       p21,
       width = 10.2,
       height = 4)
# ggsave("test.pdf", p11,
#        width = 1.7,
#        height = 3)
# ggsave("test2.pdf", p12,
#        width = 1.7,
#        height = 3)

#-----------------------------
#Nudge each of the T_i's by uniformly distributed random variables.
#The uniform distribution is U(-e, e).
#Take a range of values for e, e.g., 0, 0.001*(1/k0+k1), 0.005*(1/k0+k1), 0.01*(1/k0+k1),
#{0.05, 0.1, 0.5, 1, 1.5, 2, 5, 10, 15, 20}*(1/k0+k1)
#Call these nudge_widths.
#For each nudge_width, do 100000 samples.
# 
# set.seed(123234)
# dat_ref <- list()
# dat_nudged <- list()
# for (nudge_ratio in c(0.01, 0.1, 1, 10, 100)) {
#   params <- data.frame(k0 = 0.05,
#                        k1 = 0.2,
#                        d0 = 0.05,
#                        v0 = 0.05)
#   
#   count1 <- c()
#   count2 <- c()
#   
#   n_runs <- 1000
#   
#   for (i in 1:1000) {
#     if (!(i %% 10)) print(i)
#     time_intrvls <- params %$%
#       c(k0 = k0, k1 = k1) %>%
#       rep(., n_runs) %>%
#       map_dbl(~ rexp(1, .x))
#     
#     polymrase_arrvl <- params %$%
#       rpois(1, v0*(time_intrvls %>% sum())) %>%
#       seq(to = .) %>%
#       map_dbl(., ~ rexp(1, params$v0)) %>%
#       cumsum() %>%
#       data.frame(t = ., y = 1)
#     
#     df <- data.frame(Stop = time_intrvls %>% cumsum())
#     df$Start <- c(0, df$Stop %>% extract(., -length(.)))
#     df$type <- c("OFF", "ON") %>% rep(., n_runs) %>%
#       factor(., levels = c("ON", "OFF"))
#     
#     polymrase_arrvl %<>% subset(., t < time_intrvls %>% sum())
#     polymrase_arrvl$type <- apply(polymrase_arrvl, 1,
#                                   function(x) {
#                                     subset(df, (Start < x[1]) & (x[1] < Stop))$type %>% as.character()
#                                   })
#     
#     time_scale <- params %$%
#       sum(k0, k1) %>%
#       divide_by(1, .)
#     
#     
#     count1 <- c(polymrase_arrvl$type %>%
#                   equals(., "ON") %>%
#                   sum(),
#                 count1)
#     
#     nudge_by <- runif(nrow(polymrase_arrvl), min = -nudge_ratio*time_scale, max = nudge_ratio*time_scale)
#     polymrase_arrvl$t %<>% add(., nudge_by)
#     polymrase_arrvl$type <- apply(polymrase_arrvl, 1,
#                                   function(x) subset(df, (Start < as.numeric(x[1])) &
#                                                        (as.numeric(x[1]) < Stop))$type)
#     
#     count2 <- c(polymrase_arrvl$type %>%
#                   equals(., "ON") %>%
#                   sum(),
#                 count2)
# 
#   }
#   
#   dat_ref[[paste0("X_", nudge_ratio)]] <- count1
#   dat_nudged[[paste0("X_", nudge_ratio)]] <- count2
# }
# 
# dists <- mapply(kl.dist, dat_ref, dat_nudged)
