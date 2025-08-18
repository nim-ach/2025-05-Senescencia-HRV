## Load libraries
library(data.table)
library(brms)
library(ggplot2)

## Load data
data("senescence_d")

tmp_d <- senescence_d[, list(cd21_p_cd11c_m,
                     cd21_p_cd11c_p,
                     cd21_m_cd11c_m, 
                     cd21_m_cd11c_p,
                     age, sex, id)]

tmp_d <- melt.data.table(tmp_d, id.vars = c("age", "sex", "id"))

tmp_d[, age_cat := cut_number(age, n = 4)]

tmp_d[, age_cat := `levels<-`(age_cat, c("Q1\n(Younger)", "Q2", "Q3", "Q4\n(Older)"))]

tmp_d[, variable := `levels<-`(variable, c("CD21+ CD11c-", "CD21+ CD11c+",
                                           "CD21- CD11c-", "CD21- CD11c+"))]

p1 <- ggplot(tmp_d, aes(x = age_cat, y = value)) +
  facet_wrap(~variable, ncol = 2, scales = "free_y") + 
  geom_violin(aes(fill = age_cat), na.rm = TRUE, trim = F, width = 0.5, alpha = .75, adjust = 2) +
  geom_point(col = "white", position = position_jitter(width = 0.01), alpha = 0.5, na.rm = TRUE,
             size = 1) +
  geom_boxplot(aes(fill = age_cat), width = 0.1, na.rm = TRUE, outlier.colour = NA, alpha = .75) +
  scale_color_brewer(type = "qual", palette = 2, aesthetics = c("color", "fill")) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
  labs(x = "Age Quartiles", y = "Absolute cell counts") +
  see::theme_modern() +
  theme(legend.position = "none")

p2 <- ggplot(tmp_d, aes(x = sex, y = value)) +
  facet_wrap(~variable, ncol = 2, scales = "free_y") + 
  geom_violin(aes(fill = sex), na.rm = TRUE, trim = F, width = 0.5, alpha = .75, adjust = 2) +
  geom_point(col = "white", position = position_jitter(width = 0.01), alpha = 0.5, na.rm = TRUE,
             size = 1) +
  geom_boxplot(aes(fill = sex), width = 0.1, na.rm = TRUE, outlier.colour = NA, alpha = .75) +
  scale_color_brewer(type = "qual", palette = 2, aesthetics = c("color", "fill")) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
  labs(x = "Sex", y = "Absolute cell counts") +
  see::theme_modern() +
  theme(legend.position = "none")

fig_b_cell_dist <- 
  ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2,
                    heights = c(4,3),
                    labels = c("A.", "B."))

ggsave(plot = fig_b_cell_dist, 
       filename = "figures/fig_b_cell_dist.png",
       width = 7, height = 10, 
       device = "png", dpi = 400)
ggsave(plot = fig_b_cell_dist, 
       filename = "figures/fig_b_cell_dist.pdf",
       width = 7, height = 10)
