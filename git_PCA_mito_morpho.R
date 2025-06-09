


setwd("your worknig directory")

mt_morpho <- read_excel("your exel file.xlsx")
mt_morpho <- mt_morpho[ , c(2,3,6,7,13,14,17:188)]

mt_morpho <- mt_morpho %>% 
  mutate(cell_mark = paste(Row, Column, Field, `Object ID`, sep = "_")) %>% 
  mutate(cell_line = case_when(Column %in% c(1:2) ~ "A375",
                               Column %in% c(7:8) ~ "WM1341D")) %>% 
  mutate(condition = case_when(Column == 1 | Column == 7 ~ "zero",
                               Column == 2 | Column == 8 ~ "one")) %>% 
  mutate(condition = factor(condition, levels = c("zero", "one"))) %>% 
  arrange(condition) %>% 
  mutate(cell_line = factor(cell_line, levels = c("A375", "WM1341D"))) %>% 
  arrange(cell_line) %>% 
  mutate(AM_descrypt = paste(Row, Column))

rownames(mt_morpho) <- mt_morpho$cell_mark 

mt_morpho %>% select(condition) %>% table()

pca <- mt_morpho[, c(5:178)]

pca <- as.data.frame(sapply(pca, as.numeric))

#pca[is.na(pca)] <- 0

rownames(pca) <- mt_morpho$cell_mark
#------------------------
setwd("your working dir")
source("corr_clustering_1.R")
setwd("your worknig directory")
# -----------------------

result <- cluster_correlated_features(pca)

result$clusters[[2]]
result$cor_matrix

result$cluster_assignments
result$representatives


df <- as.data.frame(result$representatives)

mt_morpho2 <- mt_morpho %>% select(all_of(as.character(df$`result$representatives`))) 

mt_morpho2 <- as.data.frame(apply(mt_morpho2 ,2 ,as.numeric))
#mt_morpho2 <-  as.data.frame(apply(mt_morpho2, 2, round, 3))
#mt_morpho2[is.na(mt_morpho2)] <- 0

rownames(mt_morpho2) <- mt_morpho$cell_mark 

psych::cortest.bartlett(mt_morpho2)

kmo_results <- KMO(mt_morpho2)
# Overall KMO
kmo_results$MSA  # single value

# Per-variable MSA
per_variable_kmo <- as.data.frame(kmo_results$MSAi)

# Optional: Add variable names as a column
per_variable_kmo$variable <- rownames(per_variable_kmo)

# Optional: Rename column for clarity
colnames(per_variable_kmo)[1] <- "KMO_value"

variable_after_kmo <- per_variable_kmo %>% filter(KMO_value >= 0.6)

pca <- mt_morpho2 %>% select(all_of(as.character(variable_after_kmo$variable)))

pca_obj <- PCA(pca)

data_plot <- data.frame(cbind(PC1 = pca_obj$ind$coord[, 1], PC2 = pca_obj$ind$coord[, 2]))

data_plot$Groups <- mt_morpho$cell_line 
data_plot$condition <- mt_morpho$condition
data_plot$AM_descrypt <- mt_morpho$AM_descrypt
data_plot$Row <- mt_morpho$Row
data_plot$Column <- mt_morpho$Column
data_plot$Field <- mt_morpho$Field 


centroids <- data_plot %>% 
  group_by(Groups, condition) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2)) %>% 
  mutate(Groups_condi = paste(Groups, condition))

p_cell_line <- fviz_pca_ind(pca_obj, label = "none", habillage = factor(mt_morpho$cell_line), addEllipses  = T) 
p_cell_line
p_cell_data <- p_cell_line$data

identical(data_plot$PC1, p_cell_data$x)

ggplot(data_plot, aes(x = PC1, y = PC2, colour = Groups, alpha = Groups)) +
  geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Adds horizontal axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Adds vertical axis
  geom_point(data = centroids, aes(x = PC1, y = PC2), size = 2, alpha = 1, color = "black") +
  theme_minimal() +
  scale_alpha_manual(values = c("A375" = 0.3, "WM1341D" = 0.3)) +
  facet_wrap(~ interaction(Groups, condition)) +
  xlim(-10, 7) + 
  ylim(-7, 10) +
  xlab("PC1 14.5 %") +
  ylab("PC2 13.24 %") +
  theme(legend.position = "none")

ggsave("all_cells.png", width = 9, height = 10)

ggplot(centroids, aes(x = PC1, y = PC2, color = Groups, shape = Groups, fill = Groups, size = condition, alpha = condition))  + 
  geom_point(stroke = 1) +
  scale_color_manual(values = c("A375" = "green4", "WM1341D" = "red")) +
  scale_shape_manual(values = c("A375" = 22 , "WM1341D" = 24)) +
  scale_fill_manual(values = c("A375" = "green4", "WM1341D" = "red")) +
  scale_size_manual(values = c("zero" = 5, "one" = 6)) +
  scale_alpha_manual(values = c("zero" = 1, "one" = .4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  +
  xlab("PC1 14.5 %") +
  ylab("PC2 13.24 %") +
  theme_minimal() 

ggsave("all_cells_mean.png", width = 4, height = 5)

centroids_col <- data_plot %>% 
  filter(Column == 1 | Column == 2) %>% 
  group_by(Row, Column, condition) %>% 
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

ggplot(data_plot %>% filter(Groups == "A375"), aes(x = PC1, y = PC2, colour = Groups, alpha = Groups)) +
  geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Adds horizontal axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Adds vertical axis
  geom_point(data = centroids_col, aes(x = PC1, y = PC2), size = 2, alpha = 1, color = "black") +
  theme_minimal() +
  scale_alpha_manual(values = c("A375" = 0.5)) +
  facet_wrap(~condition) +
  xlim(-10, 7) + 
  ylim(-7, 10) +
  xlab("PC1 14.5 %") +
  ylab("PC2 13.24 %")


centroids_col2 <- data_plot %>% 
  filter(Column == 7 | Column == 8) %>% 
  group_by(Row, Column, condition) %>% 
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

ggplot(data_plot %>% filter(Groups == "WM1341D"), aes(x = PC1, y = PC2, colour = Groups, alpha = Groups)) +
  geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Adds horizontal axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Adds vertical axis
  geom_point(data = centroids_col2, aes(x = PC1, y = PC2), size = 2, alpha = 1, color = "black") +
  theme_minimal() +
  scale_alpha_manual(values = c("A375" = 0.5)) +
  facet_wrap(~condition) +
  xlim(-10, 7) + 
  ylim(-7, 10) +
  xlab("PC1 14.5 %") +
  ylab("PC2 13.24 %")

columns_centroids <- rbind(centroids_col, centroids_col2) %>% 
  mutate(Groups = case_when(Column == 1 | Column == 2 ~ "A375",
                            Column == 7 | Column == 8 ~ "WM1341D")) %>% 
  mutate(Groups_condi = paste(Groups, condition))

ggplot(columns_centroids, aes(x = PC1, y = PC2, color = Groups_condi, shape = Groups, fill = Groups_condi, size = condition, alpha = condition))+ 
  geom_point(stroke = 1) +
  geom_point(data = centroids, aes(x = PC1, y = PC2), size = 8, color = "black") +
  scale_color_manual(values = c("A375 zero" = "green4", "WM1341D zero" = "red", "A375 one" = "green4", "WM1341D one" = "red")) +
  scale_shape_manual(values = c("A375" = 22 , "WM1341D" = 24)) +
  scale_fill_manual(values = c("A375 zero" = "green4", "WM1341D zero" = "red", "A375 one" = "green4", "WM1341D one" = "red")) +
  scale_size_manual(values = c("zero" = 4, "one" = 4)) +
  scale_alpha_manual(values = c("zero" = 1, "one" = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  +
  xlab("PC1 14.5 %") +
  ylab("PC2 13.24 %") +
  theme_minimal() +
  theme(legend.position = "none")


ggsave("all_cells_mean_wells.png", width = 4, height = 5)


# ------------------------------------------------ for grant application
ggplot(columns_centroids, aes(x = PC1, y = PC2, color = Groups_condi, shape = Groups, fill = Groups_condi, size = condition, alpha = condition))+ 
  geom_point(stroke = 1) +
  geom_point(data = centroids, aes(x = PC1, y = PC2), size = 6, alpha = 0.5) +
  scale_color_manual(values = c("A375 zero" = "black", "WM1341D zero" = "black", "A375 one" = "red", "WM1341D one" = "red")) +
  scale_shape_manual(values = c("A375" = 22 , "WM1341D" = 24)) +
  scale_fill_manual(values = c("A375 zero" = "black", "WM1341D zero" = "black", "A375 one" = "red", "WM1341D one" = "red")) +
  scale_size_manual(values = c("zero" = 2, "one" = 2)) +
  scale_alpha_manual(values = c("zero" = 1, "one" = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")  +
  xlab("PC1 14.5 %") +
  ylab("PC2 13.24 %") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, hjust = 1)) +
  facet_wrap(~ Groups)

ggsave("PCA.png", width = 3, height = 3)







