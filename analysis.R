# ==============================================================================
# 项目：ST001683 代谢组学二次挖掘 (Secondary Mining)
# ==============================================================================

# --- 0. 环境准备 ---
packages <- c("jsonlite", "dplyr", "tidyr", "ggplot2", "stringr", "gridExtra", "ggpubr", "tibble")

library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gridExtra) # 用于拼图
library(ggpubr)    # 用于统计标注
library(tibble)
library(impute)    # 用于KNN填补
library(ggVennDiagram)
library(pheatmap)
library(pROC)

# ==============================================================================
# 第一部分：数据抓取与全量找回 (Data Fetching & Recovery)
# ==============================================================================
cat(">>> [Step 1] 开始抓取原始数据 (含盲肠/粪便)...\n")

# 1.1 获取元数据
factors_url <- "https://www.metabolomicsworkbench.org/rest/study/study_id/ST001683/factors"
metadata_raw <- bind_rows(fromJSON(factors_url))

# 1.2 获取代谢物数据
data_url <- "https://www.metabolomicsworkbench.org/rest/study/study_id/ST001683/data"
parsed_json <- fromJSON(data_url)
data_values <- bind_rows(lapply(parsed_json, function(x) x$DATA))
# 处理重名代谢物
metab_names <- make.unique(as.character(sapply(parsed_json, function(x) x$metabolite_name)))

# 1.3 构建矩阵并强制转数值
abundance_raw <- cbind(Metabolite = metab_names, data_values)
abundance_raw[,-1] <- lapply(abundance_raw[,-1], function(x) as.numeric(as.character(x)))
abundance_t <- as.data.frame(t(abundance_raw[, -1]))
colnames(abundance_t) <- metab_names

# 1.4 定义强力 ID 清洗函数 
clean_id_robust <- function(x) {
  x <- as.character(x)
  x <- gsub("^X", "", x)       # 去掉 R 自动加的 X
  x <- gsub("\\.", "", x)      # 去掉点
  x <- gsub("-", "", x)        # 去掉横杠
  x <- gsub("_", "", x)        # 去掉下划线
  x <- toupper(x)              # 转大写
  return(x)
}

# 1.5 整理与合并 
abundance_ready <- abundance_t %>%
  rownames_to_column("Raw_ID") %>%
  mutate(Join_ID = clean_id_robust(Raw_ID)) %>%
  group_by(Join_ID) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

metadata_ready <- metadata_raw %>%
  mutate(Join_ID = clean_id_robust(local_sample_id)) %>%
  select(Sample_ID = local_sample_id, Factors = factors, Join_ID)

final_df_raw <- metadata_ready %>%
  inner_join(abundance_ready, by = "Join_ID")

# 1.6 解析样本类型
final_df <- final_df_raw %>%
  mutate(
    Organ = case_when(
      grepl("serum", Factors, ignore.case=T) ~ "Serum",
      grepl("urine", Factors, ignore.case=T) ~ "Urine",
      grepl("feces|fecal", Factors, ignore.case=T) ~ "Feces",
      grepl("caecal|cecum", Factors, ignore.case=T) ~ "Caecal",
      TRUE ~ "Other"
    ),
    Grp = case_when(
      grepl("Colonization:Bt", Factors) ~ "Bt",
      grepl("Colonization:germ-free", Factors) ~ "GF",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Grp %in% c("Bt", "GF"), Organ != "Other")

cat(paste(">>> 数据整合完成。总样本数:", nrow(final_df), "\n"))
print(table(final_df$Organ, final_df$Grp))

# ==============================================================================
# 第二部分：鲁棒性预处理管道 (Robust Preprocessing)
# ==============================================================================
cat("\n>>> [Step 2] 执行数据清洗与标准化...\n")

# 提取代谢物矩阵
meta_cols <- c("Sample_ID", "Factors", "Join_ID", "Organ", "Grp")
data_mat <- final_df %>% select(-any_of(meta_cols))

# 2.1 清除全 NA 列和极低方差列
data_mat[is.na(data_mat)] <- NA # 确保格式统一
valid_cols <- colSums(!is.na(data_mat)) > 5 # 至少5个样本有值
data_mat <- data_mat[, valid_cols]

# 2.2 混合填补 (Min/2 + KNN)
imputed_mat <- data_mat
# 简单策略：如果某列缺失严重，用最小值一半填补(LOD)
for(col in names(data_mat)) {
  x <- data_mat[[col]]
  if(mean(is.na(x)) > 0.5) {
    min_val <- min(x, na.rm = TRUE)
    if(is.infinite(min_val)) min_val <- 0
    imputed_mat[is.na(x), col] <- min_val / 2
  }
}
# 剩余随机缺失用 KNN
if(sum(is.na(imputed_mat)) > 0) {
  # 捕获 KNN 可能的报错 (如果数据太稀疏)
  tryCatch({
    knn_res <- impute.knn(as.matrix(t(imputed_mat)), k = 5)
    imputed_mat <- as.data.frame(t(knn_res$data))
  }, error = function(e) {
    # 如果 KNN 失败，回退到最小值填补
    imputed_mat[is.na(imputed_mat)] <- 0 
  })
}

# 2.3 Log2 转化 
log_mat <- log2(imputed_mat + 1)
# 合并回主数据框，用于后续画 Boxplot
final_df_processed <- cbind(final_df %>% select(all_of(meta_cols)), log_mat)

# ==============================================================================
# 第三部分：代谢重塑全景图 (Landscape PCA)
# ==============================================================================
# [第 3a 部分] PCA

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr) # 核心：用于合并图例

cat("\n>>> [Step 3a] 正在生成含图例的 PCA 拼图并导出...\n")

run_pca_single <- function(organ_name) {
  # 1. 数据准备
  sub_df <- final_df_processed %>% filter(Organ == organ_name)
  pca_input <- sub_df %>% select(-any_of(meta_cols))
  pca_input <- pca_input[, sapply(pca_input, sd) > 0] # 过滤零方差
  
  # 2. Pareto Scaling
  pca_scaled <- as.data.frame(lapply(pca_input, function(x) (x - mean(x)) / sqrt(sd(x))))
  
  # 3. PCA 计算
  pca <- prcomp(pca_scaled, center = FALSE, scale. = FALSE)
  var_pct <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  
  # 4. 绘图 (注意：这里我们暂时不隐藏 legend，由后续 ggarrange 统一处理)
  plot_dat <- cbind(as.data.frame(pca$x), sub_df)
  
  ggplot(plot_dat, aes(x = PC1, y = PC2, color = Grp, fill = Grp)) +
    geom_point(size = 3.5, shape = 21, color = "white", alpha = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, linewidth = 0.3) +
    scale_color_manual(values = c("Bt" = "#E41A1C", "GF" = "#377EB8"), name = "Colonization Status") +
    scale_fill_manual(values = c("Bt" = "#E41A1C", "GF" = "#377EB8"), name = "Colonization Status") +
    labs(title = organ_name,
         x = paste0("PC1 (", var_pct[1], "%)"), 
         y = paste0("PC2 (", var_pct[2], "%)")) +
    theme_bw() + 
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

# 生成四个子图
p1 <- run_pca_single("Caecal")
p2 <- run_pca_single("Feces")
p3 <- run_pca_single("Serum")
p4 <- run_pca_single("Urine")

# ==============================================================================
# 使用 ggarrange 进行智能拼图 (自动提取公共图例)
# ==============================================================================
p_final <- ggarrange(
  p1, p2, p3, p4, 
  ncol = 2, nrow = 2, 
  common.legend = TRUE,      # 开启公共图例
  legend = "bottom",         # 图例放在底部
  labels = c("A", "B", "C", "D") # 给每个子图加上 A,B,C,D 标签 (SCI 标准)
)

# 添加总标题
p_final <- annotate_figure(p_final,
               top = text_grob("Landscape of Metabolic Remodeling by Bt Colonization", 
                               face = "bold", size = 16))

# ==============================================================================
# 高保真导出
# ==============================================================================
# 导出 PDF
ggsave("Figure_1_PCA_Final.pdf", plot = p_final, width = 10, height = 9, device = "pdf")

# 导出 300 DPI PNG
ggsave("Figure_1_PCA_Final.png", plot = p_final, width = 10, height = 9, dpi = 300)

cat(">>> 导出完成！已生成 Figure_1_PCA_Final.pdf 和 .png (含图例与 ABCD 标注)\n")

# ==============================================================================
# [第 3b 部分] 跨器官显著性全量分析 (Venn Diagram + Core List)
# 目标：证明系统性影响，提取 41 个核心分子并核对明星路径
# ==============================================================================
library(ggVennDiagram)
library(dplyr)

cat("\n>>> [Step 3b] 正在执行 Venn 交叉分析并生成投稿级图表...\n")

# 1. 提取显著差异代谢物函数 (保留你原始的 T-test 逻辑)
get_sig_mets <- function(organ_name) {
  sub_df <- final_df_processed %>% filter(Organ == organ_name)
  mets <- setdiff(colnames(sub_df), meta_cols)
  
  p_vals <- sapply(mets, function(m) {
    tryCatch(t.test(sub_df[[m]] ~ sub_df$Grp)$p.value, error = function(e) 1)
  })
  
  # 筛选 P < 0.05 的代谢物
  return(mets[p_vals < 0.05])
}

# 2. 生成四个器官的显著差异列表
venn_list <- list(
  Serum  = get_sig_mets("Serum"),
  Urine  = get_sig_mets("Urine"),
  Caecal = get_sig_mets("Caecal"),
  Feces  = get_sig_mets("Feces")
)

# 3. 绘制投稿级 Venn 图 (优化颜色、边框与标注)
p_venn_final <- ggVennDiagram(venn_list, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "#F7FBFF", high = "#084594", name = "Count") +
  # 设置黑色边框增强印刷感
  scale_color_manual(values = rep("black", 4)) +
  labs(title = "Overlap of Significant Metabolites (P < 0.05)",
       subtitle = "Demonstrating Systemic Metabolic Consistency") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14))

# 4. 高保真导出 (PDF 矢量图 + 300 DPI PNG)
ggsave("Figure_2_Venn_Overlap.pdf", plot = p_venn_final, width = 8, height = 7, device = "pdf")
ggsave("Figure_2_Venn_Overlap.png", plot = p_venn_final, width = 8, height = 7, dpi = 300)

# 5. 核心代谢物提取逻辑 (保留并增强)
core_mets <- Reduce(intersect, venn_list)
cat("\n>>> 四个器官共同显著的核心代谢物数量:", length(core_mets), "\n")

# 6. [新增] 自动保存核心分子列表为 Excel/CSV (用于 Supplementary Table S3)
core_mets_df <- data.frame(Core_Metabolites = core_mets)
write.csv(core_mets_df, "Table_S3_Core_41_Metabolites.csv", row.names = FALSE)

# 7. 明星分子路径核对 (保留你原代码的精华)
stars <- c("TRYPTOPHAN", "5-HYDROXYINDOLE", "INDOLEPROPIONIC", "N-ACETYLPUTRESCINE")
found_stars <- stars[sapply(stars, function(s) any(grepl(s, core_mets)))]

# 结果总结打印
cat("\n======================================================\n")
cat("CORE METABOLITES SUMMARY:\n")
cat("------------------------------------------------------\n")
cat("跨器官核心集内包含的明星分子:\n")
print(found_stars)
cat("\n核心集列表预览 (Top 10):\n")
print(head(core_mets, 10))
cat("======================================================\n")

print(p_venn_final)

# ==============================================================================
# ==============================================================================
# [第 3c 部分] 创新点挖掘：代谢稳态收缩全量分析
# 目标：量化证明 Bt 定植作为“代谢定海神针”，显著降低了宿主代谢的随机波动
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

cat(">>> [Step 3c] 正在执行全量 CV 变异收缩分析并生成图表...\n")

# 1. 计算变异系数 (CV = SD / Mean) - 保持原始逻辑
cv_analysis <- final_df_processed %>%
  pivot_longer(cols = -all_of(meta_cols), names_to = "Metabolite", values_to = "Value") %>%
  group_by(Organ, Grp, Metabolite) %>%
  summarise(
    sd_val = sd(Value, na.rm = TRUE),
    mean_val = mean(Value, na.rm = TRUE),
    cv = sd_val / mean_val,
    .groups = "drop"
  ) %>%
  filter(!is.na(cv), !is.infinite(cv))

# 2. 统计稳定化趋势 (保留并增强统计维度)
stability_summary <- cv_analysis %>%
  pivot_wider(id_cols = c(Organ, Metabolite), names_from = Grp, values_from = cv) %>%
  filter(!is.na(Bt), !is.na(GF)) %>%
  mutate(is_stabilized = Bt < GF) %>%
  group_by(Organ) %>%
  summarise(
    Total_Mets = n(),
    Stabilized_Percent = round(mean(is_stabilized) * 100, 1),
    # [新增] 计算中位数变异收缩率 (Reduction Intensity)
    Median_GF_CV = round(median(GF), 4),
    Median_Bt_CV = round(median(Bt), 4),
    Median_Reduction_Pct = round(median((GF - Bt) / GF, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

# 打印并保存统计结果 (用于 Results 3.1b 的写作支撑)
cat("\n>>> 代谢稳态收敛统计表 (Bt vs GF):\n")
print(stability_summary)
write.csv(stability_summary, "Table_S2_CV_Stability_Summary.csv", row.names = FALSE)

# 3. 绘制出版级对比图 (修复 linewidth 警告，优化布局)
p_cv_final <- ggplot(cv_analysis, aes(x = Grp, y = cv, fill = Grp)) +
  # 使用 linewidth 适配新版 ggplot2
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5, linewidth = 0.6) +
  # 聚焦 95% 置信区间内的数据，避免极端离散点破坏构图
  coord_cartesian(ylim = c(0, quantile(cv_analysis$cv, 0.95, na.rm = TRUE))) + 
  facet_wrap(~Organ, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Bt" = "#E41A1C", "GF" = "#377EB8")) +
  # 添加统计检验
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     label.y.npc = "top", size = 5) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Systemic Metabolic Variance Compression",
       subtitle = "Bt colonization suppresses host individual stochasticity",
       y = "Coefficient of Variation (CV)")

# 4. 高保真导出
ggsave("Figure_3_CV_Compression.pdf", plot = p_cv_final, width = 10, height = 5, device = "pdf")
ggsave("Figure_3_CV_Compression.png", plot = p_cv_final, width = 10, height = 5, dpi = 300)

cat(">>> 导出完成：Figure_3_CV_Compression.pdf 和 Table_S2_CV_Stability_Summary.csv\n")

# ==============================================================================
# [第 3d 部分] Tryptophan 路径稳定性深度分析
# ==============================================================================
library(ggplot2); library(dplyr); library(tidyr); library(ggrepel)

cat(">>> [Step 3d] 正在生成 Tryptophan 路径稳定性分析 (全量无损版)...\n")

# 1. 提取并清洗数据
trp_keywords <- "TRYPTOPHAN|INDOLE|KYNURENINE|SEROTONIN|5-HYDROXY"

trp_stability_full <- cv_analysis %>%
  filter(Organ == "Urine", grepl(trp_keywords, Metabolite, ignore.case = TRUE)) %>%
  pivot_wider(id_cols = Metabolite, names_from = Grp, values_from = cv) %>%
  # 关键：移除无法计算 CV 的代谢物 (Mean 为 0 的情况)
  filter(!is.na(Bt), !is.na(GF)) %>%
  mutate(
    # 计算缩减率，并处理分母为 0 的异常情况 (防止产生 Inf)
    Reduction_Pct = ifelse(GF == 0, 0, (GF - Bt) / GF * 100),
    Stabilized = ifelse(Bt < GF, "Yes (Stabilized)", "No (Increased Noise)")
  ) %>%
  # 彻底移除任何计算残留的 NA 或非有限数值
  filter(is.finite(Reduction_Pct))

# 2. 统计总结
trp_summary_final <- trp_stability_full %>%
  summarise(
    Total_Mets = n(),
    Highly_Stabilized = sum(Bt < 0.05), 
    Significant_Reduction = sum(Reduction_Pct > 50),
    Median_Reduction = round(median(Reduction_Pct, na.rm = TRUE), 1)
  )

# 3. 保存表格 (对齐命名规范)
write.csv(trp_stability_full, "Table_S3_Trp_Stability_Details.csv", row.names = FALSE)
cat("\n>>> 统计摘要 (完全无损):\n"); print(trp_summary_final)

# 4. 绘图 (采用完全自适应坐标系)
p_trp_dot <- ggplot(trp_stability_full, aes(x = GF, y = Bt)) +
  # 绘制 1:1 对角线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  # 散点：大小代表收缩比例
  geom_point(aes(color = Stabilized, size = abs(Reduction_Pct)), alpha = 0.7) +
  # 智能标注
  geom_text_repel(aes(label = ifelse(Reduction_Pct > 80, Metabolite, "")), 
                  size = 3.5, 
                  box.padding = 0.5, 
                  max.overlaps = Inf, # 允许所有标注排布
                  min.segment.length = 0) +
  scale_color_manual(values = c("Yes (Stabilized)" = "#E41A1C", "No (Increased Noise)" = "gray60")) +
  # 设置气泡大小范围，增强视觉对比
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw(base_size = 14) +
  # 使用 coord_equal 但不设范围，让 R 自动包含所有数据点
  coord_equal() + 
  labs(title = "Variance Compression: Tryptophan Pathway",
       subtitle = paste0("Pathway-wide stability: ", trp_summary_final$Median_Reduction, "% Median Reduction"),
       x = "GF Group Variation (CV)", y = "Bt Group Variation (CV)",
       size = "Reduction %") +
  theme(legend.position = "right", 
        panel.grid.minor = element_blank())

# 5. 高保真并行导出
ggsave("Figure_4_Trp_Stability.pdf", plot = p_trp_dot, width = 8, height = 7, device = "pdf")
ggsave("Figure_4_Trp_Stability.tiff", plot = p_trp_dot, width = 8, height = 7, dpi = 300, compression = "lzw")
ggsave("Figure_4_Trp_Stability.png", plot = p_trp_dot, width = 8, height = 7, dpi = 300)

cat(">>> 导出完成：Figure_4 相关文件。警告应已消除。\n")

# ==============================================================================
# 第四部分：跨器官空间追踪分析 (Spatial Translocation Tracking)
# 目标：追踪核心产物 (IPA) 与前体 (5-HTP) 在四种基质中的动态流转
# ==============================================================================
library(ggplot2); library(dplyr); library(tidyr); library(ggpubr)

cat(">>> [Step 3.2b] 正在执行跨器官空间追踪分析...\n")

# 1. 优化版匹配函数：精准锁定内源性代谢物 (排除 D5 等同位素)
get_target_col <- function(keyword) {
  cols <- colnames(final_df_processed)
  matches <- grep(keyword, cols, value = TRUE, ignore.case = TRUE)
  matches <- matches[!grepl("D[0-9]", matches)] # 排除同位素
  return(matches[1])
}

# 2. 核心追踪函数：支持统计输出与高质量绘图
plot_spatial_optimized <- function(keyword, title_label) {
  col_name <- get_target_col(keyword)
  if(is.na(col_name)) { cat("Warning: Not found", keyword, "\n"); return(NULL) }
  
  # 提取并排序器官
  track_dat <- final_df_processed %>%
    select(Organ, Grp, Value = all_of(col_name)) %>%
    mutate(Organ = factor(Organ, levels = c("Caecal", "Feces", "Serum", "Urine")))
  
  # 绘制高质量箱线图
  p <- ggplot(track_dat, aes(x = Organ, y = Value, fill = Grp)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.6) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 1.2, alpha = 0.4) +
    stat_compare_means(aes(group = Grp), label = "p.signif", method = "t.test", hide.ns = FALSE) +
    scale_fill_manual(values = c("Bt" = "#E41A1C", "GF" = "#377EB8")) +
    labs(title = paste(" ", title_label),
         subtitle = paste("Metabolite ID:", col_name),
         y = "Log2 Abundance", x = "") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  
  return(list(plot = p, data = track_dat, id = col_name))
}

# 3. 执行分析并获取数据
res_ipa  <- plot_spatial_optimized("INDOLEPROPIONIC", "Spatial Tracking: Indolepropionic Acid (IPA)")
res_5htp <- plot_spatial_optimized("5-HYDROXY-TRYPTOPHAN", "Spatial Tracking: Host 5-HTP Depletion")

# 4. 生成跨器官统计汇总表 (Table S5)
cat(">>> 正在生成跨器官丰度统计表...\n")
summary_stats <- bind_rows(
  res_ipa$data %>% mutate(Metabolite = res_ipa$id),
  res_5htp$data %>% mutate(Metabolite = res_5htp$id)
) %>%
  group_by(Metabolite, Organ, Grp) %>%
  summarise(Mean = mean(Value), SD = sd(Value), .groups = "drop")

write.csv(summary_stats, "Table_S5_Spatial_Tracking_Summary.csv", row.names = FALSE)

# 5. 高清导出 (PDF + TIFF 300 DPI)
ggsave("Figure_6_Spatial_IPA.pdf", res_ipa$plot, width = 8, height = 5)
ggsave("Figure_6_Spatial_IPA.tiff", res_ipa$plot, width = 8, height = 5, dpi = 300, compression = "lzw")
ggsave("Figure_6_Spatial_5HTP.pdf", res_5htp$plot, width = 8, height = 5)
ggsave("Figure_6_Spatial_5HTP.tiff", res_5htp$plot, width = 8, height = 5, dpi = 300, compression = "lzw")

print(res_ipa$plot)
print(res_5htp$plot)

# ==============================================================================
# 第五部分：尿液核心功能签名分析 (Urinary Functional Signatures)
# 目标：定量展示 Bt 定植的标志性代谢产出与细菌生物量信号
# ==============================================================================
library(ggplot2); library(dplyr); library(tidyr); library(ggpubr)

cat(">>> [Step 3.2c] 正在执行尿液明星分子特写分析...\n")

# 1. 自动定位尿液中的核心差异物 (增加精确匹配逻辑)
# fMet: 细菌负荷信号; 5-HI: 核心色氨酸产物; N-Ac-Put: 脱毒产物
key_targets <- c("N-FORMYL-METHIONINE", "5-HYDROXYINDOLE", "N-ACETYLPUTRESCINE")

# 辅助函数：抓取最干净的列名
get_clean_cols <- function(keywords) {
  all_cols <- colnames(final_df_processed)
  matched <- unlist(lapply(keywords, function(k) {
    # 找匹配项，排除同位素，优先取不带后缀的
    m <- grep(k, all_cols, value = TRUE, ignore.case = TRUE)
    m <- m[!grepl("D[0-9]", m)]
    return(m[1])
  }))
  return(matched[!is.na(matched)])
}

target_cols <- get_clean_cols(key_targets)

# 2. 提取并整理尿液数据
urine_dat <- final_df_processed %>% 
  filter(Organ == "Urine") %>%
  select(Grp, all_of(target_cols)) %>%
  pivot_longer(-Grp, names_to = "Metabolite", values_to = "Abundance") %>%
  # 规范化名称：将下划线转为空格，去掉后缀以美化图表
  mutate(Metabolite = gsub("\\.[0-9]+", "", Metabolite),
         Metabolite = gsub("_", " ", Metabolite))

# 3. 绘制出版级分面箱线图
p_signatures_final <- ggplot(urine_dat, aes(x = Grp, y = Abundance, fill = Grp)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  facet_wrap(~Metabolite, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "t.test", label = "p.signif", label.y.npc = "top") +
  scale_fill_manual(values = c("Bt" = "#E41A1C", "GF" = "#377EB8")) +
  labs(title = "Key Functional Signatures in Urinary Profiles",
       subtitle = "Non-invasive readout of microbial biomass and secondary metabolism",
       y = "Log2 Abundance", x = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())

# 4. 生成统计统计汇总表 (Table S6)
cat(">>> 正在生成尿液特征统计表...\n")
urine_stats <- urine_dat %>%
  group_by(Metabolite, Grp) %>%
  summarise(Mean = mean(Abundance), SD = sd(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Grp, values_from = c(Mean, SD)) %>%
  mutate(Log2FC = Mean_Bt - Mean_GF)

write.csv(urine_stats, "Table_S6_Urinary_Signatures_Stats.csv", row.names = FALSE)

# 5. 高清导出
ggsave("Figure_7_Urine_Signatures.pdf", p_signatures_final, width = 10, height = 5)
ggsave("Figure_7_Urine_Signatures.tiff", p_signatures_final, width = 10, height = 5, dpi = 300, compression = "lzw")

print(p_signatures_final)

# ==============================================================================
# 第六部分：代谢网络分析
# ==============================================================================
# ==============================================================================
# [第 3.4 部分] 最终全量版：协同网络 (6a) + 量效耦合 (6b)
# 修复：完美合并两部分逻辑，解决热图报错，输出 Table S8 + 高清 PDF/TIFF
# ==============================================================================
library(pheatmap); library(ggplot2); library(dplyr); library(tidyr); library(ggpubr)

cat(">>> [Step 3.4] 正在执行全量代谢协同与量效耦合分析...\n")

# ------------------------------------------------------------------------------
# PART 1: [6a] 代谢网络分析 (Heatmap)
# ------------------------------------------------------------------------------

# 1. 明星分子名单 (完全保留你的“全明星阵容”)
target_map_list <- list(
  "TRYPTOPHAN" = "Tryptophan", "KYNURENINE" = "Kynurenine", "SEROTONIN" = "Serotonin",
  "INDOLEPROPIONIC" = "IPA", "5-HYDROXYINDOLE" = "5-HI",
  "PUTRESCINE" = "Putrescine", "N-ACETYLPUTRESCINE" = "N-Ac-Put",
  "N-FORMYL-METHIONINE" = "fMet"
)

# 2. 稳健提取函数 (无损保留你的原代码逻辑)
get_curated_network_data <- function(df, group_name) {
  sub_df <- df %>% filter(Grp == group_name, Organ %in% c("Serum", "Urine"))
  extracted_mat <- data.frame(row.names = 1:nrow(sub_df))
  for (key in names(target_map_list)) {
    candidates <- grep(key, colnames(df), value = TRUE)
    candidates <- candidates[!grepl("D[0-9]", candidates)]
    if (length(candidates) > 0) extracted_mat[[target_map_list[[key]]]] <- sub_df[[candidates[1]]]
  }
  return(extracted_mat)
}

# 3. 绘图与数据导出
plot_curated_heatmap <- function(group_name, filename) {
  mat <- get_curated_network_data(final_df_processed, group_name)
  # 【核心修正】：移除方差为 0 的列（解决 hclust 报错）
  mat <- mat[, sapply(mat, sd, na.rm=TRUE) > 0, drop=FALSE]
  
  if(ncol(mat) < 2) return(NULL) # 如果全是常量则跳过
  
  # 计算 Spearman 相关性
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
  
  # 绘制并保存 PDF (300 DPI 级别精度)
  pheatmap(cor_mat, 
           main = paste0("Metabolic Network: ", group_name),
           color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = TRUE, number_format = "%.2f",
           cluster_rows = TRUE, cluster_cols = TRUE,
           cellwidth = 35, cellheight = 35, border_color = "black",
           filename = filename, width = 8, height = 7)
  return(cor_mat)
}

cor_bt <- plot_curated_heatmap("Bt", "Figure_9_Network_Bt.pdf")
cor_gf <- plot_curated_heatmap("GF", "Figure_9_Network_GF.pdf")

# 保存相关性表格 (Table S8)
if(!is.null(cor_bt)) write.csv(cor_bt, "Table_S8_Bt_Correlation_Matrix.csv")

# ------------------------------------------------------------------------------
# PART 2: [6b] 细菌负荷与功能的量效相关性 (Scatter Plot)
# ------------------------------------------------------------------------------
cat(">>> 正在生成 fMet vs 5-HI 量效耦合散点图...\n")

# 1. 精准定位列名并提取 Bt 组数据
col_fmet_raw <- grep("N-FORMYL-METHIONINE", colnames(final_df_processed), value = TRUE)[1]
col_5hi_raw  <- grep("5-HYDROXYINDOLE", colnames(final_df_processed), value = TRUE)[1]

cor_dat_scatter <- final_df_processed %>% 
  filter(Grp == "Bt", Organ %in% c("Urine", "Serum")) %>%
  select(Organ, 
         fMet_Signal = all_of(col_fmet_raw), 
         HI5_Product = all_of(col_5hi_raw))

# 2. 绘制高质量 ggscatter (无损保留你的统计参数)
p_fmet_final <- ggscatter(cor_dat_scatter, x = "fMet_Signal", y = "HI5_Product",
          add = "reg.line",                         
          conf.int = TRUE,                          
          color = "Organ", palette = "jco",         
          shape = "Organ",
          cor.coef = TRUE, cor.method = "spearman", 
          cor.coeff.args = list(label.sep = "\n", size = 5)) +
  labs(title = "Bacterial Workload vs. Systemic Output",
       subtitle = "Quantitative coupling between fMet (biomass) and 5-HI (function)",
       x = "Bacterial Initiation Signal (fMet) Log2 Intensity",
       y = "Microbial Functional Product (5-HI) Log2 Intensity") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

# 3. 高保真导出 (PDF + TIFF 300 DPI)
ggsave("Figure_10_fMet_Coupling.pdf", plot = p_fmet_final, width = 8, height = 7)
ggsave("Figure_10_fMet_Coupling.tiff", plot = p_fmet_final, width = 8, height = 7, dpi = 300, compression = "lzw")

print(p_fmet_final)
cat(">>> [Step 3.4] 全部图表与 Table S8 已生成并保存。\n")

# ==============================================================================
# 第七部分：高通量标志物扫描与置换验证 (High-Throughput AUC & Permutation)
# ==============================================================================
library(pROC); library(dplyr); library(tidyr); library(ggplot2)

cat(">>> [Step 3.2a] 正在执行全量代谢物 AUC 扫描与显著性验证...\n")

# 1. 执行全量 AUC 扫描函数 (增加变化倍数计算)
run_auc_scan <- function(df, target_organ) {
  sub_df <- df %>% filter(Organ == target_organ)
  mets <- setdiff(colnames(sub_df), meta_cols)
  
  lapply(mets, function(m) {
    # 计算 AUC
    r <- tryCatch({ roc(sub_df$Grp, sub_df[[m]], quiet=T, direction="auto") }, error=function(e) NULL)
    if(is.null(r)) return(NULL)
    
    # 计算倍数变化 (Log2FC)
    m_bt <- mean(sub_df[[m]][sub_df$Grp == "Bt"], na.rm=T)
    m_gf <- mean(sub_df[[m]][sub_df$Grp == "GF"], na.rm=T)
    
    data.frame(Metabolite = m, Organ = target_organ, AUC = as.numeric(r$auc), Log2FC = m_bt - m_gf)
  }) %>% bind_rows()
}

all_aucs_raw <- bind_rows(run_auc_scan(final_df_processed, "Urine"), 
                          run_auc_scan(final_df_processed, "Serum"))

# 2. 【核心优化】过滤人工同位素与噪声，保留生物学代谢物
# 排除包含 D5, D3, D9, INTERNAL, IS 等字符的非内源性分子
all_aucs_clean <- all_aucs_raw %>%
  filter(!grepl("D[0-9]|INTERNAL|IS_|C13|STANDARD", Metabolite, ignore.case = TRUE)) %>%
  arrange(desc(AUC))

# 3. 【核心新增】置换检验函数 (Permutation Test)
# 针对 AUC > 0.95 的 Top 分子执行 1000 次随机打乱验证
run_permutation <- function(df, met_name, organ_name, n_perm = 1000) {
  sub_df <- df %>% filter(Organ == organ_name)
  actual_auc <- as.numeric(roc(sub_df$Grp, sub_df[[met_name]], quiet=T)$auc)
  null_aucs <- replicate(n_perm, as.numeric(roc(sample(sub_df$Grp), sub_df[[met_name]], quiet=T)$auc))
  p_val <- sum(null_aucs >= actual_auc) / n_perm
  return(p_val)
}

# 4. 提取 Top 20 标志物并进行置换验证
cat(">>> 正在对 Top 标志物执行置换检验 (可能需要 1-2 分钟)...\n")
top_biomarkers <- all_aucs_clean %>% head(20)
top_biomarkers$Permutation_P <- sapply(1:nrow(top_biomarkers), function(i) {
  run_permutation(final_df_processed, top_biomarkers$Metabolite[i], top_biomarkers$Organ[i])
})

# 5. 保存 Supplementary Table S4 (全量 + 置换 P 值)
# 将 P 值格式化为学术风格
top_biomarkers <- top_biomarkers %>%
  mutate(Permutation_P = ifelse(Permutation_P == 0, "< 0.001", sprintf("%.3f", Permutation_P)))

write.csv(top_biomarkers, "Table_S4_Full_Biomarker_Performance.csv", row.names = FALSE)
cat(">>> 已保存核心标志物数据至 Table_S4_Full_Biomarker_Performance.csv\n")

# 6. 绘制: AUC 棒棒糖图 (Lollipop Plot)
p_auc_rank <- ggplot(top_biomarkers, aes(x = reorder(Metabolite, AUC), y = AUC, color = Organ)) +
  geom_segment(aes(xend = Metabolite, yend = 0.5), linewidth = 0.5, color = "gray") +
  geom_point(aes(size = abs(Log2FC)), alpha = 0.8) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.5) +
  coord_flip() +
  scale_color_manual(values = c("Serum" = "#E41A1C", "Urine" = "#377EB8")) +
  theme_bw(base_size = 12) +
  labs(title = "Performance of Key Functional Biomarkers",
       subtitle = "Top 20 metabolites ranked by AUC (excluding isotopes)",
       x = "", y = "Area Under the Curve (AUC)", size = "|Log2FC|") +
  theme(panel.grid.minor = element_blank())

# 高清导出
ggsave("Figure_5_AUC_Ranking.pdf", plot = p_auc_rank, width = 8, height = 7, device = "pdf")
ggsave("Figure_5_AUC_Ranking.png", plot = p_auc_rank, width = 8, height = 7, dpi = 300)

print(p_auc_rank)


# ==============================================================================
# 第八部分：系统代谢流重塑分析 (Systemic MCI Flux Analysis)
# 目标：通过 MCI 指数定量证明“底物窃取”与“共生脱毒”机制
# ==============================================================================
library(ggplot2); library(dplyr); library(tidyr); library(ggpubr)

cat(">>> [Step 3.3] 正在执行 MCI 代谢流机制验证分析...\n")

# 1. 锁定核心代谢物 (继承前步优化逻辑)
col_trp  <- get_col_precise("TRYPTOPHAN", exclude = "HYDROXY")
col_kyn  <- get_col_precise("KYNURENINE")
col_put  <- get_col_precise("PUTRESCINE", exclude = "ACETYL")
col_nput <- get_col_precise("N-ACETYLPUTRESCINE")
col_hi   <- get_col_precise("5-HYDROXYINDOLE")
col_ipa  <- get_col_precise("INDOLEPROPIONIC ACID")

# 2. 统一计算比率并整合器官顺序
ratio_full_df <- final_df_processed %>%
  mutate(
    # [宿主路径] IDO 活性 (Kyn / Trp)
    Ratio_IDO = .data[[col_kyn]] - .data[[col_trp]],
    # [细菌脱毒] 乙酰化效率 (N-Ac-Put / Put)
    Ratio_Detox = .data[[col_nput]] - .data[[col_put]],
    # [细菌产物] 综合吲哚产率 ( (IPA+5HI)/2 / Trp )
    Ratio_Indole = (.data[[col_hi]] + .data[[col_ipa]])/2 - .data[[col_trp]]
  ) %>%
  mutate(Organ = factor(Organ, levels = c("Caecal", "Feces", "Serum", "Urine")))

# 3. 提取 8x5 核心统计摘要表 (Table S7)
cat(">>> 正在生成跨器官 MCI 统计汇总表...\n")
summary_mci_table <- ratio_full_df %>%
  group_by(Organ, Grp) %>%
  summarise(across(starts_with("Ratio"), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

write.csv(summary_mci_table, "Table_S7_Systemic_MCI_Trajectory_Stats.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. 绘图 A：机制箱线图 (Figure 8A - 显著性验证)
# ------------------------------------------------------------------------------
p_boxplots_mci <- ratio_full_df %>%
  filter(Organ %in% c("Serum", "Urine")) %>%
  pivot_longer(cols = starts_with("Ratio"), names_to = "Type", values_to = "Val") %>%
  mutate(Type = case_when(
    Type == "Ratio_IDO" ~ "IDO Activity (Kyn/Trp)",
    Type == "Ratio_Detox" ~ "Detox Efficiency (N-Ac-Put/Put)",
    Type == "Ratio_Indole" ~ "Indole Production (Indoles/Trp)"
  )) %>%
  ggplot(aes(x = Grp, y = Val, fill = Grp)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.2) +
  facet_wrap(Organ ~ Type, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "t.test", label = "p.signif", size = 5) +
  scale_fill_manual(values = c("Bt" = "#E41A1C", "GF" = "#377EB8")) +
  theme_bw(base_size = 12) +
  labs(title = "A. Mechanism Validation: Statistical Significance",
       y = "Log2 Ratio (Product / Substrate)", x = "") +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# ------------------------------------------------------------------------------
# 5. 绘图 B：系统轨迹图 (Figure 8B - 跨器官流向)
# ------------------------------------------------------------------------------
p_trajectory_final <- summary_mci_table %>%
  pivot_longer(cols = starts_with("Ratio"), names_to = "Type", values_to = "MCI") %>%
  mutate(Type = case_when(
    Type == "Ratio_IDO" ~ "IDO Activity (Kyn/Trp)",
    Type == "Ratio_Detox" ~ "Detox Efficiency (N-Ac-Put/Put)",
    Type == "Ratio_Indole" ~ "Indole Production (Indole/Trp)"
  )) %>%
  ggplot(aes(x = Organ, y = MCI, color = Grp, group = Grp)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 4, stroke = 1.5, shape = 21, fill = "white") +
  facet_wrap(~Type, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Bt" = "#E41A1C", "GF" = "#377EB8")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  theme_bw(base_size = 12) +
  labs(title = "B. Systemic Metabolic Trajectory (MCI)",
       y = "Mean Conversion Index (MCI)", x = "") +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

# 6. 高清并行导出
ggsave("Figure_8A_MCI_Boxplots.pdf", p_boxplots_mci, width = 10, height = 6)
ggsave("Figure_8A_MCI_Boxplots.tiff", p_boxplots_mci, width = 10, height = 6, dpi = 300, compression = "lzw")
ggsave("Figure_8B_MCI_Trajectory.pdf", p_trajectory_final, width = 8, height = 8)
ggsave("Figure_8B_MCI_Trajectory.tiff", p_trajectory_final, width = 8, height = 8, dpi = 300, compression = "lzw")

print(p_boxplots_mci)
print(p_trajectory_final)


# ==============================================================================
# 第九部分：系统级通路富集分析 (FGSEA)
# ==============================================================================
library(fgsea); library(dplyr); library(ggplot2); library(tidyr)

cat(">>> [Step 3.5] 正在执行系统级通路富集分析 (最终修正版)...\n")

# 1. 准备 Rank 数据 (保持原代码逻辑)
pathway_dat <- final_df_processed %>% filter(Organ == "Urine")
all_mets <- setdiff(colnames(pathway_dat), c("Sample_ID", "Factors", "Join_ID", "Organ", "Grp"))
stats_list <- numeric()

for(m in all_mets) {
  try({
    res <- t.test(pathway_dat[[m]] ~ pathway_dat$Grp)
    stats_list[m] <- res$statistic
  }, silent=TRUE)
}

# 消除并列 (Jittering)
set.seed(42)
stats_list <- stats_list + rnorm(length(stats_list), sd = 1e-10)
ranks <- sort(stats_list, decreasing = TRUE)

# 2. 运行 FGSEA (使用你定义好的 pathways)
fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 2, maxSize = 500)

# 3. 【核心修正】：处理数据以便导出和绘图
# A. 将 leadingEdge 从 List 转为字符串，否则 write.csv 会报错
# B. 创建显式的 Direction 字符列，消除 ggplot 颜色警告
fgsea_final <- as.data.frame(fgseaRes) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", "),
         Direction = ifelse(NES > 0, "Up_in_Bt", "Down_in_Bt")) %>%
  arrange(desc(NES))

# 保存数据表 (Table S9)
write.csv(fgsea_final, "Table_S9_Pathway_Enrichment_Results.csv", row.names = FALSE)

# 4. 绘图：采用最健壮的显式字符映射
p_fgsea_final <- ggplot(fgsea_final, aes(x = reorder(pathway, NES), y = NES, fill = Direction)) +
  geom_col(color = "black", linewidth = 0.3) + # 增加黑边，更有出版感
  coord_flip() +
  # 明确指定颜色，不再依赖逻辑判断
  scale_fill_manual(values = c("Up_in_Bt" = "red", "Down_in_Bt" = "blue"),
                    labels = c("Down in Bt (GF Enrich)", "Up in Bt (Bt Enrich)")) +
  labs(title = "Pathway Enrichment Analysis (Urine)",
       subtitle = "Global Metabolic Shifts",
       x = "", y = "Normalized Enrichment Score (NES)",
       fill = "Enrichment Direction") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(face = "bold"))

# 5. 高保真导出 (PDF + 300 DPI TIFF/PNG)
ggsave("Figure_11_FGSEA_Pathways.pdf", plot = p_fgsea_final, width = 9, height = 6, device = "pdf")
ggsave("Figure_11_FGSEA_Pathways.tiff", plot = p_fgsea_final, width = 9, height = 6, dpi = 300, compression = "lzw")
ggsave("Figure_11_FGSEA_Pathways.png", plot = p_fgsea_final, width = 9, height = 6, dpi = 300)

# 打印显示
print(p_fgsea_final)

cat(">>> [Step 3.5] 全部工作圆满完成！\n")
