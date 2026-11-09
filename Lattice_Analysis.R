# ============================================================
# Lattice 6x6 Analysis Script (Exact ANOVA Table like SAS)
# ============================================================

if (!require("agricolae")) install.packages("agricolae")
library(agricolae)

# 1. تنظیم مسیر و خواندن داده‌ها
work_dir <- "F:/Lattice/Chikpea"
setwd(work_dir)
lattice_data <- read.csv("Data.csv", stringsAsFactors = FALSE)

lattice_data$Rep <- as.factor(lattice_data$Rep)
lattice_data$Block <- as.factor(lattice_data$Block)
lattice_data$ENTRY <- as.factor(lattice_data$ENTRY)

traits <- c("VGA", "DHE", "DMA", "PLH", "DTP", "TKW", 
            "Resistant", "YLD", "SiP", "NMB", "NSB", "HFP", "PIP")

# ============================================================
# تابع ساخت جدول دقیق تجزیه واریانس لاتیس (مطابق عکس)
# ============================================================
generate_lattice_anova <- function(y, rep, block, trt, trait_name) {
  
  # مدل 1: RCBD (بدون بلوک) -> برای تیمار تصحیح نشده و خطای RCBD
  m_unadj <- lm(y ~ rep + trt)
  a_unadj <- anova(m_unadj)
  
  # مدل 2: Lattice (با بلوک) -> برای تیمار تصحیح شده، بلوک و خطای داخل بلوک
  m_adj <- lm(y ~ rep + block + trt)
  a_adj <- anova(m_adj)
  
  v <- length(unique(trt))
  r <- length(unique(rep))
  k <- length(unique(block)) / r
  
  # استخراج درجات آزادی (df)
  df_rep <- a_adj["rep", "Df"]
  df_unadj_trt <- a_unadj["trt", "Df"]
  df_adj_trt <- a_adj["trt", "Df"]
  df_block <- a_adj["block", "Df"]
  df_rcbd_err <- a_unadj["Residuals", "Df"]
  df_intra_err <- a_adj["Residuals", "Df"]
  
  # استخراج میانگین مربعات (MS)
  ms_rep <- a_adj["rep", "Mean Sq"]
  ms_unadj_trt <- a_unadj["trt", "Mean Sq"]
  ms_adj_trt <- a_adj["trt", "Mean Sq"]
  ms_block <- a_adj["block", "Mean Sq"]
  ms_rcbd_err <- a_unadj["Residuals", "Mean Sq"]
  ms_intra_err <- a_adj["Residuals", "Mean Sq"]
  
  # محاسبه خطای موثر (Effective Error) - فرمول کاکرین و کاکس
  loss <- ((v-1) * (ms_rcbd_err - ms_intra_err)) / ms_rcbd_err
  ms_eff_err <- ms_intra_err * (v-1) / (v-1 - loss)
  if(is.na(ms_eff_err) || ms_eff_err == ms_rcbd_err) ms_eff_err <- ms_intra_err
  
  # استخراج P-value برای ستاره‌گذاری
  p_unadj_trt <- a_unadj["trt", "Pr(>F)"]
  p_adj_trt <- a_adj["trt", "Pr(>F)"]
  p_block <- a_adj["block", "Pr(>F)"]
  
  get_stars <- function(p) {
    if(is.na(p)) return("")
    if(p < 0.001) return("***")
    if(p < 0.01) return("**")
    if(p < 0.05) return("*")
    if(p < 0.1) return(".")
    return("ns")
  }
  
  # ساخت جدول نهایی
  res_table <- data.frame(
    Source = c("Replication", 
               "Treatment (unadjusted)", 
               "Treatment (adjusted)", 
               "Block within Replication", 
               "Error (Effective)", 
               "RCBD Error", 
               "Intra-block Error"),
    df = c(df_rep, df_unadj_trt, df_adj_trt, df_block, df_intra_err, df_rcbd_err, df_intra_err),
    Mean_Sq = round(c(ms_rep, ms_unadj_trt, ms_adj_trt, ms_block, ms_eff_err, ms_rcbd_err, ms_intra_err), 3),
    Sig = c("", 
            get_stars(p_unadj_trt), 
            get_stars(p_adj_trt), 
            get_stars(p_block), 
            "", "", ""),
    stringsAsFactors = FALSE
  )
  
  # نام‌گذاری ستون MS به اسم صفت
  colnames(res_table)[3] <- paste0(trait_name, " (MS)")
  return(res_table)
}

# ============================================================
# اجرای حلقه برای تمام صفات
# ============================================================
for (trait in traits) {
  cat("\n\n==================================================\n")
  cat("ANOVA Table for Trait:", trait, "\n")
  cat("==================================================\n")
  
  # حذف NA ها
  y_vec <- suppressWarnings(as.numeric(lattice_data[[trait]]))
  keep <- !is.na(y_vec)
  df_sub <- lattice_data[keep, ]
  y_sub <- y_vec[keep]
  
  # تولید جدول
  anova_table <- generate_lattice_anova(
    y = y_sub, 
    rep = df_sub$Rep, 
    block = df_sub$Block, 
    trt = df_sub$ENTRY, 
    trait_name = trait
  )
  
  # چاپ در کنسول
  print(anova_table, row.names = FALSE)
  
  # ذخیره به صورت فایل CSV
  write.csv(anova_table, file.path(work_dir, paste0(trait, "_Full_ANOVA.csv")), row.names = FALSE)
}

cat("\n✔ تمام جدول‌های ANOVA کامل در پوشه", work_dir, "ذخیره شدند.\n")