# ============================================================
# Lattice 6x6 Analysis Script (With CV% Added)
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
# تابع ساخت جدول دقیق تجزیه واریانس لاتیس (با CV%)
# ============================================================
generate_lattice_anova <- function(y, rep, block, trt, trait_name) {
  
  # ساخت دیتافریم داخلی
  df_model <- data.frame(
    Y = as.numeric(y),
    Rep = as.factor(rep),
    Block = as.factor(block),
    Trt = as.factor(trt)
  )
  
  # میانگین کل صفت
  grand_mean <- mean(df_model$Y, na.rm = TRUE)
  
  # ساخت بلوک آشیانه‌شده درون تکرار
  df_model$RepBlock <- interaction(df_model$Rep, df_model$Block, drop = TRUE)
  
  # مدل 1: RCBD
  m_unadj <- lm(Y ~ Rep + Trt, data = df_model)
  a_unadj <- anova(m_unadj)
  
  # مدل 2: Lattice (با بلوک آشیانه‌شده)
  m_adj <- lm(Y ~ Rep + RepBlock + Trt, data = df_model)
  a_adj <- anova(m_adj)
  
  v <- length(unique(df_model$Trt))
  r <- length(unique(df_model$Rep))
  
  # استخراج درجات آزادی (df)
  df_rep <- a_adj["Rep", "Df"]
  df_unadj_trt <- a_unadj["Trt", "Df"]
  df_adj_trt <- a_adj["Trt", "Df"]
  df_block <- a_adj["RepBlock", "Df"]  
  df_rcbd_err <- a_unadj["Residuals", "Df"] 
  df_intra_err <- a_adj["Residuals", "Df"]  
  
  # استخراج میانگین مربعات (MS)
  ms_rep <- a_adj["Rep", "Mean Sq"]
  ms_unadj_trt <- a_unadj["Trt", "Mean Sq"]
  ms_adj_trt <- a_adj["Trt", "Mean Sq"]
  ms_block <- a_adj["RepBlock", "Mean Sq"]
  ms_rcbd_err <- a_unadj["Residuals", "Mean Sq"]
  ms_intra_err <- a_adj["Residuals", "Mean Sq"]
  
  # محاسبه CV% (بر اساس خطای داخل بلوکی و میانگین کل)
  cv_val <- (sqrt(ms_intra_err) / grand_mean) * 100
  cv_str <- sprintf("%.2f", cv_val) # فرمت کردن به 2 رقم اعشار
  
  # محاسبه خطای موثر
  loss <- ((v-1) * (ms_rcbd_err - ms_intra_err)) / ms_rcbd_err
  ms_eff_err <- ms_intra_err * (v-1) / (v-1 - loss)
  if(is.na(ms_eff_err) || ms_eff_err == ms_rcbd_err) ms_eff_err <- ms_intra_err
  
  # استخراج P-value
  p_unadj_trt <- a_unadj["Trt", "Pr(>F)"]
  p_adj_trt <- a_adj["Trt", "Pr(>F)"]
  p_block <- a_adj["RepBlock", "Pr(>F)"]
  
  get_stars <- function(p) {
    if(is.na(p)) return("")
    if(p < 0.001) return("***")
    if(p < 0.01) return("**")
    if(p < 0.05) return("*")
    if(p < 0.1) return(".")
    return("ns")
  }
  
  # ساخت جدول نهایی (اضافه شدن ستون CV)
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
    CV_percent = c("", "", "", "", "", "", cv_str), # قرار دادن CV در ردیف خطای داخل بلوکی
    Sig = c("", 
            get_stars(p_unadj_trt), 
            get_stars(p_adj_trt), 
            get_stars(p_block), 
            "", "", ""),
    stringsAsFactors = FALSE
  )
  
  colnames(res_table)[3] <- paste0(trait_name, " (MS)")
  colnames(res_table)[4] <- "CV (%)"
  
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

cat("\n✔ تمام جدول‌های ANOVA کامل به همراه CV% در پوشه", work_dir, "ذخیره شدند.\n")