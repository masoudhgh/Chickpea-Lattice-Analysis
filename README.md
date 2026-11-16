# Chickpea Lattice Analysis (طرح لاتیس نخود)

This repository contains R scripts for analyzing a 6x6 simple lattice design (طرح لاتیس ساده ۶×۶) for Chickpea traits. 
این مخزن شامل کدهای R برای تجزیه و تحلیل طرح لاتیس ساده ۶×۶ برای صفات نخود است.

## 📌 Features (ویژگی‌ها)
- **Nested Block Design:** Correctly models `Rep/Block` to estimate lattice ANOVA.
- **Comprehensive Output:** Calculates Mean Squares (MS), F-values, and p-values for all traits.
- **CV% Calculation:** Automatically calculates the Coefficient of Variation based on the intra-block error.
- **Wide Format Export:** Combines all traits into a single, easy-to-read CSV table (`Combined_ANOVA_Table.csv`).

## 🛠️ Requirements (پیش‌نیازها)
You need R installed on your system with the following packages:
```R
install.packages(c("agricolae", "tidyr", "reshape2"))
