# Data Analysis in R

 **In reference to : Lentil-Analysis_2026.R file**
---

# Lentil Diversity Screen Analysis in R

**From raw phenotyping data to statistics, visualisation, and clustering**

This tutorial walks through a complete R workflow for analysing a lentil diversity screen under **high and low nitrogen treatments**, covering:

* Data cleaning and genotype harmonisation
* Summary statistics (mean, SD, SEM)
* Bar plots and scatter plots with biologically meaningful highlighting
* Trait–trait relationships (biomass vs nodulation, %NDFA)
* Heatmap visualisation and genotype clustering

The workflow is designed for **reproducible analysis** and **publication-quality figures**.

---

## 1. Setup and Libraries

We rely on the tidyverse ecosystem for data handling, `ggplot2` for visualisation, and `ComplexHeatmap` for clustering.

```r
library(tidyverse)
library(ggplot2)
library(readxl)
library(writexl)
library(ComplexHeatmap)
library(circlize)
```

---

## 2. Data Import and Initial Filtering

Raw data are read from Excel and filtered to retain key traits and metadata:

* **Genotype**
* **Nitrogen treatment (High / Low)**
* **Shoot dry weight**
* **Active nodules**
* **%NDFA**

```r
raw_lent <- read_xlsx("Lentil_Final_v2.xlsx", sheet = "all data")

raw_filt <- raw_lent[, c(11, 15, 20, 21, 22, 33)] %>%
  rename(Genotype = Genotype...15) %>%
  mutate(
    Genotype = factor(Genotype),
    Trt = factor(Trt),
    `Dry Weight (g)` = as.numeric(`Dry Weight (g)`),
    Count_active = as.numeric(Count_active),
    `%Ndfa` = as.numeric(`%Ndfa`)
  )
```

---

## 3. Genotype Name Harmonisation

Some genotypes are recorded with replicate suffixes (e.g. `Sym9_1`, `Sym9_3`) or multiple internal IDs. These are collapsed into biologically meaningful groups.

```r
raw_filt <- raw_filt %>%
  mutate(
    Genotype_Name = str_replace(Genotype_Name, "(Sym9|Sym19).*", "\\1"),
    Genotype = case_when(
      Genotype %in% paste0("G", 307:313) ~ "Sym9",
      Genotype %in% paste0("G", 314:320) ~ "Sym19",
      TRUE ~ Genotype
    )
  )
```

---

## 4. Summary Statistics by Genotype and Treatment

Traits are reshaped to long format and summarised to calculate:

* Mean
* Standard deviation
* Sample size
* Standard error of the mean (SEM)

```r
long_dat <- reshape2::melt(
  raw_filt,
  id.vars = c("Genotype_Name", "Genotype", "Trt"),
  na.rm = TRUE
)

summary_stats <- long_dat %>%
  group_by(Genotype_Name, Genotype, Trt, variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    sem = sd / sqrt(n),
    .groups = "drop"
  )
```

These statistics are saved for downstream use and transparency.

---

## 5. Sorting Genotypes by Performance

To aid interpretation, genotypes are ordered by **biomass under low or high nitrogen**, then reused consistently across plots.

```r
sorted_lowN <- summary_stats %>%
  filter(Trt == "Low", variable == "Dry Weight (g)") %>%
  arrange(desc(mean)) %>%
  pull(Genotype)
```

---

## 6. Biomass Bar Plots (High vs Low N)

Bar plots show genotype means ± SEM, either:

* **Combined** (both treatments together), or
* **Faceted** by nitrogen treatment

```r
ggplot(
  summary_stats %>% filter(variable == "Dry Weight (g)"),
  aes(x = Genotype, y = mean, fill = Trt)
) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  theme_bw() +
  labs(y = "Biomass (g DW)", x = "Genotype")
```

---

## 7. Comparing Performance Across Nitrogen Levels

Data are reshaped to wide format to directly compare **Low vs High N biomass** per genotype.

```r
biomass <- summary_stats %>%
  filter(variable == "Dry Weight (g)") %>%
  select(Genotype, Trt, mean) %>%
  pivot_wider(names_from = Trt, values_from = mean)
```

A scatter plot highlights **controls and commercial cultivars** while keeping all other genotypes visible.

---

## 8. Trait–Trait Relationships

### Biomass vs Active Nodules

### Biomass vs %NDFA

Relationships are visualised separately for High and Low N using faceting, with linear models overlaid.

```r
ggplot(lentils_wide,
       aes(x = mean_Dry_Weight, y = mean_Count_active)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~Trt) +
  theme_bw()
```

This allows rapid assessment of treatment-specific trait coupling.

---

## 9. Heatmap and Genotype Clustering

### Z-score Scaling

Biomass values are scaled **column-wise** (within treatment) to highlight relative genotype performance.

```r
mat <- biomass %>% select(Low, High) %>% as.matrix()
rownames(mat) <- biomass$Genotype

mat_scaled <- scale(mat)
```

---

### Heatmap Generation

Genotypes are clustered using Ward’s method, with treatments shown as column groups.

```r
Heatmap(
  mat_scaled,
  name = "Biomass (z-score)",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_columns = FALSE,
  row_km = 16,
  column_split = c("Low N", "High N")
)
```

---

## 10. Extracting Clusters

Cluster membership is extracted directly from the heatmap object and exported for downstream biological interpretation.

```r
clusters <- row_order(heatmap_object)

cluster_df <- map_df(
  names(clusters),
  ~data.frame(
    Genotype = rownames(mat)[clusters[[.x]]],
    Cluster = paste0("Cluster_", .x)
  )
)

write_xlsx(cluster_df, "Biomass_clusters.xlsx")
```

---

## Summary

This workflow demonstrates how to:

* Clean and harmonise complex genotype datasets
* Generate reproducible summary statistics
* Create interpretable, publication-quality figures
* Explore trait relationships under contrasting environments
* Identify genotype performance clusters using heatmaps

The approach is modular and can be readily adapted to other crops, traits, or experimental designs.

---


