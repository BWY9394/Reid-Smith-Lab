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
#more efficient way to group rename from same DF, code does effectively the same thing
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

Some genotypes are recorded with replicate suffixes (e.g. `Sym9_1`, `Sym9_3`) or multiple internal IDs. These are collapsed into biologically meaningful groups. E.g.

```r
#First let's make clean up all the Sym9_1 _3 etc into just Sym9 and the like
raw_filt_lent_symadj <- raw_filt_lent %>%
  mutate(Genotype_Name = str_replace(Genotype_Name,
                                     "(Sym9|Sym19).*",
                                     "\\1")
         )

#Additional problem, he has multiple genotypes named per Genotype_Name, i.e Sym9 can be G
raw_filt_lent_symadj <- raw_filt_lent_symadj %>%
  mutate(    Genotype = str_replace(Genotype,
                                    "(G307|G308|G309|G310|G311|G312|G313).*",
                                    "Sym9")
        )

raw_filt_lent_symadj <- raw_filt_lent_symadj %>%
  mutate(    Genotype = str_replace(Genotype,
                                    "(G314|G315|G316|G317|G318|G319|G320).*",
                                    "Sym19")
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
combined_melt= reshape2::melt(raw_filt_lent_symadj,id.vars=c("Genotype_Name","Genotype","Trt"),na.rm=TRUE)#na.rm removes NA values, groups by the stated factors
combined_melt$value <- as.numeric(combined_melt$value)

control_lentils_average_stats <- plyr::ddply(
  combined_melt,
  c("Genotype_Name", "Genotype","Trt","variable"),
  summarise,
  mean = mean(value, na.rm = TRUE),
  sd   = sd(value, na.rm = TRUE),
  n    = sum(!is.na(value)),
  sem  = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
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
lentils_wide <- control_lentils_average_stats %>%
  pivot_wider(
    id_cols = c(Genotype_Name, Genotype, Trt), #id_cols → the identifiers you want to keep as rows (Genotype.number, Genotype.name, N2.Treatment).
    names_from = variable, #names_from = variable → the factor column whose levels become new column names.
    values_from = c(mean, sd, n, sem) #values_from = c(mean, sd, n, sem) → the measurement columns that get spread into those new variable-specific columns.
  )

#Biomass
biomass <- control_lentils_average_stats %>%
  filter(variable == "Dry Weight (g)") %>%
  dplyr::select(Genotype, Genotype_Name, Trt, mean) %>%
  pivot_wider(
    names_from = Trt,
    values_from = mean
  )

```

```r
Scatter_biomass_prelim <- ggplot(biomass, aes(x = Low, y = High)) +
  geom_point(color = "black") +
  geom_text(aes(label = Genotype), color="black",check_overlap = FALSE, size = 3, vjust = -0.5) +
  #scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size=12) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x  = element_text(size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    strip.text   = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    #legend.position = "none" # <-- hides the legend  # <-- hides the legend
  ) +
  labs(title = "Lentil diversity screen: Shoot dry weight (g)")
  Scatter_biomass_prelim #View it
```
Next, generate some groupings so that you can highlight controls and commercial cultivars, as well as any other genotypes you may find interesting

```r
biomass <- biomass %>%
  mutate(
    `Commercial Cultivars/Peas` = case_when(
      Genotype %in% c("Sym9", "Sym19") ~ "Control",
      #Genotype %in% c("GWHATEVER")     ~ "High N-responsive and low N resilient #examle for editting in groupings for interesting genotypes.
      Genotype %in% c("G300")          ~ "Commercial",
      TRUE                             ~ "Other"
    )
  )

```


Now, we then generate a scatter plot highlighting **controls and commercial cultivars** while keeping all other genotypes visible.

```r
# now you can plot
Scatter_biomass <- ggplot(biomass, aes(x = Low, y = High)) +
  geom_point(aes(color = `Commercial Cultivars/Peas`)) +
  # label ONLY Control + Commercial
  geom_text(
    data = subset(biomass,`Commercial Cultivars/Peas` %in% c("Control", "Commercial")),
    aes(label = Genotype, color = `Commercial Cultivars/Peas`),
    size = 3, vjust = -0.5) +
  #p rovide color scale
  scale_color_manual(
    values = c("Control" = "blue",
               "Commercial" = "red",
               "Other" = "black") ) +
    theme_bw(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(title = "Lentil diversity screen: Shoot dry weight (g)")
Scatter_biomass #View plot
```

---

## 8. Trait–Trait Relationships

### Biomass vs Active Nodules
```r
lentils_wide <- lentils_wide %>%
  mutate(
    `Commercial Cultivars/Peas` = case_when(
      Genotype %in% c("Sym9", "Sym19") ~ "Control",
      Genotype %in% c("G300")          ~ "Commercial",
      TRUE                             ~ "Other"
    )
  )

#Now, do scatter plots for DW vs NDFA
colnames(lentils_wide)
DW_vsactNods <- ggplot(lentils_wide,aes(x = `mean_Dry Weight (g)`, y = mean_Count_active)) +
  geom_point(aes(color = `Commercial Cultivars/Peas`)) +
  # label ONLY Control + Commercial
  geom_text(
    data = subset(lentils_wide,`Commercial Cultivars/Peas` %in% c("Control", "Commercial")),
    aes(label = Genotype, color = `Commercial Cultivars/Peas`),
    size = 3, vjust = -0.5) +
  #p rovide color scale
  scale_color_manual(
    values = c("Control" = "blue",
               "Commercial" = "red",
               "Other" = "black") ) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm", se = TRUE,color="orange",fill="orange",linewidth=0.5) +
  facet_wrap(~Trt) +
  # ylim(-5,NA)+
  #xlim(-5,NA)+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x  = element_text(size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    strip.text   = element_text(size = 11, face = "bold",hjust =0.01),
    strip.background = element_blank(),#element_rect(fill = "gray90", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    #legend.position = "none"   # <-- hides the legend
  ) +
  labs(title = "Lentil diversity screen: Shoot dry weight (g) vs. active nodules for High and Low N")

DW_vsactNods

```

### Biomass vs %NDFA
```r
DW_vs_NDFA <- ggplot(lentils_wide,aes(x = `mean_Dry Weight (g)`, y = `mean_%Ndfa`)) +
  geom_point(aes(color = `Commercial Cultivars/Peas`)) +
  # label ONLY Control + Commercial
  geom_text(
    data = subset(lentils_wide,`Commercial Cultivars/Peas` %in% c("Control", "Commercial")),
    aes(label = Genotype, color = `Commercial Cultivars/Peas`),
    size = 3, vjust = -0.5) +
  #p rovide color scale
  scale_color_manual(
    values = c("Control" = "blue",
               "Commercial" = "red",
               "Other" = "black") ) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm", se = TRUE,color="orange",fill="orange",linewidth=0.5) +
  facet_wrap(~Trt) +
  # ylim(-5,NA)+
  #xlim(-5,NA)+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x  = element_text(size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    strip.text   = element_text(size = 11, face = "bold",hjust =0.01),
    strip.background = element_blank(),#element_rect(fill = "gray90", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    #legend.position = "none"   # <-- hides the legend
  ) +
  labs(title = "Lentil diversity screen: Shoot dry weight (g) vs. %NDFA for High and Low N")

DW_vs_NDFA

```

Relationships are visualised separately for High and Low N using faceting, with linear models overlaid. Example code trend:

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

Interpreting clusters:

| Cluster type         | Biomass pattern across N treatments                    | Biological interpretation                                                                | Potential use                                                       |
| -------------------- | ------------------------------------------------------ | ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| High–High            | High biomass under both Low and High N                 | Genotypes with strong overall growth and stable performance across nitrogen availability | Candidate parents for breeding broadly adapted, high-yielding lines |
| Low–Low              | Low biomass under both Low and High N                  | Generally poor performers, possibly limited by growth capacity or symbiotic efficiency   | Lower priority lines; useful as contrasts in physiological studies  |
| Low–High             | Low biomass under Low N, high biomass under High N     | Nitrogen-responsive genotypes that rely strongly on external N supply                    | Useful for identifying traits linked to fertiliser responsiveness   |
| High–Low             | High biomass under Low N, reduced biomass under High N | Potentially nitrogen-efficient or symbiosis-optimised genotypes                          | High interest for low-input or sustainable systems                  |
| Intermediate / mixed | Moderate or variable biomass across treatments         | Genotypes with context-dependent or plastic responses to nitrogen                        | Targets for studying genotype × environment interactions            |




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


