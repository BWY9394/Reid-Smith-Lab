library("gridExtra")
library("PCAtools")
library("scales")
library("ggplot2")
library("plyr")
library(tidyr)
library("ggfortify")
library(dplyr)
library(gridExtra)
library(stringr)
library(ggtext)
#stats
library("multcomp")
library(multcompView)
library(emmeans)


#wd
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Consultation/")
setwd("C:/Users/BWeeY/OneDrive - LA TROBE UNIVERSITY/Consultation/")

#load in datasets
raw_lent <- readxl::read_xlsx("Meysam/Lentil_Final_v2.xlsx",sheet="all data")
colnames(raw_lent)
#Need DW, NDFA, Act Nods, for high and low N for factors of Genotype (Genotype Name +Genotype) and Trt (High/Low)
raw_filt_lent <- raw_lent[,c(11,15,20,21,22,33)]
raw_filt_lent$`Dry Weight (g)` <-  na_if(raw_filt_lent$`Dry Weight (g)`, "") # a bit pointless as the as.numeric functions coerces characters into NAs
raw_filt_lent <- dplyr::rename(raw_filt_lent, Genotype = Genotype...15)
raw_filt_lent$Genotype <- as.factor(raw_filt_lent$Genotype)
raw_filt_lent$Trt <- as.factor(raw_filt_lent$Trt)
raw_filt_lent$`Dry Weight (g)` <- as.numeric(raw_filt_lent$`Dry Weight (g)`)
raw_filt_lent$Count_active <- as.numeric(raw_filt_lent$Count_active )
raw_filt_lent$`%Ndfa` <- as.numeric(raw_filt_lent$`%Ndfa` )


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

colnames(raw_filt_lent_symadj)
#do stats
combined_melt= reshape2::melt(raw_filt_lent_symadj,id.vars=c("Genotype_Name","Genotype","Trt"),na.rm=TRUE)#na.rm removes NA values, groups by the stated factors
combined_melt$value <- as.numeric(combined_melt$value)
#stats, mean, sd, se, n
control_lentils_average_stats <- plyr::ddply(
  combined_melt,
  c("Genotype_Name", "Genotype","Trt","variable"),
  summarise,
  mean = mean(value, na.rm = TRUE),
  sd   = sd(value, na.rm = TRUE),
  n    = sum(!is.na(value)),
  sem  = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
)



#Great save now
writexl::write_xlsx(
  control_lentils_average_stats,
  "Meysam/2026/Lentil_Avge_Stats_2026.xlsx"
)
control_lentils_average_stats <- readxl::read_xlsx("Meysam/2026/Lentil_Avge_Stats_2026.xlsx")
str(control_lentils_average_stats)
control_lentils_average_stats$Genotype <- as.factor(control_lentils_average_stats$Genotype)
control_lentils_average_stats$Trt <- as.factor(control_lentils_average_stats$Trt)
control_lentils_average_stats$variable <- as.factor(control_lentils_average_stats$variable)


#Ok, now do some sorting
levels(control_lentils_average_stats$Trt)
levels(control_lentils_average_stats$variable)


lentil_average_stats_sortHigh  <- control_lentils_average_stats %>%
  dplyr::mutate(Genotype = factor(Genotype,
                                         levels = control_lentils_average_stats %>%
                                           filter(Trt == "High", variable == "Dry Weight (g)") %>%
                                           arrange(desc(mean)) %>%
                                           pull(Genotype)
  )
  )

#not working very well
#lentil_average_stats_sortHigh <- lentil_average_stats_sortHigh %>%
 # mutate(
  #  Genotype_label = case_when(
   #   Genotype_Name == "G228"  ~ "<span style='color:red'>G228</span>",
    #  Genotype_Name == "Sym9"  ~ "<span style='color:blue'>Sym9</span>",
     # Genotype_Name == "Sym19" ~ "<span style='color:purple'>Sym19</span>",
      #TRUE ~ as.character(Genotype)
  #  )
#  )

#View(lentil_average_stats_sortHigh)

lentil_average_stats_sortLow  <- control_lentils_average_stats %>%
  dplyr::mutate(Genotype = factor(Genotype,
                                         levels = control_lentils_average_stats %>%
                                           filter(Trt == "Low", variable == "Dry Weight (g)") %>%
                                           arrange(desc(mean)) %>%
                                           pull(Genotype)
  )
  )


#now plot barplots

#First let's do the low N sorted bargraphs
#combined treatments, sorted by low N in one graph
Combined_lowN_biomass_lowN_sort <- ggplot(subset(lentil_average_stats_sortLow, variable %in% "Dry Weight (g)"), aes(x = Genotype, y = mean, fill = Trt)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  #facet_wrap(~ Trt, ncol = 1, scales = "fixed") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Lentil Biomass sorted by low N") +
  scale_fill_brewer(palette = "Dark2") +  # Colorblind-friendly
  theme_bw(base_size = 12) +
  # ylim(NA,NA)+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(size = 2.5, face = "plain", angle=80,hjust=1,vjust=0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title = element_text(hjust = 0.2, size = 14, face = "bold")
  )

Combined_lowN_biomass_lowN_sort

#combined treatments, sorted by low N in one graph, faceted into 2 graphs
lowN_biomass_lowN_sort_facet <- ggplot(subset(lentil_average_stats_sortLow, variable %in% "Dry Weight (g)"), aes(x = Genotype, y = mean, fill = Trt)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  facet_wrap(~ Trt, ncol = 1, scales = "free",labeller = "label_both") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Lentil Biomass sorted by low N, combined") +
  scale_fill_brewer(palette = "Dark2") +  # Colorblind-friendly
  theme_bw(base_size = 12) +
  # ylim(NA,NA)+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(size = 2.5, face = "plain", angle=80,hjust=1,vjust=0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title = element_text(hjust = 0.2, size = 14, face = "bold")
  )

lowN_biomass_lowN_sort_facet
str(lentil_average_stats_sortHigh)
# Now let's do the high N sorted bargraphs
# combined treatments, sorted by high N in one graph
Combined_highN_biomass_HighN_sort <- ggplot(subset(lentil_average_stats_sortHigh, variable %in% "Dry Weight (g)"), aes(x = Genotype, y = mean, fill = Trt)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  #facet_wrap(~ Trt, ncol = 1, scales = "fixed") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Lentil Biomass sorted by high N") +
  scale_fill_brewer(palette = "Dark2") + # Colorblind-friendly
  theme_bw(base_size = 12) +
  # ylim(NA,NA)+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(size = 2.0, face = "plain", angle=80,hjust=1,vjust=0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title = element_text(hjust = 0.2, size = 14, face = "bold")
  )

Combined_highN_biomass_HighN_sort

#combined treatments, sorted by high N in one graph, faceted into 2 graphs
highN_biomass_highN_sort_facet <- ggplot(subset(lentil_average_stats_sortHigh, variable %in% "Dry Weight (g)"), aes(x = Genotype, y = mean, fill = Trt)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  facet_wrap(~ Trt, ncol = 1, scales = "free",labeller = "label_both") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Lentil Biomass sorted by high N, combined") +
  scale_fill_brewer(palette = "Dark2") +  # Colorblind-friendly
  theme_bw(base_size = 12) +
  # ylim(NA,NA)+
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(size = 2.5, face = "plain", angle=80,hjust=1,vjust=0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title = element_text(hjust = 0.2, size = 14, face = "bold")
  )

highN_biomass_highN_sort_facet
ggsave(Combined_highN_biomass_HighN_sort,file="./Meysam/2026/testcolors.pdf",units="cm",height=24,width=25)



#Ok great, reproduced, now moving on.
str(control_lentils_average_stats)
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


#first let's explore the data
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
ggsave(Scatter_biomass_prelim, file="./Meysam/2026/Scatterplot_Lent_all_labels.pdf",units="cm",height=24,width=25)
getwd()
#great, now let's say we want to highlight individual cultivars
biomass <- biomass %>%
  mutate(
    `Commercial Cultivars/Peas` = case_when(
      Genotype %in% c("Sym9", "Sym19") ~ "Control",
      Genotype %in% c("G300")          ~ "Commercial",
      TRUE                             ~ "Other"
    )
  )
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

Scatter_biomass

ggsave(Scatter_biomass,file="./Meysam/2026/Scatterplot_Lent.pdf",units="cm",height=25,width=30)

#Ok, now let's do what meysam requested. DW vs Act Nods
#now do a preview plot
ggplot(lentils_wide,  aes(x = `mean_Dry Weight (g)`, y = mean_Count_active)) +
  geom_point(color = "black", size = 1.5) +
  geom_text(aes(label = Genotype), color = "black", size = 3, vjust = -0.5) +
  facet_wrap(~ Trt) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    strip.text   = element_text(size = 11, face = "bold", hjust = 0),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Lentil diversity screen: Shoot dry weight (g) vs. active nodules for High and Low N",
    x = "Mean shoot dry weight (g)",
    y = "Mean active nodules"
  )
#Decide what genotype numbers you want to retain
#Then, add to the below code
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
ggsave(DW_vsactNods,file="./Meysam/2026/Lentils_DW_vs_ActNods.pdf",units="cm",height=20,width=30)


# DW vs NDFA
#now do a preview plot
ggplot(lentils_wide,  aes(x = `mean_Dry Weight (g)`, y = `mean_%Ndfa`)) +
  geom_point(color = "black", size = 1.5) +
  geom_text(aes(label = Genotype), color = "black", size = 3, vjust = -0.5) +
  facet_wrap(~ Trt) +
  theme_bw(base_size = 12) +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    strip.text   = element_text(size = 11, face = "bold", hjust = 0),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Lentil diversity screen: Shoot dry weight (g) vs. active nodules for High and Low N",
    x = "Mean shoot dry weight (g)",
    y = "Mean active nodules"
  )

#Decide what genotype numbers you want to retain
#Then, add to the below code
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
ggsave(DW_vs_NDFA,file="./Meysam/2026/Lentils_DW_vs_NDFA.pdf",units="cm",height=20,width=30)


#Ok, time for heatmaps
View(biomass)
str(biomass)
biomassfilt <- biomass[c(1,3,4)]
str(biomassfilt)

#prepare matrix
mat <- biomassfilt %>%
 dplyr:: select(High, Low) %>%
  as.matrix()

rownames(mat) <- biomassfilt$Genotype

# Z-score per genotype (row-wise)
mat_scale <- as.data.frame(apply(mat, 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))) #column wise scoring, better for comparing between genotypes per Treatment.
mat_scale <- as.matrix(mat_scale)

library(ComplexHeatmap)
library(circlize)
colfun2 = circlize::colorRamp2(c(-2,0,2), #gives range of the scale
                               c("#408ea4","lightyellow1","#c5242a")) #for zscoring
col_order <- data.frame(
  Column = colnames(mat_scaled),         # e.g., "Low", "High"
  Grouping = c("Low N", "High N"),         # your desired group labels
  stringsAsFactors = FALSE
)
View(col_order)

dev.off()
set.seed(123)
pdf("Biomass_zscore_HM.pdf", height=10, width =8)
hf=Heatmap(
  mat_scale,
  use_raster=TRUE,
  name = "Biomass\n(z-score)",
  col = colfun2,
  #cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_km = 16,
  row_km_repeats = 100,
  cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST. ALSO HAVE IT OFF TO EXTRACT CLUSTERS VERY IMPORTANT
  show_row_names = TRUE,
  #row_split = row_km, #split by cluster numbe
  clustering_method_rows = "ward.D",
  row_dend_reorder = TRUE,
  row_names_gp = gpar(fontsize = 2),
  row_title_gp = gpar(font = 2,fontsize=6),
  row_gap = unit(2, "mm"),
  column_split = factor(col_order$Grouping),
  column_title = ("Nitrogen treatment"),
  #row_title = "Genotypes"
 # show_row_dend = TRUE
)
hf
HM=draw(hf)
dev.off()
library(magrittr)
#extract clusters

r.dend <- row_dend(HM)  #If needed, extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x)) #check/confirm size gene clusters
# Extract clusters as a data.frame
clusterlist = row_order(mat_scale)

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(geneID = rownames(mat_scale[rcl.list[[i]],]), #rownames(tf.log)[clusterlist[[i]]], #if cluster somehow only has one gene. https://www.biostars.org/p/465304/
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

View(clu_df)
writexl::write_xlsx(
  clu_df,
  "Meysam/2026/clustered z-score biomass (per Treatment across genotypes).xlsx"
)

