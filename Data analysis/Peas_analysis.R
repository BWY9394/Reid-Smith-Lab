# Set WD
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Consultation/")
setwd("C:/Users/BWeeY/OneDrive - LA TROBE UNIVERSITY/Consultation/")

#load required packages
#PCA Analysis and plotting
#BiocManager::install("pcaExplorer")
library("pcaExplorer")
library("gridExtra")
library("PCAtools")
library("scales")
library("ggplot2")
#library("sleuth")
library("plyr")
library(tidyr)
library("ggfortify")
library(dplyr)
library(gridExtra)
library(stringr)
#stats
library("multcomp")
library(multcompView)
library(emmeans)


#load in data, clean and verify
lentil_biomass_HighN <- read.csv("Meysam/GRDC milestone 105_Pea_Nov 2025_HighN.csv",header=TRUE)
lentil_biomass_lowN <- read.csv("Meysam/GRDC milestone 105_Pea_Nov 2025_LowN.csv",header=TRUE)
colnames(lentil_biomass_HighN)
colnames(lentil_biomass_lowN)

combined <- rbind(lentil_biomass_HighN,lentil_biomass_lowN)
head(combined$Genotype.number)

# Add the prefix to the 'Name' column
combined$Genotype.number <- paste0("P", combined$Genotype.number)
head(combined$Genotype.number)

#Make factor
combined$Genotype.number <- as.factor(combined$Genotype.number)
combined$N2.Treatment <- as.factor(combined$N2.Treatment)

#ok now just keep what you want, which is 
colnames(combined)
#Genotype.number,Genotype.name,N2.Treatment,Dry.weight..g.,Count_active.nodule,Count_senesced.nodule,Count_immature.nodule
combined <- combined[,c("Genotype.number","Genotype.name","N2.Treatment","Dry.weight..g.","Count_active.nodule","Count_senesced.nodule","Count_immature.nodule")] #1920
 

combined <- dplyr::filter(
  combined,
  !stringr::str_detect(Genotype.name, "(Sym|sym)") #1830, correct
)



#Intermediary step, which is to collapse all syms into just one per genotype for all "batches"
combined_filt_lent <- dplyr::filter(
  combined,
  stringr::str_detect(Genotype.name, "(Sym|sym)")
)

# Simplify Genotype.name by keeping only "Sym9" or "Sym19" and removing any trailing text
#“If Genotype.name starts with Sym9 or Sym19, remove everything after that and keep only Sym9 or Sym19
#\\1 refers to capture group1,Capture what you want to keep, match everything else, then replace with \\1.
#( ) → defines a capture group
#\\1, \\2, \\3, … → refer to capture groups in the replacement
#Double backslash is needed because R escapes strings first
combined_filt_lent$Genotype.name <- str_replace(combined_filt_lent$Genotype.name, "(sym9|sym19|Sym9|Sym19).*", "\\1")

#Uppercase, remove spaces
combined_filt_lent$Genotype.name <- combined_filt_lent$Genotype.name |>
  stringr::str_to_title() |>
  stringr::str_replace_all("\\s+", "")


View(combined_filt_lent)

combined_control_melt= reshape2::melt(combined_filt_lent,id.vars=c("Genotype.name","N2.Treatment"),na.rm=TRUE)#na.rm removes NA

head(combined_control_melt$variable)

combined_control_melt <- dplyr::filter(
  combined_control_melt,
  !stringr::str_detect(variable, "Genotype.number")
)

combined_control_melt$value <- as.numeric(combined_control_melt$value)


control_peas_average_stats <- plyr::ddply(
  combined_control_melt,
  c("Genotype.name", "N2.Treatment","variable"),
  summarise,
  mean = mean(value, na.rm = TRUE),
  sd   = sd(value, na.rm = TRUE),
  n    = sum(!is.na(value)),
  sem  = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
)

View(control_peas_average_stats)

colnames(control_peas_average_stats)
levels(control_peas_average_stats$N2.Treatment)
control_peas_average_stats$N2.Treatment <- factor(control_peas_average_stats$N2.Treatment,levels=c("3.5 mM ","1 mM"))
levels(control_peas_average_stats$N2.Treatment)
control_peas_average_stats$variable <- as.factor(control_peas_average_stats$variable)
control_peas_average_stats$Genotype.number <- control_peas_average_stats$Genotype.name
control_peas_average_stats$Genotype.number <- as.factor(control_peas_average_stats$Genotype.number)
colnames(control_peas_average_stats)
control_peas_average_stats <- control_peas_average_stats[,c(8,1:7)]

#


#Make Long
combined_melt= reshape2::melt(combined,id.vars=c("Genotype.number","Genotype.name","N2.Treatment"),na.rm=TRUE)#na.rm removes NA values, groups by week and treatment for long data format <-
combined_melt$value <- as.numeric(combined_melt$value)

str(combined_melt)


#PBA prefix highlight as different color per PBA cultivar in the faceted graph

pea_average_stats <- plyr::ddply(
  combined_melt,
  c("Genotype.number", "Genotype.name", "N2.Treatment","variable"),
  summarise,
  mean = mean(value, na.rm = TRUE),
  sd   = sd(value, na.rm = TRUE),
  n    = sum(!is.na(value)),
  sem  = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
)
View(pea_average_stats)
write.csv(pea_average_stats,file="./Meysam/Pea_Avge_Stats.csv",row.names=FALSE)
pea_average_stats <- read.csv("Meysam/Pea_Avge_Stats.csv",header = TRUE)

colnames(pea_average_stats)
levels(pea_average_stats$N2.Treatment)
pea_average_stats$N2.Treatment <- factor(pea_average_stats$N2.Treatment,levels=c("3.5 mM ","1 mM"))
levels(pea_average_stats$N2.Treatment)
pea_average_stats$variable <- as.factor(pea_average_stats$variable)
pea_average_stats$Genotype.number <- as.factor(pea_average_stats$Genotype.number)

#Combine
control_peas_average_stats_combined <- rbind(pea_average_stats,control_peas_average_stats)

write.csv(control_peas_average_stats_combined,file="./Meysam/Pea_Avge_Stats_Combined.csv",row.names=FALSE)
control_peas_average_stats_combined <- read.csv("Meysam/Pea_Avge_Stats_Combined.csv",header = TRUE)



pea_average_stats_sortHigh  <- control_peas_average_stats_combined %>%
  dplyr::mutate(Genotype.number = factor(Genotype.number,
                           levels = control_peas_average_stats_combined %>%
                             filter(N2.Treatment == "3.5 mM ", variable == "Dry.weight..g.") %>%
                             arrange(desc(mean)) %>%
                             pull(Genotype.number)
  )
  )

pea_average_stats_sortLow  <- control_peas_average_stats_combined %>%
  dplyr::mutate(Genotype.number = factor(Genotype.number,
                                         levels = control_peas_average_stats_combined %>%
                                           filter(N2.Treatment == "1 mM", variable == "Dry.weight..g.") %>%
                                           arrange(desc(mean)) %>%
                                           pull(Genotype.number)
  )
  )
              

View(pea_average_stats_sortLow)

#Biomass low N sort first
#P242-249 needs to be separate colors for peas for both barcharts and scatterplots
#Sym9 and Sym19

Combined_lowN_biomass_lowN_sort <- ggplot(subset(pea_average_stats_sortLow, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  #facet_wrap(~ Trt, ncol = 1, scales = "fixed") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by low N, combined") +
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
ggsave(Combined_lowN_biomass_lowN_sort,file="./Meysam/Combined_lowN_biomass_lowN_sort2.pdf",width=27.5,height=12,units="cm")

lowN_biomass_lowN_sort_facet <- ggplot(subset(pea_average_stats_sortLow, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  facet_wrap(~ N2.Treatment, ncol = 1, scales = "free",labeller = "label_both") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by low N, combined") +
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
ggsave(lowN_biomass_lowN_sort_facet,file="./Meysam/Peas_lowN_biomass_lowN_sort_facet2.pdf",width=27.5,height=12*2,units="cm")

print(pea_average_stats_sortHigh$variable)

Combined_lowN_actnod_lowN_sort <- ggplot(subset(pea_average_stats_sortLow, variable %in% "Count_active.nodule"), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  #facet_wrap(~ Trt, ncol = 1, scales = "fixed") +
  labs(y = "Active Nodule Count", x = "Genotype",title="Active Nodule Count sorted by , combined") +
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


Combined_lowN_actnod_lowN_sort


#Scatterplots

str(pea_average_stats)
peas_wide <- control_peas_average_stats_combined %>%
   pivot_wider(
    id_cols = c(Genotype.number, Genotype.name, N2.Treatment), #id_cols → the identifiers you want to keep as rows (Genotype.number, Genotype.name, N2.Treatment).
    names_from = variable, #names_from = variable → the factor column whose levels become new column names.
    values_from = c(mean, sd, n, sem) #values_from = c(mean, sd, n, sem) → the measurement columns that get spread into those new variable-specific columns.
  )

peas_wide$N2.Treatment <- factor(peas_wide$N2.Treatment,levels=c("1 mM", "3.5 mM "))


#Scatterplot
# keep only the biomass variable
biomass <- control_peas_average_stats_combined %>%
  filter(variable == "Dry.weight..g.") %>%
  dplyr::select(Genotype.number, Genotype.name, N2.Treatment, mean) %>%
  pivot_wider(
    names_from = N2.Treatment,
    values_from = mean
  ) %>%
  rename(Low = `1 mM`, High = `3.5 mM `)

Scatter_pea <- ggplot(biomass, aes(x = Low, y = High)) +
  geom_point(aes(shape = Genotype.number %in% c("P109","P181"),color = Genotype.number %in% c("P109", "P181"))) +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
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
    legend.position = "none"   # <-- hides the legend
  ) +
  labs(title = "Peas diversity screen: Shoot dry weight (g)")


Scatter_pea

ggsave(Scatter_pea,file="./Meysam/Scatterplot_Peas.pdf",units="cm",height=25,width=30)

#Biomass High N vs Nodule Count
str(peas_wide)
colnames(peas_wide)
str(pea_average_stats)
HighN_BiomassvsActNod <- ggplot(subset(peas_wide,N2.Treatment %in% "3.5 mM "), aes(x = mean_Dry.weight..g., y = mean_Count_active.nodule)) +
  geom_point(aes(shape = Genotype.number %in% c("P109","P181"),
    color = Genotype.number %in% c("P109", "P181"))) +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
 # ylim(-5,NA)+
  #xlim(-5,NA)+
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
   legend.position = "none"   # <-- hides the legend
  ) +
  labs(title = "Peas diversity screen: Shoot dry weight (g) vs. active nodules for High N")
HighN_BiomassvsActNod


#Let's do something combined
Combined_vsactNods <- ggplot(peas_wide,aes(x = mean_Dry.weight..g., y = mean_Count_active.nodule)) +
  geom_point(aes(shape = Genotype.number %in% c("Sym9","Sym119"),
                 color = Genotype.number %in% c("Sym9", "Sym119"))) +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "glm", se = TRUE,color="orange",fill="orange",linewidth=0.5) +
  facet_wrap(~N2.Treatment) +
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
    legend.position = "none"   # <-- hides the legend
  ) +
  labs(title = "Peas diversity screen: Shoot dry weight (g) vs. active nodules")

Combined_vsactNods

#Filter out for above to keep only the controls
Combined_vsactNods <- ggplot(
  peas_wide,
  aes(x = mean_Dry.weight..g., y = mean_Count_active.nodule)
) +
  geom_point(
    aes(
      shape = Genotype.number %in% c("P242","P243","P244","P245","P246","P247","P248","P249","Sym9","Sym119"),
      colour = Genotype.number %in% c("P242","P243","P244","P245","P246","P247","P248","P249","Sym9","Sym119"),
      alpha=0.3
    )
  ) +
  geom_text(
    data = dplyr::filter(
      peas_wide,
      Genotype.number %in% c("P242","P243","P244","P245","P246","P247","P248","P249","Sym9","Sym119")
    ),
    aes(label = Genotype.number),
    size = 3,
    vjust=-0.5,
    check_overlap = FALSE
  ) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm", se = TRUE, color = "orange", fill = "orange", linewidth = 0.5) +
  facet_wrap(~N2.Treatment) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x  = element_text(size = 7.5, face = "bold"),
    axis.text.y  = element_text(size = 7.5, face = "bold"),
    strip.text   = element_text(size = 11, face = "bold", hjust = 0.01),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "gray80"),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(title = "Lentil diversity screen: Shoot dry weight (g) vs. active nodules")

Combined_vsactNods
ggsave(Combined_vsactNods,file="./Meysam/Peas_DW_vs_ActNods_notext2.pdf",units="cm",height=20,width=30)

Combined_vssenesNods <- ggplot(peas_wide,aes(x = mean_Dry.weight..g., y = mean_Count_senesced.nodule)) +
  geom_point(aes(shape = Genotype.number %in% c("P109","P181"),
                 color = Genotype.number %in% c("P109", "P181"))) +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "glm", se = TRUE,color="orange",fill="orange",linewidth=0.5) +
  facet_wrap(~N2.Treatment) +
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
    legend.position = "none"   # <-- hides the legend
  ) +
  labs(title = "Peas diversity screen: Shoot dry weight (g) vs. senescing nodules")

Combined_vssenesNods

CombinedBiomasNods <- grid.arrange(Combined_vsactNods,Combined_vssenesNods)

ggsave(CombinedBiomasNods,file="./Meysam/CombinedBiomasNods.pdf",units="cm",height=50,width=50)


HighN_BiomassvsActNod



ggplot(peas_wide, aes(x = mean_Dry.weight..g., y = mean_Count_active.nodule)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~N2.Treatment) +
  theme_bw()

HighN_BiomassvsActNod +
  annotate("text", x = 200, y = 3, label = paste0("r = ", round(cor_val, 2)), size = 5)

HighN_BiomassvsActNod +
  #geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1)+
  geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 1,alpha=0.5)

HighN_BiomassvsSenesNod <- ggplot(subset(peas_wide,N2.Treatment %in% "3.5 mM "), aes(x = mean_Dry.weight..g., y = mean_Count_senesced.nodule)) +
  geom_point(aes(color = Genotype.number %in% c("P109", "P181"))) +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("TRUE" = "lightblue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
  # ylim(-5,NA)+
  #xlim(-5,NA)+
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
    legend.position = "none"   # <-- hides the legend
  ) +
  labs(title = "Peas diversity screen: Shoot dry weight (g) vs. senesced nodules for High N")

HighN_BiomassvsSenesNod

CombinedHighN <- grid.arrange(HighN_BiomassvsActNod,HighN_BiomassvsSenesNod,ncol=1)


warnings()
ggplot(
  subset(peas_wide, N2.Treatment %in% "3.5 mM " & 
           !is.na(mean_Dry.weight..g.) & 
           !is.na(mean_Count_active.nodule)),
  aes(x = mean_Count_active.nodule, y = mean_Dry.weight..g.)
) +
  geom_point() +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  theme_bw() +
  labs(title = "Peas diversity screen: Shoot dry weight (g) vs. active nodules for High N",
       x = "Active nodules (mean count)",
       y = "Shoot dry weight (g)")











Combined_highN_biomass_highN_sort <- ggplot(subset(pea_average_stats_sortHigh, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  #facet_wrap(~ Trt, ncol = 1, scales = "fixed") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by High N, faceted") +
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

Combined_highN_biomass_highN_sort
ggsave(Combined_highN_biomass_highN_sort,file="./Meysam/Combined_highN_biomass_highN_sort.pdf",width=27.5,height=12,units="cm")

HighN_biomass_highN_sort_facet <- ggplot(subset(pea_average_stats_sortHigh, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  facet_wrap(~ N2.Treatment, ncol = 1, scales = "free",labeller = "label_both") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by High N, faceted") +
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

HighN_biomass_highN_sort_facet
ggsave(HighN_biomass_highN_sort_facet,file="./Meysam/HighN_biomass_highN_sort_facet.pdf",width=27.5,height=12*2,units="cm")






























Combined_lowN_biomass_facet <- ggplot(subset(pea_average_stats_sortLow, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean,fill=N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  facet_grid(~ N2.Treatment, scales = "fixed") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by low N,faceted") +
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

Combined_lowN_biomass_facet
ggsave(Combined_lowN_biomass_facet,file="test.pdf",units="cm",height=25,width=70)






Combined_HighN

ggsave(Combined,file="HighNSorted_Combined_DW.pdf",width=50,height=25,units="cm")


