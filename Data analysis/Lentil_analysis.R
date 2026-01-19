#Lentil_analysis with more response variables
install.packages("ggtext")

#libraries
#load required packages
#PCA Analysis and plotting
#BiocManager::install("pcaExplorer")
#library("pcaExplorer")
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



#rawdata
lentil_biomass_HighN <- read.csv("Meysam/GRDC milestone 105_Lentil_Nov 2025_HighN.csv",header=TRUE)
lentil_biomass_lowN <- read.csv("Meysam/GRDC milestone 105_Lentil_Nov 2025_LowN.csv",header=TRUE)
colnames(lentil_biomass_HighN)
colnames(lentil_biomass_lowN)

combined <- rbind(lentil_biomass_HighN,lentil_biomass_lowN)


#processdata
head(combined$Genotype.number) #add a L
combined$Genotype.number <- paste0("L", combined$Genotype.number)
head(combined$Genotype.number)
combined$Genotype.number <- as.factor(combined$Genotype.number)
combined$N2.Treatment <- as.factor(combined$N2.Treatment)
colnames(combined)
#Genotype.number,Genotype.name,N2.Treatment,Dry.weight..g.,Count_active.nodule,Count_senesced.nodule,Count_immature.nodule
combined <- combined[,c("Genotype.number","Genotype.name","N2.Treatment","Dry.weight..g.","Count_active.nodule","Count_senesced.nodule","Count_immature.nodule")]



#now get rid of peas # keep syms
combined_filt <-combined %>% 
  dplyr::filter(!str_detect(Genotype.name, "Sym"))

#Get the Sym values
combined_control <- (dplyr::filter(combined,str_detect(Genotype.name, "Sym")))
combined_control_summarizeSyms <- summarise(combined_control,)

combined_control$Genotype.name <- str_replace(combined_control$Genotype.name, "(Sym9|Sym19).*", "\\1")

combined_control_melt= reshape2::melt(combined_control,id.vars=c("Genotype.name","N2.Treatment"),na.rm=TRUE)#na.rm removes NA

combined_control_melt$value <- as.numeric(combined_control_melt$value)

control_lentils_average_stats <- plyr::ddply(
  combined_control_melt,
  c("Genotype.name", "N2.Treatment","variable"),
  summarise,
  mean = mean(value, na.rm = TRUE),
  sd   = sd(value, na.rm = TRUE),
  n    = sum(!is.na(value)),
  sem  = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
)

control_lentils_average_stats$Genotype.number <- control_lentils_average_stats$Genotype.name
control_lentils_average_stats <- control_lentils_average_stats[,c(8,1:7)]

colnames(control_lentils_average_stats)
levels(control_lentils_average_stats$N2.Treatment)
control_lentils_average_stats$N2.Treatment <- factor(control_lentils_average_stats$N2.Treatment,levels=c("3.5 mM ","1 mM"))
levels(control_lentils_average_stats$N2.Treatment)
control_lentils_average_stats$variable <- as.factor(control_lentils_average_stats$variable)
control_lentils_average_stats$Genotype.number <- as.factor(control_lentils_average_stats$Genotype.number)
View(control_lentils_average_stats)



#now some stats
combined_melt= reshape2::melt(combined_filt,id.vars=c("Genotype.number","Genotype.name","N2.Treatment"),na.rm=TRUE)#na.rm removes NA values, groups by week and treatment for long data format
combined_melt$value <- as.numeric(combined_melt$value)

lentils_average_stats <- plyr::ddply(
  combined_melt,
  c("Genotype.number", "Genotype.name", "N2.Treatment","variable"),
  summarise,
  mean = mean(value, na.rm = TRUE),
  sd   = sd(value, na.rm = TRUE),
  n    = sum(!is.na(value)),
  sem  = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
)
View(lentils_average_stats)
write.csv(lentils_average_stats,file="./Meysam/Lentil_Avge_Stats.csv",row.names=FALSE)
lentils_average_stats <- read.csv("Meysam/Pea_Avge_Stats.csv",header = TRUE)

#reproduce bar graphs (maybelater)
colnames(lentils_average_stats)
levels(lentils_average_stats$N2.Treatment)
lentils_average_stats$N2.Treatment <- factor(lentils_average_stats$N2.Treatment,levels=c("3.5 mM ","1 mM"))
levels(lentils_average_stats$N2.Treatment)
lentils_average_stats$variable <- as.factor(lentils_average_stats$variable)
lentils_average_stats$Genotype.number <- as.factor(lentils_average_stats$Genotype.number)

#Ok now we add back Syms at meysams'request

Combinedlentils_average_stats <- rbind(lentils_average_stats,control_lentils_average_stats) #2456
Combinedlentils_average_stats <- (dplyr::filter(Combinedlentils_average_stats,!stringr::str_detect(variable, "Genotype.number")))
write.csv(Combinedlentils_average_stats,file="./Meysam/Lentil_Avge_Stats_Combined.csv",row.names=FALSE)
Combinedlentils_average_stats <- read.csv("Meysam/Lentil_Avge_Stats_Combined.csv",header = TRUE)



lentil_average_stats_sortHigh  <- Combinedlentils_average_stats %>%
  dplyr::mutate(Genotype.number = factor(Genotype.number,
                                         levels = lentils_average_stats %>%
                                           filter(N2.Treatment == "3.5 mM ", variable == "Dry.weight..g.") %>%
                                           arrange(desc(mean)) %>%
                                           pull(Genotype.number)
  )
  )

lentil_average_stats_sortLow  <- Combinedlentils_average_stats %>%
  dplyr::mutate(Genotype.number = factor(Genotype.number,
                                         levels = Combinedlentils_average_stats %>%
                                           filter(N2.Treatment == "1 mM", variable == "Dry.weight..g.") %>%
                                           arrange(desc(mean)) %>%
                                           pull(Genotype.number)
  )
  )


# Get unique genotype labels
genos <- unique(lentil_average_stats_sortLow$Genotype.number)

# Assign colors: L1 = red, others = black
geno_colors <- setNames(
  ifelse(genos == "L228","Sym9","Sym19", "red", "blue","purple","Transparent"),
  genos
)

geno_colors <- setNames(
  case_when(
    genos == "L228" ~ "red",
    genos == "Sym9" ~ "blue",
    genos == "Sym19" ~ "purple",
    TRUE ~ "black"
  ),
  genos
)


print(geno_colors)
#L228 to color separately in faceted graph


# Add this to your plot:
scale_y_discrete(
  labels = function(x) {
    paste0("<span style='color:", geno_colors[x], ";'>", x, "</span>")
  }
) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 10, face = "bold")
  )

#Now start plotting


Combined_lowN_biomass_lowN_sort <- ggplot(subset(lentil_average_stats_sortLow, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  #facet_wrap(~ Trt, ncol = 1, scales = "fixed") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by low N, combined (Lentils)") +
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

lowN_biomass_lowN_sort_facet <- ggplot(subset(lentil_average_stats_sortLow, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  facet_wrap(~ N2.Treatment, ncol = 1, scales = "free",labeller = "label_both") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by low N, combined (Lentils)") +
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

ggsave(lowN_biomass_lowN_sort_facet,file="./Meysam/Lents_Combined_lowN_biomass_lowN_sort_included.pdf",width=27.5,height=12*2,units="cm")



Combined_highN_actnod_highN_sort <- ggplot(subset(lentil_average_stats_sortHigh, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  #facet_wrap(~ Trt, ncol = 1, scales = "fixed") +
  labs(y = "Active Nodule Count", x = "Genotype",title="Active Nodule Count sorted by hi , combined (Lentils") +
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

Combined_highN_actnod_highN_sort

highN_biomass_highN_sort_facet <- ggplot(subset(lentil_average_stats_sortLow, variable %in% "Dry.weight..g."), aes(x = Genotype.number, y = mean, fill = N2.Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.2,colour = "black") +
  facet_wrap(~ N2.Treatment, ncol = 1, scales = "free",labeller = "label_both") +
  labs(y = "Biomass (g DW)", x = "Genotype",title="Biomass sorted by high N, combined (Lentils") +
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




#scatterplots
str(Combinedlentils_average_stats)
lentils_wide <- Combinedlentils_average_stats %>%
  pivot_wider(
    id_cols = c(Genotype.number, Genotype.name, N2.Treatment), #id_cols → the identifiers you want to keep as rows (Genotype.number, Genotype.name, N2.Treatment).
    names_from = variable, #names_from = variable → the factor column whose levels become new column names.
    values_from = c(mean, sd, n, sem) #values_from = c(mean, sd, n, sem) → the measurement columns that get spread into those new variable-specific columns.
  )

biomass <- Combinedlentils_average_stats %>%
  filter(variable == "Dry.weight..g.") %>%
  dplyr::select(Genotype.number, Genotype.name, N2.Treatment, mean) %>%
  pivot_wider(
    names_from = N2.Treatment,
    values_from = mean
  ) %>%
  rename(Low = `1 mM`, High = `3.5 mM `)

Scatter_lent <- ggplot(biomass, aes(x = Low, y = High)) +
  geom_point(aes(shape = factor(Genotype.number %in% c("L300","Sym9","Sym119")),color = Genotype.number %in% c("L300","Sym9","Sym119"))) +
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
    legend.position = "none" # <-- hides the legend  # <-- hides the legend
  ) +
  labs(title = "Lentil diversity screen: Shoot dry weight (g)")

Scatter_lent
ggsave(Scatter_lent,file="./Meysam/Scatterplot_Lent.pdf",units="cm",height=25,width=30)

Combined_vsactNods <- ggplot(lentils_wide,aes(x = mean_Dry.weight..g., y = mean_Count_active.nodule)) +
  geom_point(aes(shape = Genotype.number %in% c("L300","Sym9","Sym119"),
                 color = Genotype.number %in% c("L300","Sym9","Sym119"))) +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm", se = TRUE,color="orange",fill="orange",linewidth=0.5) +
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
  labs(title = "Lentil diversity screen: Shoot dry weight (g) vs. active nodules")

Combined_vsactNods

#Filter out for above to keep only the controls
Combined_vsactNods <- ggplot(
  lentils_wide,
  aes(x = mean_Dry.weight..g., y = mean_Count_active.nodule)
) +
  geom_point(
    aes(
      shape = Genotype.number %in% c("L228","Sym9","Sym119"),
      color = Genotype.number %in% c("L228","Sym9","Sym119")
    )
  ) +
  geom_text(
    data = dplyr::filter(
      lentils_wide,
      Genotype.number %in% c("L228","Sym9","Sym119")
    ),
    aes(label = Genotype.number),
    size = 3,
    vjust = -0.5
  ) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm", se = TRUE, color = "orange", fill = "orange", linewidth = 0.5) +
  facet_wrap(~N2.Treatment, scales="fixed") +
  xlim(0,2.2)+
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

ggsave(Combined_vsactNods,file="./Meysam/Lentils_DW_vs_ActNods_notext3.pdf",units="cm",height=20,width=30)


Combined_vssenesNods <- ggplot(lentils_wide,aes(x = mean_Dry.weight..g., y = mean_Count_senesced.nodule)) +
  geom_point(aes(shape = Genotype.number %in% c("L300"),
                 color = Genotype.number %in% c("L300"))) +
  geom_text(aes(label = Genotype.number), check_overlap = FALSE, size = 3, vjust = -0.5) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black")) +
  theme_bw(base_size = 12) +
  geom_smooth(method = "lm", se = TRUE,color="orange",fill="orange",linewidth=0.5) +
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
  labs(title = "Lentil diversity screen: Shoot dry weight (g) vs. senescing nodules")

Combined_vssenesNods
colnames(lentils_wide)
levels(lentils_wide$N2.Treatment)
model <- lm(mean_Count_senesced.nodule ~ mean_Dry.weight..g., data = subset(lentils_wide, N2.Treatment %in% "3.5 mM "))
summary(model)$r.squared

CombinedBiomasNods <- grid.arrange(Combined_vsactNods,Combined_vssenesNods)

ggsave(CombinedBiomasNods,file="./Meysam/CombinedBiomasNodsLentil.pdf",units="cm",height=50,width=50)
