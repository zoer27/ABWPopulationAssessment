#Manuscript Figures
#Zoe Rand
#Note that current file paths should be adjusted to where model results are stored
library(tidyverse)
library(posterior) #manipulating chains
library(patchwork) #plotting
library(scales) #plotting
library(ggforce) #plotting
library(paletteer) # for colors


# Read in Data ------------------------------------------------------------
#catch data
catch_dat<-read_csv("data/ABW_catch_series_2024.csv")

#abundance data
abund_dat<-read_csv("data/Abund_Est.csv")
#remove 6 from JARPA
abund_dat<-abund_dat[abund_dat$N > 6,]
JSV_dat<-read_csv("data/JSV.csv")
photo_id<-read_csv("data/photoID_abund.csv")


#turn CVS into confidence intervals
abund_dat<-abund_dat[abund_dat$N > 6,]
abund_dat$C<-exp(1.96 * sqrt(log(1+abund_dat$CV^2)))
abund_dat$Lower<-abund_dat$N/abund_dat$C
abund_dat$Upper<-abund_dat$N*abund_dat$C

photo_id$C<-exp(1.96 * sqrt(log(1+photo_id$CV^2)))
photo_id$Lower<-photo_id$N/photo_id$C
photo_id$Upper<-photo_id$N*photo_id$C

#other constants
#number of draws
basemod_ndraws<-4000
abund_years_J<-rep(c(abund_dat$Year, photo_id$Year[1], JSV_dat$Year), basemod_ndraws)
abund_years_B<-rep(c(abund_dat$Year, photo_id$Year[1]), basemod_ndraws)
abund_years<-tibble(Year = c(abund_years_J, abund_years_B), model = c(rep("JSV", length(abund_years_J)), rep("base", length(abund_years_B))))
source_lab_J<-rep(c(rep("SOWER", 3), rep("JARPA", 9), "Photo ID", rep("JSV", 21)), basemod_ndraws)
source_lab_B<-rep(c(rep("SOWER", 3), rep("JARPA", 9), "Photo ID"), basemod_ndraws)
source_labs<-tibble(source = c(source_lab_J, source_lab_B), model = c(rep("JSV", length(source_lab_J)), rep("base", length(source_lab_B))))

# Read in Model Fits ------------------------------------------------------

#Base model
#model fits were saved without burn in
B_fit1<-read.csv("Code/Results/Fit_BaseModel_uniformR_112.csv")
B_fit2<-read.csv("Code/Results/Fit_BaseModel_calldensityR_112.csv")
B_fit3<-read.csv("Code/Results/Fit_BaseModel_uniformR_239.csv")
B_fit4<-read.csv("Code/Results/Fit_JSVModel_uniformR_112.csv")
B_fit5<-read.csv("Code/Results/Fit_BaseModel_nophoto.csv")
B_fit6<-read.csv("Code/Results/Fit_BaseModel_noJARPA.csv")
B_fit7<-read.csv("Code/Results/Fit_BaseModel_noSOWER.csv")
B_fit8<-read.csv("Code/Results/Fit_BaseModel_justSOWER.csv")
B_fit9<-read.csv("Code/Results/Fit_BaseModel_fixed_s.csv")


#convert to draws mat
draws_mat_uR_112<-as_draws_matrix(B_fit1)
draws_mat_cR_112<-as_draws_matrix(B_fit2)
draws_mat_uR_239<-as_draws_matrix(B_fit3)
draws_mat_JSV<-as_draws_matrix(B_fit4)
draws_mat_nophoto<-as_draws_matrix(B_fit5)
draws_mat_noJARPA<-as_draws_matrix(B_fit6)
draws_mat_noSOWER<-as_draws_matrix(B_fit7)
draws_mat_justSOWER<-as_draws_matrix(B_fit8)
draws_mat_fixeds<-as_draws_matrix(B_fit9)



#convert to list
draws_list<-list(draws_mat_uR_112, 
                 draws_mat_JSV, 
                 draws_mat_cR_112, 
                 draws_mat_uR_239, 
                 draws_mat_fixeds,
                 draws_mat_nophoto, 
                 draws_mat_noJARPA, 
                 draws_mat_noSOWER, 
                 draws_mat_justSOWER
)


#posterior predictive for base model
post_pred_base<-read.csv("code/results/PostPred_BaseModel_uniformR_112.csv")
post_pred_base_mat<-as_draws_matrix(post_pred_base)

# Helper Functions ------------------------------------------------------
#subset draws by parameter
subset_draws_par<-function(draws_mat, par, mod){
  post<-subset_draws(draws_mat, par)
  posttab<-post %>% as_draws_matrix() %>% as_tibble() %>%
    mutate(model = mod)
  return(posttab)
}

#get population trajectory and quantiles
Npopkeys<-paste0("Npop", 1:121)
subset_npop<-function(drawsmat, mod){
  post<-subset_draws(drawsmat, Npopkeys)
  pop_quants<-apply(post, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) %>% as_tibble() %>%
    add_column(est = c("2.5%", "50%", "97.5%"), .before = 1) %>%
    pivot_longer(-(est), names_to = "par", values_to = "N") %>%
    add_column(year = rep(1904:2024, 3), .before = 1) %>% 
    pivot_wider(names_from = "est", values_from = "N") %>% 
    add_column(model = mod)
  return(pop_quants)
}

#add additional CVs to confidence intervals
addCV_func<-function(addCV_branch, addCV_ham, mod){
  AddCVs<-c(rep(addCV_branch, 3), rep(addCV_ham, 9))
  AddCVsC<-exp(1.96 * sqrt(log(1+abund_dat$CV^2 + AddCVs^2)))
  out<-abund_dat %>% add_column(model = mod) %>% 
    mutate(C = AddCVsC, AddLower = N/C, AddUpper = N*C)
  return(out)
}

#get out posterior predictive
post_out<-function(post_pred_df, md){
  if(md == "JSV"){
    out<-post_pred_df %>% dplyr::select(-c("iteration", "chain")) %>% 
      dplyr::select("N.1.", "N.2.", "N.3.", "N.4.", "N.5.", "N.6.", "N.7.", "N.8.", "N.9.", "N.10.", "N.11.", "N.12.", "photo_n.1.", 
                    "sightings_early.1.", "sightings_early.2.", "sightings_early.3.", "sightings_early.4.", "sightings_early.5.", "sightings_early.6.", "sightings_early.7.", 
                    "sightings_early.8.", "sightings_early.9.", "sightings_early.10.", "sightings_early.11.", "sightings_late.1.", "sightings_late.2.", "sightings_late.3.","sightings_late.4.", 
                    "sightings_late.5.", "sightings_late.6.", "sightings_late.7.", "sightings_late.8.", "sightings_late.9.", "sightings_late.10.") %>%
      pivot_longer(everything(), names_to = "est", values_to = "value") %>%                                                                        
      add_column(year = abund_years$Year[abund_years$model==md], .before = 1) %>%
      add_column(source = source_labs$source[source_labs$model == md], .before = 1)
  }else{
    out<-post_pred_df %>% dplyr::select(-c("iteration", "chain")) %>% 
      dplyr::select("N.1.", "N.2.", "N.3.", "N.4.", "N.5.", "N.6.", "N.7.", "N.8.", "N.9.", "N.10.", "N.11.", "N.12.", "photo_n.1.") %>%
      pivot_longer(everything(), names_to = "est", values_to = "value") %>%                                                                        
      add_column(year = abund_years$Year[abund_years$model==md], .before = 1) %>%
      add_column(source = source_labs$source[source_labs$model == md], .before = 1)
  }
  
  return(out)
}

# Colors ------------------------------------------------------------------
#color palette for sensitivities
pal<-c("#000000", '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#BF5700', '#cab2d6','#6a3d9a')


bckgrd<-"transparent"

model_names_txt<-c(
  "Base", 
  "JSV",
  "Call density r",
  "Adjusted \u03B8  ",
  "Fixed survival",
  "No photo ID",
  "No JARPA",
  "No SOWER", 
  "Just SOWER" 
)

model_names<-factor(model_names_txt, levels = c("Base", 
                                                 "JSV",
                                                 "Call density r",
                                                 "Adjusted \u03B8  ",
                                                 "Fixed survival",
                                                 "No photo ID",
                                                 "No JARPA",
                                                 "No SOWER", 
                                                 "Just SOWER"))
levels(model_names)
#name color palette
names(pal)<-model_names

#color palette for data sources
data_pal <- c("#133C55", "#59A5D8", "#91E5F6")

# Catch data plot--------------------------------------------------------------
head(catch_dat)
catch_plot<-ggplot(catch_dat[-1,]) + geom_col(aes(x = Year, y = ABW_catch), fill = "deepskyblue4", color = "white") + 
  theme_classic() + 
  scale_x_continuous(expand =c(0,0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Whaling season", y = "Catch (#s)")
catch_plot
ggsave("Figures/Manuscript/CatchPlot.png",catch_plot,  width = 4, height = 3, units = "in", dpi = 600)


# Plotting Functions ------------------------------------------------------

plot_R<-function(r_post, prior, r_prior, r_prior_dens, base, r_post_base){
  if(base){
    rpostplot<-ggplot() + 
      geom_density(data = r_post, aes(x = r, y = after_stat(density), group = model),fill = "gray80", 
                   key_glyph = draw_key_path, adjust = 1.2)
  }
    if(prior){
  rpostplot<- rpostplot + geom_density(aes(x = r_prior, y = r_prior_dens), linetype = "dotted") 
    }
  if(!base){
    rpostplot<-ggplot() + 
      geom_density(data = r_post, aes(x = r, y = after_stat(density), group = model, color = model, linewidth = model),fill = NA, 
                   key_glyph = draw_key_path, adjust = 1.2)
    rpostplot<- rpostplot + geom_density(data = r_post_base, aes(x = r, y = after_stat(density), group = model, color = model, linewidth = model), fill = NA, key_glyph = draw_key_path)
  }
  rpostplot<- rpostplot + 
    theme_classic() + 
    labs(x = expression(paste("Intrinsic growth (", r[max], ")")), y = "density", color = "Model", linewidth = "Model") + 
    scale_color_manual(values = pal) + 
    #scale_fill_manual(values = pal) +
    scale_linewidth_manual(values = c(1, 1, 1, 1, 1, 1, 1, 1, 2)) +
    guides(color = guide_legend(reverse = TRUE), linewidth = guide_legend(reverse = TRUE)) +
    scale_x_continuous(expand = c(0,0), limits  = c(0,0.12), breaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.11)) + 
    scale_y_continuous(expand = c(0,0.1)) + theme(panel.background = element_rect(fill = bckgrd), 
                                              plot.background = element_rect(fill = bckgrd),
                                              panel.border = element_blank(),
                                              #text = element_text(size = 18, family = "Arial"),
                                              plot.tag.position = c(0, 1), 
                                              #plot.tag = element_text(vjust = 1.3, hjust = -1), 
                                              axis.text.y = element_blank(), 
                                              axis.ticks.y = element_blank(), 
                                              axis.title.x = element_text(vjust=1))
                                              #legend.position = "none")
  print(rpostplot)
  return(rpostplot)
}

plot_K<-function(K_post, base, K_post_base){
  if(base){
    kpostplot<-ggplot() + 
      geom_density(data = K_post, aes(x = K, y = after_stat(density), group = model),  fill = "gray80",
                   key_glyph = draw_key_path, 
                   adjust = 1.2) 
  }
  if(!base){
    kpostplot<-ggplot() + 
      geom_density(data = K_post, aes(x = K, y = after_stat(density), group = model, 
                                      color = model, linewidth = model), key_glyph = draw_key_path, 
                   adjust = 1.2) 
    kpostplot<- kpostplot + geom_density(data = K_post_base, aes(x = K, y = after_stat(density), group = model, 
                                                                color = model, linewidth = model), key_glyph = draw_key_path)
  }
    kpostplot<-kpostplot + theme_classic() + labs(x = "Carrying capacity (K, # thousands)", y = "density", color = "Model", linewidth = "Model") + 
    scale_x_continuous(expand = c(0.01,0.01), breaks = seq(150000, 400000, by = 20000), 
                       label = seq(150, 400, by = 20)) + 
    #coord_cartesian(xlim = c(140000, 400000)) +
    scale_y_continuous(expand = c(0.002,0)) +
    scale_color_manual(values = pal) +
    scale_linewidth_manual(values = c(1, 1, 1, 1, 1, 1, 1, 1, 2)) +
    guides(color = guide_legend(reverse = TRUE), linewidth = guide_legend(reverse = TRUE)) +
    theme(
      #text = element_text(size = 18, family = "Arial"),
      axis.title.x = element_text(vjust=1),
      plot.tag.position = c(0, 1), 
      #plot.tag = element_text(vjust = 1.3, hjust = -1), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      #axis.title.y = element_blank(),
      legend.position = "none")
  print(kpostplot)
  return(kpostplot)
}

plot_Nmin<-function(Nmin_post, prior, Nmin_prior, Nmin_prior_dens, base, Nmin_post_base){
  if(base){
    nminpostplot<-ggplot() + 
      geom_density(data = Nmin_post, aes(x = Nmin, y = after_stat(density), group = model), 
                                         fill = "gray80",
                   key_glyph = draw_key_path, adjust = 1.2)
  }
  if(prior){
    nminpostplot<- nminpostplot + geom_density(aes(x = Nmin_prior, y = Nmin_prior_dens), linetype = "dotted")
  }
  if(!base){
    nminpostplot<-ggplot() + 
      geom_density(data = Nmin_post, aes(x = Nmin, y = after_stat(density), group = model, 
                                         color = model, linewidth = model), key_glyph = draw_key_path, adjust = 1.2)
    nminpostplot<- nminpostplot + geom_density(data = Nmin_post_base, aes(x = Nmin, y = after_stat(density), group = model, 
                                                                        color = model, linewidth = model), key_glyph = draw_key_path)
  }
  nminpostplot <- nminpostplot + theme_classic() + 
    labs(x = "Minimum population size (#s)", y = "density", color = "Model", linewidth = "Model") + 
    scale_x_continuous(expand = c(0.03,0.03), limits = c(50, 2000))+ 
    scale_y_continuous(expand = c(0.0001,0)) +
    scale_color_manual(values = pal) +
    scale_linewidth_manual(values = c(1, 1, 1, 1, 1, 1, 1, 1, 2)) +
    guides(color = guide_legend(reverse = TRUE), linewidth = guide_legend(reverse = TRUE)) +
    theme(
      #text = element_text(size = 18, family = "Arial"),
      axis.title.x = element_text(vjust=1),
      plot.tag.position = c(0, 1), 
      #plot.tag = element_text(vjust = 1.3, hjust = -1), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      #axis.title.y = element_blank(),
      legend.position = "none")
  print(nminpostplot)
  return(nminpostplot)
}


plot_addCV<-function(add_CVpost, prior, addCV_prior, addCV_prior_dens, source){
  if(source == "JARPA"){
    colnames(add_CVpost)<-gsub("addCV_Ham", "addCV", colnames(add_CVpost))
    title <- "Additional CV (JARPA)"
    tag <- "(f)"
  }
  if(source == "SOWER"){
    colnames(add_CVpost)<-gsub("addCV_Branch", "addCV", colnames(add_CVpost))
    title <- "Additional CV (SOWER)"
    tag <- "(e)"
  }
  addCVplot<-ggplot() +
    geom_density(data = add_CVpost, aes(x = addCV, y = after_stat(density), group = as.factor(model)), 
                 alpha = 1, key_glyph = draw_key_path, adjust = 1.2, fill = "gray80")
  if(prior){
      addCVplot <- addCVplot + geom_line(aes(x = addCV_prior, y = addCV_prior_dens, group = 1), linetype = "dotted")
    }
  addCVplot <- addCVplot + theme_classic() + labs(x = title, y = "density", color = "Model", linewidth = "Model") +
    scale_color_manual(values = pal) +
    scale_linewidth_manual(values = c(1, 1, 1, 1, 1, 1, 1)) +
    guides(color = guide_legend(reverse = TRUE), linewidth = guide_legend(reverse = TRUE)) +
    scale_x_continuous(expand = c(0.01,0.01), limits = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
    scale_y_continuous(expand = c(0.0001,0)) +
    theme(panel.background = element_rect(fill = bckgrd),
          plot.background = element_rect(fill = bckgrd),
          panel.border = element_blank(),
          #text = element_text(size = 18, family = "Arial"),
          plot.tag.position = c(0, 1),
          #plot.tag = element_text(vjust = 1.3, hjust = -1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(vjust=1),
          #axis.title.y = element_blank(),
          legend.position = "none")
  print(addCVplot)
  return(addCVplot)
}

plot_depletion<-function(Nmin_div_k_post){
  depl_plot<-ggplot() + 
    geom_density(data = Nmin_div_k_post, aes(x = Nmin_div_K, y = after_stat(density),
                                       group = as.factor((model))),alpha = 1, 
                 key_glyph = draw_key_path, adjust = 1.2, fill = "gray80") + 
    theme_classic() + labs(x = "Minimum population/K", y = "density", color = "Model") + 
    scale_color_manual(values = pal) +
    scale_linewidth_manual(values = c(1, 1, 1, 1, 1, 1, 1, 2))+
    guides(color = guide_legend(reverse = TRUE), linewidth = guide_legend(reverse = TRUE)) +
    scale_x_continuous(expand = c(0.00002,0.000013), breaks = seq(0, 0.007, by = 0.001), labels =c(0, seq(0.001, 0.007, by = 0.001))) + 
    scale_y_continuous(expand = c(0.003,0)) +
    theme(#axis.text.x = element_text(angle = 90),
      #text = element_text(size = 18, family = "Arial"),
      axis.title.x = element_text(vjust=1),
      plot.tag.position = c(0, 1), 
      #plot.tag = element_text(vjust = 1.3, hjust = -1), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      legend.position = "none")
  print(depl_plot)
  return(depl_plot)
}

plot_N2024<-function(N2024_post, base, N2024_post_base){
  if(base){
    N2024plot<-ggplot() + 
      geom_density(data = N2024_post, aes(x = N2024, y = after_stat(density), 
                                          group = model), fill = "gray80", 
                   key_glyph = draw_key_path, adjust = 1.2)
  }
  if(!base){
    N2024plot<-ggplot() + 
      geom_density(data = N2024_post, aes(x = N2024, y = after_stat(density), 
                                          group = model, color = model,
                                          linewidth = model),alpha = 1, 
                   key_glyph = draw_key_path, adjust = 1.2)
    N2024plot <- N2024plot + geom_density(data = N2024_post_base, aes(x = N2024, y = after_stat(density), 
                                                                      group = model, color = model,
                                                                      linewidth = model), alpha = 1, key_glyph = draw_key_path)
  }
    N2024plot<- N2024plot + theme_classic() + 
    labs(x = "Projected population in 2024 (#s)", y = "density", color = "Model", linewidth = "Model") + 
    scale_x_continuous(limits = c(1000, 14000), expand = c(0.02,0.02), breaks = seq(1000, 14000, by = 2000)) + 
    scale_color_manual(values = pal) +
    scale_linewidth_manual(values = c(1, 1, 1, 1, 1, 1, 1, 1, 2)) +
    guides(color = guide_legend(reverse = TRUE), linewidth = guide_legend(reverse = TRUE)) +
    scale_y_continuous(expand = c(0.0003,0)) +
    theme(#axis.text.x = element_text(angle = 90),
      #text = element_text(size = 18, family = "Arial"),
      axis.title.x = element_text(vjust=1),
      plot.tag.position = c(0, 1), 
      #plot.tag = element_text(vjust = 1.3, hjust =-0.2), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      #axis.title.y = element_blank(),
      legend.position = "none")
  print(N2024plot)
  return(N2024plot)
}
# Base Model --------------------------------------------------------------
#Posterior densities
#R

base_r_post<-subset_draws_par(draws_mat_uR_112, "r", "Base")
Base_r_prior<-runif(10000, 0, 0.115)
Base_r_prior_dens<-dunif(Base_r_prior, 0, 0.115)

Base_rpostplot<- plot_R(base_r_post, prior = TRUE, r_prior = Base_r_prior, r_prior_dens = Base_r_prior_dens, base = TRUE, r_post_base = NULL)

Base_rpostplot


#K
base_K_post<-subset_draws_par(draws_mat_uR_112, "K", "Base")
base_Kpostplot<-plot_K(base_K_post, base = TRUE, K_post_base = NULL)

#Nmin
base_Nmin_post<-subset_draws_par(draws_mat_uR_112, "Nmin", "Base")
base_Nmin_prior<-runif(10000, log(106), log(2000))
base_Nmin_prior_exp<-exp(base_Nmin_prior)
base_Nmin_prior_dens<-dunif(base_Nmin_prior, log(106), log(2000))
base_Nmin_prior_dens2<-base_Nmin_prior_dens/sum(base_Nmin_prior_dens) #normalize density

base_nmin_postplot<-plot_Nmin(base_Nmin_post, prior = TRUE, Nmin_prior = base_Nmin_prior_exp, Nmin_prior_dens = base_Nmin_prior_dens2 + 0.0005, base = TRUE, Nmin_post_base = NULL)

#add CV Branch
base_addCV_Branch_post<-subset_draws_par(draws_mat_uR_112, "addCV_Branch", "Base")
base_addCV_prior<-rexp(10000, 1)
base_addCV_prior_dens<-dexp(base_addCV_prior, 1)

base_addCV_Branch_postplot<-plot_addCV(base_addCV_Branch_post, prior = TRUE, addCV_prior = base_addCV_prior, addCV_prior_dens = base_addCV_prior_dens, source = "SOWER")

#add CV Ham
base_addCV_Ham_post<-subset_draws_par(draws_mat_uR_112, "addCV_Ham", "Base")
base_addCV_Ham_prior<-rexp(10000, 1)
base_addCV_Ham_prior_dens<-dexp(base_addCV_Ham_prior, 1)
base_addCV_Ham_postplot<-plot_addCV(base_addCV_Ham_post, prior = TRUE, addCV_prior = base_addCV_Ham_prior, addCV_prior_dens = base_addCV_Ham_prior_dens, source = "JARPA")


#Depletion (Nmin/K)
base_depl_post<-subset_draws_par(draws_mat_uR_112, "Nmin_div_K", "Base")
base_depl_plot<-plot_depletion(base_depl_post)

#N2024
base_N2024_post<-subset_draws_par(draws_mat_uR_112, "N2024", "Base")
base_N2024_plot<-plot_N2024(base_N2024_post, base = TRUE, N2024_post_base = NULL)


#plot together
base_plots<-Base_rpostplot + base_Kpostplot + base_nmin_postplot + base_depl_plot + base_addCV_Branch_postplot + base_addCV_Ham_postplot + base_N2024_plot +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect_y", tag_level = "new") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") & 
  theme(plot.tag.location = "panel", plot.tag.position = c(0.03, 1.05), text = element_text(size = 14, family = "Arial"), legend.position = "none")

base_plots

#save
ggsave("Figures/Manuscript/BaseModelPosteriors.png", base_plots, width = 8, height = 10, units = "in", dpi = 600)

#population trajectory
Npop_post<-subset_npop(draws_mat_uR_112, "Base")
write_csv(Npop_post, "./Code/Results/BaseModel_Npop_Posteriors.csv")
#annotated
FullPop_annot<-ggplot(Npop_post) + 
  geom_ribbon(aes(x = year, ymin= `2.5%`, ymax = `97.5%`),
              alpha = 0.3, fill = "deepskyblue4") + 
  geom_line(aes(x = year, y = `50%`), color = "deepskyblue4") +
  annotate("point", x = c(1904, 1973, 2024), y = c(Npop_post$`50%`[1],Npop_post$`50%`[69], Npop_post$`50%`[121]), 
           color = "deepskyblue4", size = 3) +
  annotate("text", x = c(1913, 1973, 2023), y = c(Npop_post$`50%`[1]+68000, Npop_post$`50%`[69] + 52000, Npop_post$`50%`[121] + 52000), 
           label = formatC(c(Npop_post$`50%`[1], Npop_post$`50%`[69], Npop_post$`50%`[121]), big.mark = ",", digits = 0, format = "f"), color = "deepskyblue4", size = 4, hjust = 0.5, fontface = "bold") +
  annotate("text", x = c(1914.5, 1973, 2023), y = c(Npop_post$`50%`[1]+30000,  Npop_post$`50%`[69] + 24000, Npop_post$`50%`[121] + 25000), 
           label = c("Estimated\npre-whaling\npopulation size", "Estimated\nminimum\npopulation size", "Current\npopulation\nsize"), color = "black", size = 3.3, hjust = 0.5) +
  scale_x_continuous(expand = c(0.01,0),limits = c(1904, 2029), breaks = seq(1904,2024, by = 10)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 310000), breaks = seq(0,300000, by = 50000), label = comma) +
  coord_cartesian(clip = "off") +
  theme_classic() + labs(x = "Year", y = "Antarctic blue whale population size (#s)") + 
  theme(text = element_text(size = 12, family = "Arial"), 
        plot.tag.position = c(0, 1), 
        plot.tag = element_text(hjust = -3), 
        #legend.position = "none", 
        axis.text.x = element_text(angle = 0))

FullPop_annot

ggsave("Figures/Manuscript/Poptraj.png", FullPop_annot, width = 6, height = 4, units = "in", dpi = 600)

#model fit
FullPop<-ggplot(Npop_post) + 
  geom_ribbon(aes(x = year, ymin= `2.5%`, ymax = `97.5%`),
              alpha = 0, linetype = "dashed", color = "black") + 
  geom_line(aes(x = year, y = `50%`)) +
  scale_x_continuous(expand = c(0,0),limits = c(1904, 2026.5), breaks = seq(1904,2024, by = 10)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 310000), breaks = seq(0,300000, by = 50000), label = comma) +
  theme_classic() + labs(x = "Season", y = "Population size") + 
  theme(text = element_text(size = 18, family = "Arial"), 
        plot.tag.position = c(0, 1), 
        plot.tag = element_text(hjust = -3), 
        #legend.position = "none", 
        axis.text.x = element_text(angle = 90))

#add abundance estimates

#add in additional CVs
abund_dat_plot<-addCV_func(median(base_addCV_Branch_post$addCV_Branch), median(base_addCV_Ham_post$addCV_Ham), "base")

#add missing year for JARPA
abund_dat_plot<-abund_dat_plot %>%
  add_row(Year = 1995, N = NA, CV = NA, Lower = NA, Upper = NA, Source = "Hamabe 2023", model = "base",
          C = NA, AddLower = NA, AddUpper = NA)

#zoomed in population plot with circumpolar abundance
Abund_zoom<-FullPop + coord_cartesian(xlim = c(1978, 2025.5), ylim = c(0, 6000), expand = FALSE) +
  scale_y_continuous(breaks = seq(0,6000,1000)) + 
  scale_color_manual(values = rep("black", 6)) +
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source != "Hamabe 2023",], aes(x = Year, ymin = AddLower, ymax = AddUpper), color = data_pal[1], linewidth = 0.5) +
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source != "Hamabe 2023",], aes(x = Year, ymin = Lower, ymax = Upper), color = data_pal[1], linewidth = 1.3) +
  geom_point(data = abund_dat_plot[abund_dat_plot$Source != "Hamabe 2023",], aes(x = Year, y = N), color = data_pal[1], size = 2) + 
  geom_point(data = photo_id[2,], aes(x = Year, y = N), color = data_pal[3], shape = 15, size = 3) +
  geom_linerange(data = photo_id[2,], aes(x = Year, ymin = Lower, ymax = Upper), color = data_pal[3], linewidth = 1.3) +
  labs(tag = "a)") #+ theme(plot.tag = element_text(hjust = -3))
Abund_zoom

#JARPA abundance
q_ham_post<-subset_draws_par(draws_mat_uR_112, "q_Ham", "Base")
qN_post<-subset_npop(draws_mat_uR_112*q_ham_post$q_Ham, "Base")

Adjusted_pop<-ggplot(qN_post) +
  geom_ribbon(aes(x = year, ymin= `2.5%`, ymax = `97.5%`), color = "black",
              alpha = 0, linetype = "dashed") + 
  geom_line(aes(x = year, y = `50%`)) +
  scale_x_continuous(expand = c(0,0),limits = c(1904, 2026.5), breaks = seq(1904,2024, by = 10)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 100000), breaks = seq(0,100000, by = 50000), label = comma) +
  #scale_color_manual(values = pal2) +
  theme_classic() + labs(x = "Season", y = "Adjusted population size (q*N)") + 
  theme(text = element_text(size = 18, family = "Arial"), 
        plot.tag.position = c(0, 1), 
        plot.tag = element_text(hjust = -3), 
        #legend.position = "none", 
        axis.text.x = element_text(angle = 90))
Adjusted_pop 

Jarpa_zoom<-Adjusted_pop + coord_cartesian(xlim = c(1978, 2025), ylim = c(0, 3000), expand = FALSE) +
  scale_y_continuous(breaks = seq(0,3000,500)) +
  scale_color_manual(values = rep("black", 6)) +
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source == "Hamabe 2023",], aes(x = Year, ymin = AddLower, ymax = AddUpper), color = data_pal[2], linewidth = 0.5) +
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source == "Hamabe 2023",], aes(x = Year, ymin = Lower, ymax = Upper), color = data_pal[2], linewidth = 1.3) +
  geom_point(data = abund_dat_plot[abund_dat_plot$Source == "Hamabe 2023",], aes(x = Year, y = N), color = data_pal[2], shape = 17, size = 3) + 
  labs(tag = "b)") #+ theme(plot.tag = element_text(hjust = -2.5))
Jarpa_zoom


#posterior predictive
PPout_base<-post_out(post_pred_base, "base")

PP_S_base<-PPout_base %>% filter(source == "SOWER") %>%
  ggplot() + geom_violin(aes(x = as.factor(year), y = value), scale = "width", width = 0.4, fill = data_pal[1], color = "gray70") +
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source != "Hamabe 2023",], aes(x = as.factor(Year), ymin = AddLower, ymax = AddUpper), color = "white", linewidth = 0.5) +
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source != "Hamabe 2023",], aes(x = as.factor(Year), ymin = Lower, ymax = Upper), color = "white", linewidth = 1.3) +
  geom_point(data = abund_dat_plot[abund_dat_plot$Source != "Hamabe 2023",], aes(x = as.factor(Year), y = N), color = "white") +
  theme_classic() + 
  coord_cartesian(ylim = c(0, 5000), expand = FALSE) +
  labs(x = "Year", y = "Abundance estimate") +
  theme(legend.position = "none", text = element_text(size = 18, family = "Arial")) +
 ggtitle("c) SOWER")
PP_S_base
PPout_base$Year<-factor(PPout_base$year, levels = sort(abund_dat_plot$Year[-c(1:3)]))
PP_J_base<-PPout_base%>% filter(source == "JARPA") %>%
  ggplot() + geom_violin(aes(x = Year, y = value), scale = "width", width = 0.5, fill = data_pal[2], alpha = 1, color = "gray70") +
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source == "Hamabe 2023",], aes(x = as.factor(Year), ymin = AddLower, ymax = AddUpper), color = "white", linewidth = 0.5) + 
  geom_linerange(data = abund_dat_plot[abund_dat_plot$Source == "Hamabe 2023",], aes(x = as.factor(Year), ymin = Lower, ymax = Upper), color = "white", linewidth = 1.3) +
  geom_point(data = abund_dat_plot[abund_dat_plot$Source == "Hamabe 2023",], aes(x = as.factor(Year), y = N), shape = 17, color = "white") +
  theme_classic() + 
  coord_cartesian(ylim = c(0, 3000), expand = FALSE) +
  labs(x = "Year", y = "Abundance estimate") +
  scale_x_discrete(drop = FALSE) + 
  theme(legend.position = "none", text = element_text(size = 18, family = "Arial"), 
        axis.title.y = element_blank()) +
  ggtitle("e) JARPA")

PP_J_base

PP_P_base<-PPout_base %>% filter(source == "Photo ID") %>% 
  ggplot() + geom_violin(aes(x = as.factor(year), y = value), scale = "width", width = 0.2, fill = data_pal[3], alpha = 1, color = "gray70") +
  geom_linerange(data = photo_id[2,], aes(x = as.factor(Year), ymin = Lower, ymax = Upper), color = "white") +
  geom_point(data = photo_id[2,], aes(x = as.factor(Year), y = N), shape = 8, color = "white") +
  theme_classic() + 
  coord_cartesian(ylim = c(0, 12000), expand = FALSE) +
  labs(x = "Year", y = "Abundance estimate") +
  theme(legend.position = "none", text = element_text(size = 18, family = "Arial"), 
        axis.title.y = element_blank(), 
        title = element_text(size = 18)) +
  ggtitle("d) Photo \n    ID")

PP_P_base

#plot together
lay<-"
AABCCCCCC"
PP_base_tog<-PP_S_base + PP_P_base + PP_J_base +
  plot_layout(design= lay) & theme(plot.tag = element_text(hjust = -2))
PP_base_tog

#add poseteriors on top
Post_traj_abund<-Abund_zoom + Jarpa_zoom
model_fit<-Post_traj_abund / PP_base_tog

model_fit
#save
ggsave("Figures/Manuscript/post_posterior_predictive.png", model_fit, width = 12, height = 9, units = "in", dpi = 600)

#harvest rate
harvest_rates<-Npop_post %>% filter(year< 1974) %>% 
  mutate(harvest_lwr = catch_dat$ABW_catch[-1]/`97.5%`,
         harvest_med = catch_dat$ABW_catch[-1]/`50%`,
         harvest_upr = catch_dat$ABW_catch[-1]/`2.5%`)

ggplot(harvest_rates) + 
  geom_ribbon(aes(x = year, ymin= harvest_lwr, ymax = harvest_upr),
              alpha = 0.3, fill = "deepskyblue4") + 
  geom_line(aes(x = year, y = harvest_med)) +
  scale_x_continuous(expand = c(0.01,0),limits = c(1904, 1973), breaks = seq(1904,1973, by = 10)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.5), breaks = seq(0,0.5, by = 0.1), label = percent_format(accuracy = 1)) +
  coord_cartesian(clip = "off") +
  theme_classic() + labs(x = "Whaling season", y = "Harvest rate") + 
  theme(text = element_text(size = 12, family = "Arial"), 
        plot.tag.position = c(0, 1), 
        plot.tag = element_text(hjust = -3))
ggsave("Figures/Manuscript/harvest_rate.png", width = 6, height = 4, units = "in", dpi = 600)

summary(harvest_rates$harvest_med)
# Sensitivities -----------------------------------------------------------
#plot R
sens_r_post<-mapply(subset_draws_par, draws_list, model_names_txt, MoreArgs = list(par = "r"), SIMPLIFY = FALSE)
sens_r_post<-bind_rows(sens_r_post)
base_r_post<-sens_r_post %>%
  filter(model == "Base") 
sens_r_post<-sens_r_post %>% 
  mutate(model = factor(model, levels = model_names)) %>%
  filter(model != "Base")

 
sens_r_postplot<-plot_R(sens_r_post, prior = FALSE, r_prior = NULL, r_prior_dens = NULL, base = FALSE, r_post_base = base_r_post)

#plot K

sens_K_post<-mapply(subset_draws_par, draws_list, model_names_txt, MoreArgs = list(par = "K"), SIMPLIFY = FALSE)
sens_K_post<-bind_rows(sens_K_post)
base_K_post<-sens_K_post %>%
  filter(model == "Base")
sens_K_post<-sens_K_post %>%
  filter(model != "Base") %>%
  mutate(model = factor(model, levels = model_names[-c(10,11)]))
sens_K_postplot<-plot_K(sens_K_post, base = FALSE, K_post_base = base_K_post)

#minimum N
sens_Nmin_post<-mapply(subset_draws_par, draws_list, model_names_txt, MoreArgs = list(par = "Nmin"), SIMPLIFY = FALSE)
sens_Nmin_post<-bind_rows(sens_Nmin_post)
base_Nmin_post<-sens_Nmin_post %>%
  filter(model == "Base")
sens_Nmin_post<-sens_Nmin_post %>%
  filter(model != "Base") %>%
  mutate(model = factor(model, levels = model_names))

sens_Nmin_postplot<-plot_Nmin(sens_Nmin_post, prior = FALSE, Nmin_prior = NULL, Nmin_prior_dens = NULL, base = FALSE, Nmin_post_base = base_Nmin_post)

#Population in 2024
sens_N2024_post<-mapply(subset_draws_par, draws_list, model_names_txt, MoreArgs = list(par = "N2024"), SIMPLIFY = FALSE)
sens_N2024_post<-bind_rows(sens_N2024_post)
base_N2024_post<-sens_N2024_post %>%
  filter(model == "Base")
sens_N2024_post<-sens_N2024_post %>%
  filter(model != "Base") %>%
  mutate(model = factor(model, levels = model_names))
sens_N2024_plot<-plot_N2024(sens_N2024_post, base = FALSE, N2024_post_base = base_N2024_post)


#plotting together
sens_plots<-sens_r_postplot + sens_K_postplot + sens_Nmin_postplot + sens_N2024_plot +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect_y", tag_level = "new") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") & 
  theme(plot.tag.location = "panel", plot.tag.position = c(0.05, 1), text = element_text(size = 14, family = "Arial"))
sens_plots

#save
ggsave("Figures/Manuscript/sensitivities_posteriors.png", width = 8.5, height = 5, units = "in", dpi = 600)


