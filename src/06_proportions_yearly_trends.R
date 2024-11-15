# yearly differences between END, EXO proportions

#### 1. Based on species richness ####
##### 1.1. endemic species proportion of all species #####
END_richness_prop_trends<-count_data %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_richness = n_distinct(MF)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(END_richness_prop = 
                  n_distinct(MF[establishmentmeans == "END"])/n_distinct(MF)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = END_richness_prop[project == "BALA2"] - END_richness_prop[project == "BALA1"],
              B1_B3_val_diff = END_richness_prop[project == "BALA3"] - END_richness_prop[project == "BALA1"],
              B2_B3_val_diff = END_richness_prop[project == "BALA3"] - END_richness_prop[project == "BALA2"],
              B1_B2_std_diff = (END_richness_prop[project == "BALA2"] - END_richness_prop[project == "BALA1"])/END_richness_prop[project == "BALA1"],
              B1_B3_std_diff = (END_richness_prop[project == "BALA3"] - END_richness_prop[project == "BALA1"])/END_richness_prop[project == "BALA1"],
              B2_B3_std_diff = (END_richness_prop[project == "BALA3"] - END_richness_prop[project == "BALA2"])/END_richness_prop[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

END_richness_prop_trends %>%
    select (site_code, sampling_method, 
            B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change,
            B1_B2__std_yearly_change, B1_B3__std_yearly_change, B2_B3__std_yearly_change) %>% 
    melt(.) %>% 
    separate(variable, c("comparison", "difftype"), sep = "__") %>% 
    ggplot(., aes(x = comparison, y = value, color = comparison))+
    geom_pirate()+
    theme_minimal()+
    # stat_compare_means(method = "anova", label.y = 10)+        # Add global anova p-value
    # stat_compare_means(label = "p.signif", method = "t.test",
    #                    ref.group = ".all.") +     
    facet_wrap(difftype~sampling_method, scales = "free")

ggsave("../Plots/BALA_1_3_END_richness_prop_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_END_richness_prop_differences_method.png", width = 9, height = 6)


mean_END_richness_prop_changes<-END_richness_prop_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_END_richness_prop_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_END_richness_prop_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_END_richness_prop_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_END_richness_prop_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "END_richness_prop", spp_set  = "ALL", dat))

all_trend_results<-rbind(all_trend_results, dat)


##### 1.2. exotic species proportion of all species #####
EXO_richness_prop_trends<-count_data %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_richness = n_distinct(MF)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(EXO_richness_prop = 
                  n_distinct(MF[establishmentmeans == "EXO"])/n_distinct(MF)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = EXO_richness_prop[project == "BALA2"] - EXO_richness_prop[project == "BALA1"],
              B1_B3_val_diff = EXO_richness_prop[project == "BALA3"] - EXO_richness_prop[project == "BALA1"],
              B2_B3_val_diff = EXO_richness_prop[project == "BALA3"] - EXO_richness_prop[project == "BALA2"],
              B1_B2_std_diff = (EXO_richness_prop[project == "BALA2"] - EXO_richness_prop[project == "BALA1"])/EXO_richness_prop[project == "BALA1"],
              B1_B3_std_diff = (EXO_richness_prop[project == "BALA3"] - EXO_richness_prop[project == "BALA1"])/EXO_richness_prop[project == "BALA1"],
              B2_B3_std_diff = (EXO_richness_prop[project == "BALA3"] - EXO_richness_prop[project == "BALA2"])/EXO_richness_prop[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

EXO_richness_prop_trends %>%
    select (site_code, sampling_method, 
            B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change,
            B1_B2__std_yearly_change, B1_B3__std_yearly_change, B2_B3__std_yearly_change) %>% 
    melt(.) %>% 
    separate(variable, c("comparison", "difftype"), sep = "__") %>% 
    ggplot(., aes(x = comparison, y = value, color = comparison))+
    geom_pirate()+
    theme_minimal()+
    # stat_compare_means(method = "anova", label.y = 10)+        # Add global anova p-value
    # stat_compare_means(label = "p.signif", method = "t.test",
    #                    ref.group = ".all.") +     
    facet_wrap(difftype~sampling_method, scales = "free")

ggsave("../Plots/BALA_1_3_EXO_richness_prop_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_EXO_richness_prop_differences_method.png", width = 9, height = 6)


mean_EXO_richness_prop_changes<-EXO_richness_prop_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_EXO_richness_prop_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_EXO_richness_prop_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_EXO_richness_prop_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_EXO_richness_prop_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "EXO_richness_prop", spp_set  = "ALL", dat))

all_trend_results<-rbind(all_trend_results, dat)



#### 2. Based on abundances ####
##### 2.1. the proportion of endemic species abundances to the overall abundances #####
END_abu_prop_trends<-count_data %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_abundance = sum(total_abundance)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(END_abundance_prop = 
                  sum(total_abundance[establishmentmeans == "END"])/sum(total_abundance)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = END_abundance_prop[project == "BALA2"] - END_abundance_prop[project == "BALA1"],
              B1_B3_val_diff = END_abundance_prop[project == "BALA3"] - END_abundance_prop[project == "BALA1"],
              B2_B3_val_diff = END_abundance_prop[project == "BALA3"] - END_abundance_prop[project == "BALA2"],
              B1_B2_std_diff = (END_abundance_prop[project == "BALA2"] - END_abundance_prop[project == "BALA1"])/END_abundance_prop[project == "BALA1"],
              B1_B3_std_diff = (END_abundance_prop[project == "BALA3"] - END_abundance_prop[project == "BALA1"])/END_abundance_prop[project == "BALA1"],
              B2_B3_std_diff = (END_abundance_prop[project == "BALA3"] - END_abundance_prop[project == "BALA2"])/END_abundance_prop[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

END_abu_prop_trends %>%
    select (site_code, sampling_method, 
            B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change,
            B1_B2__std_yearly_change, B1_B3__std_yearly_change, B2_B3__std_yearly_change) %>% 
    melt(.) %>% 
    separate(variable, c("comparison", "difftype"), sep = "__") %>% 
    ggplot(., aes(x = comparison, y = value, color = comparison))+
    geom_pirate()+
    theme_minimal()+
    # stat_compare_means(method = "anova", label.y = 10)+        # Add global anova p-value
    # stat_compare_means(label = "p.signif", method = "t.test",
    #                    ref.group = ".all.") +     
    facet_wrap(difftype~sampling_method, scales = "free")

ggsave("../Plots/BALA_1_3_END_abu_prop_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_END_abu_prop_differences_method.png", width = 9, height = 6)


mean_END_abu_prop_changes<-END_abu_prop_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_END_abu_prop_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_END_abu_prop_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_END_abu_prop_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_END_abu_prop_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "END_abu_prop", spp_set  = "ALL", dat))

all_trend_results<-rbind(all_trend_results, dat)


##### 2.2. the proportion of exotic species abundances to the overall abundances #####
EXO_abu_prop_trends<-count_data %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_abundance = sum(total_abundance)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(EXO_abundance_prop = 
                  sum(total_abundance[establishmentmeans == "EXO"])/sum(total_abundance)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = EXO_abundance_prop[project == "BALA2"] - EXO_abundance_prop[project == "BALA1"],
              B1_B3_val_diff = EXO_abundance_prop[project == "BALA3"] - EXO_abundance_prop[project == "BALA1"],
              B2_B3_val_diff = EXO_abundance_prop[project == "BALA3"] - EXO_abundance_prop[project == "BALA2"],
              B1_B2_std_diff = (EXO_abundance_prop[project == "BALA2"] - EXO_abundance_prop[project == "BALA1"])/EXO_abundance_prop[project == "BALA1"],
              B1_B3_std_diff = (EXO_abundance_prop[project == "BALA3"] - EXO_abundance_prop[project == "BALA1"])/EXO_abundance_prop[project == "BALA1"],
              B2_B3_std_diff = (EXO_abundance_prop[project == "BALA3"] - EXO_abundance_prop[project == "BALA2"])/EXO_abundance_prop[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

EXO_abu_prop_trends %>%
    select (site_code, sampling_method, 
            B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change,
            B1_B2__std_yearly_change, B1_B3__std_yearly_change, B2_B3__std_yearly_change) %>% 
    melt(.) %>% 
    separate(variable, c("comparison", "difftype"), sep = "__") %>% 
    ggplot(., aes(x = comparison, y = value, color = comparison))+
    geom_pirate()+
    theme_minimal()+
    # stat_compare_means(method = "anova", label.y = 10)+        # Add global anova p-value
    # stat_compare_means(label = "p.signif", method = "t.test",
    #                    ref.group = ".all.") +     
    facet_wrap(difftype~sampling_method, scales = "free")

ggsave("../Plots/BALA_1_3_EXO_abu_prop_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_EXO_abu_prop_differences_method.png", width = 9, height = 6)


mean_EXO_abu_prop_changes<-EXO_abu_prop_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_EXO_abu_prop_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_EXO_abu_prop_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_EXO_abu_prop_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_EXO_abu_prop_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "EXO_abu_prop", spp_set  = "ALL", dat))

all_trend_results<-rbind(all_trend_results, dat)



#### 2. Based on biomasses ####
##### 2.1. the proportion of endemic species' biomass to the overall biomass #####
END_biomass_prop_trends<-biomass_data %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_biomass = sum(adult_mass)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(END_biomass_prop = 
                  sum(total_biomass[establishmentmeans == "END"])/sum(total_biomass)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = END_biomass_prop[project == "BALA2"] - END_biomass_prop[project == "BALA1"],
              B1_B3_val_diff = END_biomass_prop[project == "BALA3"] - END_biomass_prop[project == "BALA1"],
              B2_B3_val_diff = END_biomass_prop[project == "BALA3"] - END_biomass_prop[project == "BALA2"],
              B1_B2_std_diff = (END_biomass_prop[project == "BALA2"] - END_biomass_prop[project == "BALA1"])/END_biomass_prop[project == "BALA1"],
              B1_B3_std_diff = (END_biomass_prop[project == "BALA3"] - END_biomass_prop[project == "BALA1"])/END_biomass_prop[project == "BALA1"],
              B2_B3_std_diff = (END_biomass_prop[project == "BALA3"] - END_biomass_prop[project == "BALA2"])/END_biomass_prop[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

END_biomass_prop_trends %>%
    select (site_code, sampling_method, 
            B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change,
            B1_B2__std_yearly_change, B1_B3__std_yearly_change, B2_B3__std_yearly_change) %>% 
    melt(.) %>% 
    separate(variable, c("comparison", "difftype"), sep = "__") %>% 
    ggplot(., aes(x = comparison, y = value, color = comparison))+
    geom_pirate()+
    theme_minimal()+
    # stat_compare_means(method = "anova", label.y = 10)+        # Add global anova p-value
    # stat_compare_means(label = "p.signif", method = "t.test",
    #                    ref.group = ".all.") +     
    facet_wrap(difftype~sampling_method, scales = "free")

ggsave("../Plots/BALA_1_3_END_biomass_prop_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_END_biomass_prop_differences_method.png", width = 9, height = 6)


mean_END_biomass_prop_changes<-END_biomass_prop_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_END_biomass_prop_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_END_biomass_prop_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_END_biomass_prop_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_END_biomass_prop_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "END_biomass_prop", spp_set  = "ALL", dat))

all_trend_results<-rbind(all_trend_results, dat)


##### 2.2. the proportion of exotic species' biomass to the overall biomass #####
EXO_biomass_prop_trends<-biomass_data %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_biomass = sum(adult_mass)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(EXO_biomass_prop = 
                  sum(total_biomass[establishmentmeans == "EXO"])/sum(total_biomass)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = EXO_biomass_prop[project == "BALA2"] - EXO_biomass_prop[project == "BALA1"],
              B1_B3_val_diff = EXO_biomass_prop[project == "BALA3"] - EXO_biomass_prop[project == "BALA1"],
              B2_B3_val_diff = EXO_biomass_prop[project == "BALA3"] - EXO_biomass_prop[project == "BALA2"],
              B1_B2_std_diff = (EXO_biomass_prop[project == "BALA2"] - EXO_biomass_prop[project == "BALA1"])/EXO_biomass_prop[project == "BALA1"],
              B1_B3_std_diff = (EXO_biomass_prop[project == "BALA3"] - EXO_biomass_prop[project == "BALA1"])/EXO_biomass_prop[project == "BALA1"],
              B2_B3_std_diff = (EXO_biomass_prop[project == "BALA3"] - EXO_biomass_prop[project == "BALA2"])/EXO_biomass_prop[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

EXO_biomass_prop_trends %>%
    select (site_code, sampling_method, 
            B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change,
            B1_B2__std_yearly_change, B1_B3__std_yearly_change, B2_B3__std_yearly_change) %>% 
    melt(.) %>% 
    separate(variable, c("comparison", "difftype"), sep = "__") %>% 
    ggplot(., aes(x = comparison, y = value, color = comparison))+
    geom_pirate()+
    theme_minimal()+
    # stat_compare_means(method = "anova", label.y = 10)+        # Add global anova p-value
    # stat_compare_means(label = "p.signif", method = "t.test",
    #                    ref.group = ".all.") +     
    facet_wrap(difftype~sampling_method, scales = "free")

ggsave("../Plots/BALA_1_3_EXO_biomass_prop_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_EXO_biomass_prop_differences_method.png", width = 9, height = 6)


mean_EXO_biomass_prop_changes<-EXO_biomass_prop_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_EXO_biomass_prop_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_EXO_biomass_prop_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_EXO_biomass_prop_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_EXO_biomass_prop_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "EXO_biomass_prop", spp_set  = "ALL", dat))

all_trend_results<-rbind(all_trend_results, dat)
