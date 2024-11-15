# yearly differences between species richnesses


#### 1. for all species ####
all_richness_trends<-count_data %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_richness = n_distinct(MF)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(richness = sum(total_richness)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = richness[project == "BALA2"] - richness[project == "BALA1"],
              B1_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA1"],
              B2_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA2"],
              B1_B2_std_diff = (richness[project == "BALA2"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B1_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B2_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA2"])/richness[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

all_richness_trends %>%
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

ggsave("../Plots/BALA_1_3_all_richness_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_all_richness_differences_method.png", width = 9, height = 6)


mean_all_richness_changes<-all_richness_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_all_richness_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_all_richness_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_all_richness_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_all_richness_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "richness", spp_set  = "ALL", dat))

all_trend_results<-rbind(all_trend_results, dat)





#### 2. for endemic species ####
END_richness_trends<-count_data %>% 
    filter(establishmentmeans=="END") %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_richness = n_distinct(MF)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(richness = sum(total_richness)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = richness[project == "BALA2"] - richness[project == "BALA1"],
              B1_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA1"],
              B2_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA2"],
              B1_B2_std_diff = (richness[project == "BALA2"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B1_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B2_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA2"])/richness[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

END_richness_trends %>%
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

ggsave("../Plots/BALA_1_3_END_richness_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_END_richness_differences_method.png", width = 9, height = 6)


mean_END_richness_changes<-END_richness_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_END_richness_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_END_richness_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_END_richness_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_END_richness_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "richness", spp_set  = "END", dat))

all_trend_results<-rbind(all_trend_results, dat)


#### 3. for native species ####
NAT_richness_trends<-count_data %>% 
    filter(establishmentmeans=="NAT") %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_richness = n_distinct(MF)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(richness = sum(total_richness)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = richness[project == "BALA2"] - richness[project == "BALA1"],
              B1_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA1"],
              B2_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA2"],
              B1_B2_std_diff = (richness[project == "BALA2"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B1_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B2_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA2"])/richness[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

NAT_richness_trends %>%
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

ggsave("../Plots/BALA_1_3_NAT_richness_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_NAT_richness_differences_method.png", width = 9, height = 6)


mean_NAT_richness_changes<-NAT_richness_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_NAT_richness_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_NAT_richness_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_NAT_richness_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_NAT_richness_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T)%>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "richness", spp_set  = "NAT", dat))

all_trend_results<-rbind(all_trend_results, dat)



#### 4. for exotic species ####
EXO_richness_trends<-count_data %>% 
    filter(establishmentmeans=="EXO") %>% 
    group_by(project, year, island, fragment_name, site_code, sampling_method,MF, establishmentmeans) %>% 
    summarise(total_richness = n_distinct(MF)) %>% 
    group_by(site_code, project, year, sampling_method) %>% 
    summarise(richness = sum(total_richness)) %>% 
    group_by(site_code, sampling_method) %>% 
    summarise(B1_B2_val_diff = richness[project == "BALA2"] - richness[project == "BALA1"],
              B1_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA1"],
              B2_B3_val_diff = richness[project == "BALA3"] - richness[project == "BALA2"],
              B1_B2_std_diff = (richness[project == "BALA2"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B1_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA1"])/richness[project == "BALA1"],
              B2_B3_std_diff = (richness[project == "BALA3"] - richness[project == "BALA2"])/richness[project == "BALA2"],
              B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
              B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
              B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]) %>% 
    mutate(B1_B2__yearly_change = B1_B2_val_diff/B1_B2_year_diff,
           B1_B3__yearly_change = B1_B3_val_diff/B1_B3_year_diff,
           B2_B3__yearly_change = B2_B3_val_diff/B2_B3_year_diff,
           B1_B2__std_yearly_change = B1_B2_std_diff/B1_B2_year_diff,
           B1_B3__std_yearly_change = B1_B3_std_diff/B1_B3_year_diff,
           B2_B3__std_yearly_change = B2_B3_std_diff/B2_B3_year_diff)

EXO_richness_trends %>%
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

ggsave("../Plots/BALA_1_3_EXO_richness_differences_method.pdf", width = 15, height = 10)
ggsave("../Plots/BALA_1_3_EXO_richness_differences_method.png", width = 9, height = 6)


mean_EXO_richness_changes<-EXO_richness_trends %>%
    ungroup() %>%
    select(sampling_method, site_code,B1_B2__yearly_change, B1_B3__yearly_change, B2_B3__yearly_change) %>% 
    group_by(sampling_method, site_code) %>% 
    summarise_all(mean) %>% 
    mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>% 
    group_by(sampling_method) %>% 
    select(-site_code) 


p_val<-mean_EXO_richness_changes %>% 
    summarise(across(where(is.numeric), is_zero)) %>% 
    melt(.)

B1_B2<-mean_EXO_richness_changes %>%     
    cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B1_B3<-mean_EXO_richness_changes %>%     
    cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

B2_B3<-mean_EXO_richness_changes %>%     
    cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T)%>% 
    select(`.y.`, sampling_method, effsize, magnitude, conf.low, conf.high) %>% 
    mutate(comparison = substr(.y., 1, 5)) %>% 
    select(comparison, sampling_method, effsize, magnitude, conf.low, conf.high)

dat<-rbind(B1_B2, B1_B3, B2_B3)

dat$p_val<-p_val$value

(dat<-data.frame(measure = "richness", spp_set  = "EXO", dat))

all_trend_results<-rbind(all_trend_results, dat)

head(all_trend_results)
