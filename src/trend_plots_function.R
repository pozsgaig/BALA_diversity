# calculates trend values for each measure, after filtering for method and biogeographic 
# origin

trend_calc <- function(imp_dat, measure, magn_alpha = NULL, 
                       cols = c("#d55e00","#0072b2", "cornsilk4")) {
    # imp_dat<-all_dat$ALL
    if (is.null(magn_alpha)){
	magn_alpha <-
    structure(c(0.1, 0.3, 0.6, 1), 
	names = c("negligible", "small", "moderate", "large"))}
	
    sum_dat<-imp_dat %>%
        
        {
            if (measure == "abundance") {
                summarise(., val = sum(total_abundance))
            }
            
            else if (measure == "richness") {
                summarise(., val = n_distinct(MF))
            }
            
            else if (measure == "biomass") {
                summarise(., val = sum(adult_mass))
            }
            
        } %>%
        
        group_by(site_code, project, year, sampling_method) %>%
        summarise(val = sum(val)) %>%
        group_by(site_code, sampling_method) %>%
        summarise(
            B1_B2_val_diff = val[project == "BALA2"] - val[project == "BALA1"],
            B1_B3_val_diff = val[project == "BALA3"] - val[project == "BALA1"],
            B2_B3_val_diff = val[project == "BALA3"] - val[project == "BALA2"],
            B1_B2_std_diff = (val[project == "BALA2"] - val[project == "BALA1"]) /
                val[project == "BALA1"],
            B1_B3_std_diff = (val[project == "BALA3"] - val[project == "BALA1"]) /
                val[project == "BALA1"],
            B2_B3_std_diff = (val[project == "BALA3"] - val[project == "BALA2"]) /
                val[project == "BALA2"],
            B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
            B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
            B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]
        ) %>%
        mutate(
            B1_B2__yearly_change = B1_B2_val_diff / B1_B2_year_diff,
            B1_B3__yearly_change = B1_B3_val_diff / B1_B3_year_diff,
            B2_B3__yearly_change = B2_B3_val_diff / B2_B3_year_diff,
            B1_B2__std_yearly_change = B1_B2_std_diff / B1_B2_year_diff,
            B1_B3__std_yearly_change = B1_B3_std_diff / B1_B3_year_diff,
            B2_B3__std_yearly_change = B2_B3_std_diff / B2_B3_year_diff
        ) %>%
        ungroup() %>%
        select(
            sampling_method,
            site_code,
            B1_B2__yearly_change,
            B1_B3__yearly_change,
            B2_B3__yearly_change
        ) %>%
        group_by(sampling_method, site_code) %>%
        summarise_all(mean) %>%
        mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>%
        group_by(sampling_method) %>%
        select(-site_code)
    
    
    B1_B2 <- sum_dat %>%
        cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>%
        select(`.y.`,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>%
        mutate(comparison = substr(.y., 1, 5)) %>%
        select(comparison,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>% 
        mutate(is_signif = conf.low  > 0 & conf.high  > 0 |
                           conf.low < 0 & conf.high < 0) %>% 
        mutate(col = ifelse(is_signif, ifelse(effsize>0, cols[1], cols[2]), cols[3])) %>% 
        mutate(alpha = ifelse(is_signif, magn_alpha[magnitude], 0.1))
    
    B1_B3 <- sum_dat %>%
        cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>%
        select(`.y.`,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>%
        mutate(comparison = substr(.y., 1, 5)) %>%
        select(comparison,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>% 
        mutate(is_signif = conf.low  > 0 & conf.high  > 0 |
                   conf.low < 0 & conf.high < 0) %>% 
        mutate(col = ifelse(is_signif, ifelse(effsize>0, cols[1], cols[2]), cols[3]))%>% 
        mutate(alpha = ifelse(is_signif, magn_alpha[magnitude], 0.1))
    
    B2_B3 <- sum_dat %>%
        cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>%
        select(`.y.`,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>%
        mutate(comparison = substr(.y., 1, 5)) %>%
        select(comparison,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>% 
        mutate(is_signif = conf.low  > 0 & conf.high  > 0 |
                   conf.low < 0 & conf.high < 0) %>% 
        mutate(col = ifelse(is_signif, ifelse(effsize>0, cols[1], cols[2]), cols[3]))%>% 
        mutate(alpha = ifelse(is_signif, magn_alpha[magnitude], 0.1))
    
    
    rbind(B1_B2, B1_B3, B2_B3)
    
    
    
}

# asseblies the complex figure for trends
trend_plots <- function(dat, measure = c("abundance", "biomass", "richness"), method,
                        cols = c("#d55e00","#0072b2", "cornsilk4"), origin_cols = c("#FF0000", "#00A08A", "darkgrey", "#F98400"))
    {
        # dat = count_data
        # measure="Richness"
        # method = "Beating"
        
        all_dat <- lapply(c("END", "NAT", "EXO"), function(x) {
            dat %>%
                group_by(
                    project,
                    year,
                    island,
                    fragment_name,
                    site_code,
                    sampling_method,
                    MF,
                    establishmentmeans
                ) %>%
                filter(sampling_method == method & establishmentmeans == x)
            
        })
        
        names(all_dat) <- c("END", "NAT", "EXO")
        all_dat[["ALL"]] <- dat %>%
            group_by(
                project,
                year,
                island,
                fragment_name,
                site_code,
                sampling_method,
                MF,
                establishmentmeans
            ) %>%
            filter(sampling_method == method)
        
        trend_dat<-lapply(all_dat, trend_calc, measure)
        
        gp3 <- data.frame(x = c(-4.5, -4.5, -3),
                          y = c(1, -3.5, -3.5))
        
        gp2 <- data.frame(x = c(-4, -4, -3),
                          y = c(1, -3, -3))
        
        gp1 <- data.frame(x = c(-3.5, -3.5, -3),
                          y = c(1, -2.5, -2.5))
        
        side_arrows <- list(END = gp1, NAT = gp2, EXO = gp3)
        horiz_arrows <-
            structure(c(3.5, 3, 2.5), names = c("END", "NAT", "EXO"))
        

        
        # reordering factor to match colours
        
        # Create barplots
        plot1 <- dat %>%
            filter(project == "BALA1" & sampling_method == method) %>%
            group_by(establishmentmeans) %>%
            {if (measure == "abundance"){summarise(., val = sum(total_abundance))}
             else if (measure == "richness") {summarise(., val = n_distinct(MF))}
             else if (measure == "biomass"){summarise(., val = sum(adult_mass))}} %>% 
            ggplot(., aes(x = "", y = val, fill = establishmentmeans)) +
            geom_bar(position = "stack", stat = "identity") +
            scale_fill_manual(values = origin_cols) +
            theme_void() +
            theme(legend.position = "none")
        
        plot2 <- dat %>%
            filter(project == "BALA2" & sampling_method == method) %>%
            group_by(establishmentmeans) %>%
            {if (measure == "abundance"){summarise(., val = sum(total_abundance))}
                else if (measure == "richness") {summarise(., val = n_distinct(MF))}
                else if (measure == "biomass"){summarise(., val = sum(adult_mass))}} %>% 
            ggplot(., aes(x = "", y = val, fill = establishmentmeans)) +
            geom_bar(position = "stack", stat = "identity") +
            scale_fill_manual(values = origin_cols) +
            theme_void() +
            theme(legend.position = "none")
        
        plot3 <- dat %>%
            filter(project == "BALA3" & sampling_method == method) %>%
            group_by(establishmentmeans) %>%
            {if (measure == "abundance"){summarise(., val = sum(total_abundance))}
                else if (measure == "richness") {summarise(., val = n_distinct(MF))}
                else if (measure == "biomass"){summarise(., val = sum(adult_mass))}} %>% 
            ggplot(., aes(x = "", y = val, fill = establishmentmeans)) +
            geom_bar(position = "stack", stat = "identity") +
            scale_fill_manual(values = origin_cols) +
            theme_void() +
            theme(legend.position = "none")
        
        
        main_plot <- ggplot() +
            xlim(-6, 6) + ylim(-6, 6) +
            ggtitle(paste(method, measure)) +
            theme_void()+
            theme(plot.title = element_text(hjust = 0.5))
            
        
        
        main_plot +
            # arrow between B1-B3
            geom_arrow_segment(
                aes(
                    x = -3,
                    y = -3,
                    xend = -1,
                    yend = -3
                ),
                colour = unlist(trend_dat$ALL[trend_dat$ALL$comparison=="B1_B3", "col"]),
                length = 1,
                arrow_head = arrow_head_wings(offset = 50, inset = 40),
                linewidth = 20,
                alpha = unlist(trend_dat$ALL[trend_dat$ALL$comparison=="B1_B3", "alpha"]),
                lineend = "butt"
            ) +
            
            # arrow between B2-B3
            geom_arrow_segment(
                aes(
                    x = 3,
                    y = -3,
                    xend = 1,
                    yend = -3
                ),
                colour = unlist(trend_dat$ALL[trend_dat$ALL$comparison=="B2_B3", "col"]),
                length = 1,
                arrow_head = arrow_head_wings(offset = 50, inset = 40),
                linewidth = 20,
                alpha = unlist(trend_dat$ALL[trend_dat$ALL$comparison=="B2_B3", "alpha"]),
                lineend = "butt"
            ) +
            
            # arrow between B1-B2
            geom_arrow_segment(
                aes(
                    x = 0,
                    y = 3,
                    xend = 3,
                    yend = 3
                ),
                colour = unlist(trend_dat$ALL[trend_dat$ALL$comparison=="B1_B2", "col"]),
                length = 1,
                arrow_head = arrow_head_wings(offset = 50, inset = 40),
                linewidth = 20,
                alpha = unlist(trend_dat$ALL[trend_dat$ALL$comparison=="B1_B2", "alpha"]),
                lineend = "butt",
                linetype = 2
            ) +
            
            annotation_custom(
                ggplotGrob(plot1),
                xmin = -5,
                xmax = -3,
                ymin = 1,
                ymax = 5
            ) +
            annotate("text",
                     label = "B1",
                     x = -4,
                     y = 5.1) +
            annotation_custom(
                ggplotGrob(plot2),
                xmin = 3,
                xmax = 5,
                ymin = 1,
                ymax = 5
            ) +
            annotate("text",
                     label = "B2",
                     x = 4,
                     y = 5.1) +
            annotation_custom(
                ggplotGrob(plot3),
                xmin = -1,
                xmax = 1,
                ymin = -1,
                ymax = -5
            ) +
            annotate("text",
                     label = "B3",
                     x = 0,
                     y = -0.9) +
            
            # origin arrow between B1-B2
            lapply(names(trend_dat)[1:3], function(x) {
                #x = "END"
                
                origin_dat <- trend_dat[[x]] %>%
                    filter(comparison == "B1_B2")
                
                geom_segment(
                    aes(x = -3,
                        y = horiz_arrows[x],
                        xend = 0,
                        yend = horiz_arrows[x]),
                    colour = origin_dat$col,
                    linewidth = 3,
                    alpha =origin_dat$alpha,
                    lineend = "butt"
                )}) +
            
            # origin points between B1-B2 start    
            lapply(names(trend_dat)[1:3], function(x) {geom_point(
                    aes(x = -3,
                        y = horiz_arrows[x]),
                    colour = origin_cols[x])
            }) +
            # origin arrow between B1-B2 stop
            lapply(names(trend_dat)[1:3], function(x) {geom_point(
                aes(x = 0,
                    y = horiz_arrows[x]),
                colour = origin_cols[x])
            }) +
            
            # origin arrow between B1-B3
            lapply(names(trend_dat)[1:3], function(x) {
                
                origin_dat <- trend_dat[[x]] %>%
                    filter(comparison == "B1_B3")
                

                geom_path(
                    aes(x = side_arrows[[x]]$x,
                        y = side_arrows[[x]]$y),
                    colour = origin_dat$col,
                    linewidth = 3,
                    alpha = origin_dat$alpha,
                    lineend = "butt"
                )
            }) +
            
            # origin points between B1-B3  
            lapply(names(trend_dat)[1:3], function(x) {geom_point(
                aes(x = side_arrows[[x]]$x,
                    y = side_arrows[[x]]$y),
                colour = origin_cols[x])
            }) +
            
            # origin arrow between B2-B3
            lapply(names(trend_dat)[1:3], function(x) {
                # x = "END"
                origin_dat <- trend_dat[[x]] %>%
                    filter(comparison == "B2_B3")

                geom_path(
                    aes(x = side_arrows[[x]]$x * -1,
                        y = side_arrows[[x]]$y),
                    colour = origin_dat$col,
                    linewidth = 3,
                    alpha = origin_dat$alpha,
                    lineend = "butt"
                )}) +
             
            # origin points between B2-B3         
            lapply(names(trend_dat)[1:3], function(x) {geom_point(
                        aes(x = side_arrows[[x]]$x* -1,
                            y = side_arrows[[x]]$y),
                        colour = origin_cols[x])
                    })
        
        
    }

### functions for calculating and plotting endemic and introduced trends

trend_calc_prop <- function(imp_dat, measure, prop, magn_alpha = NULL, 
                            cols = c("#d55e00","#0072b2", "cornsilk4")) {
    # imp_dat<-all_dat
    # prop = "END"
    if (is.null(magn_alpha)){
        magn_alpha <-
            structure(c(0.1, 0.3, 0.6, 1), 
                      names = c("negligible", "small", "moderate", "large"))}
    sum_dat<-imp_dat %>%
        
        {
            if (measure == "abundance") {
                summarise(., val = sum(total_abundance))
            }
            
            else if (measure == "richness") {
                summarise(., val = n_distinct(MF))
            }
            
            else if (measure == "biomass") {
                summarise(., val = sum(adult_mass))
            }
            
        } %>%
        
        group_by(site_code, project, year, sampling_method) %>%
        
        {
            if (measure == "richness") {
                summarise(., val = 
                              n_distinct(MF[establishmentmeans == prop])/n_distinct(MF)) 
            }
            
            else if (measure == "abundance") {
                summarise(., val = 
                              sum(val[establishmentmeans == prop])/sum(val))
            }
            
            else if (measure == "biomass") {
                summarise(., val = 
                              sum(val[establishmentmeans == prop])/sum(val))
            }
            
        } %>%
        group_by(site_code, sampling_method) %>%
        summarise(
            B1_B2_val_diff = val[project == "BALA2"] - val[project == "BALA1"],
            B1_B3_val_diff = val[project == "BALA3"] - val[project == "BALA1"],
            B2_B3_val_diff = val[project == "BALA3"] - val[project == "BALA2"],
            B1_B2_std_diff = (val[project == "BALA2"] - val[project == "BALA1"]) /
                val[project == "BALA1"],
            B1_B3_std_diff = (val[project == "BALA3"] - val[project == "BALA1"]) /
                val[project == "BALA1"],
            B2_B3_std_diff = (val[project == "BALA3"] - val[project == "BALA2"]) /
                val[project == "BALA2"],
            B1_B2_year_diff = year[project == "BALA2"] - year[project == "BALA1"],
            B1_B3_year_diff = year[project == "BALA3"] - year[project == "BALA1"],
            B2_B3_year_diff = year[project == "BALA3"] - year[project == "BALA2"]
        ) %>%
        mutate(
            B1_B2__yearly_change = B1_B2_val_diff / B1_B2_year_diff,
            B1_B3__yearly_change = B1_B3_val_diff / B1_B3_year_diff,
            B2_B3__yearly_change = B2_B3_val_diff / B2_B3_year_diff,
            B1_B2__std_yearly_change = B1_B2_std_diff / B1_B2_year_diff,
            B1_B3__std_yearly_change = B1_B3_std_diff / B1_B3_year_diff,
            B2_B3__std_yearly_change = B2_B3_std_diff / B2_B3_year_diff
        ) %>%
        ungroup() %>%
        select(
            sampling_method,
            site_code,
            B1_B2__yearly_change,
            B1_B3__yearly_change,
            B2_B3__yearly_change
        ) %>%
        group_by(sampling_method, site_code) %>%
        summarise_all(mean) %>%
        mutate(across(where(is.numeric), \(x) scale(x, center = F))) %>%
        group_by(sampling_method) %>%
        select(-site_code)
    
    
    B1_B2 <- sum_dat %>%
        cohens_d(B1_B2__yearly_change ~ 1, mu = 0, ci = T) %>%
        select(`.y.`,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>%
        mutate(comparison = substr(.y., 1, 5)) %>%
        select(comparison,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>% 
        mutate(is_signif = conf.low  > 0 & conf.high  > 0 |
                   conf.low < 0 & conf.high < 0) %>% 
        mutate(col = ifelse(is_signif, ifelse(effsize>0, cols[1], cols[2]), cols[3])) %>% 
        mutate(alpha = ifelse(is_signif, magn_alpha[magnitude], 0.1))
    
    B1_B3 <- sum_dat %>%
        cohens_d(B1_B3__yearly_change ~ 1, mu = 0, ci = T) %>%
        select(`.y.`,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>%
        mutate(comparison = substr(.y., 1, 5)) %>%
        select(comparison,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>% 
        mutate(is_signif = conf.low  > 0 & conf.high  > 0 |
                   conf.low < 0 & conf.high < 0) %>% 
        mutate(col = ifelse(is_signif, ifelse(effsize>0, cols[1], cols[2]), cols[3]))%>% 
        mutate(alpha = ifelse(is_signif, magn_alpha[magnitude], 0.1))
    
    B2_B3 <- sum_dat %>%
        cohens_d(B2_B3__yearly_change ~ 1, mu = 0, ci = T) %>%
        select(`.y.`,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>%
        mutate(comparison = substr(.y., 1, 5)) %>%
        select(comparison,
               sampling_method,
               effsize,
               magnitude,
               conf.low,
               conf.high) %>% 
        mutate(is_signif = conf.low  > 0 & conf.high  > 0 |
                   conf.low < 0 & conf.high < 0) %>% 
        mutate(col = ifelse(is_signif, ifelse(effsize>0, cols[1], cols[2]), cols[3]))%>% 
        mutate(alpha = ifelse(is_signif, magn_alpha[magnitude], 0.1))
    
    
    rbind(B1_B2, B1_B3, B2_B3)
    
    
    
}


trend_plots_prop <- function(dat, measure = c("abundance", "richness", "biomass"), 
                             method, cols = c("#d55e00","#0072b2", "cornsilk4"),
                             origin_cols = c("#00A08A", "#FF0000", "#F98400", "darkgrey"))
{
    # dat = count_data
    # measure="Abundance"
    # method = "Beating"
    # prop = "END"
    # rm(prop)
    
    
    all_dat <- dat %>%
        group_by(
            project,
            year,
            island,
            fragment_name,
            site_code,
            sampling_method,
            MF,
            establishmentmeans
        ) %>%
        filter(sampling_method == method)
    
    trend_dat<-lapply(c("END", "EXO"), function(x) {
        trend_calc_prop(all_dat, measure = measure, 
                        prop = x)})
    names(trend_dat)<- c("END", "EXO")
    
    gp3 <- data.frame(x = c(-4.5, -4.5, -3),
                      y = c(1, -3.5, -3.5))
    
    gp2 <- data.frame(x = c(-4, -4, -3),
                      y = c(1, -3, -3))
    
    gp1 <- data.frame(x = c(-3.5, -3.5, -3),
                      y = c(1, -2.5, -2.5))
    
    side_arrows <- list(END = gp1, NAT = gp2, EXO = gp3)
    horiz_arrows <-
        structure(c(3.5, 3, 2.5), names = c("END", "NAT", "EXO"))
    
    dat$establishmentmeans<-factor(dat$establishmentmeans, 
                                   levels = c("END", "NAT", "EXO", "UnK"))
    
    # reordering factor to match colours
    
    # Create barplots
    plot1 <- dat %>%
        filter(project == "BALA1" & sampling_method == method) %>%
        group_by(establishmentmeans) %>%
        {if (measure == "abundance"){summarise(., val = sum(total_abundance))}
            else if (measure == "richness") {summarise(., val = n_distinct(MF))}
            else if (measure == "biomass"){summarise(., val = sum(adult_mass))}} %>% 
        mutate(val = val/sum(val)) %>% 
        filter(establishmentmeans %in% c("END", "EXO")) %>% 
        ggplot(., aes(x = "", y = val, fill = establishmentmeans)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_fill_manual(values = origin_cols) +
        theme_void() +
        theme(legend.position = "none")
    
    plot2 <- dat %>%
        filter(project == "BALA2" & sampling_method == method) %>%
        group_by(establishmentmeans) %>%
        {if (measure == "abundance"){summarise(., val = sum(total_abundance))}
            else if (measure == "richness") {summarise(., val = n_distinct(MF))}
            else if (measure == "biomass"){summarise(., val = sum(adult_mass))}} %>% 
        mutate(val = val/sum(val)) %>% 
        filter(establishmentmeans %in% c("END", "EXO")) %>% 
        ggplot(., aes(x = "", y = val, fill = establishmentmeans)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_fill_manual(values = origin_cols) +
        theme_void() +
        theme(legend.position = "none")
    
    plot3 <- dat %>%
        filter(project == "BALA3" & sampling_method == method) %>%
        group_by(establishmentmeans) %>%
        {if (measure == "abundance"){summarise(., val = sum(total_abundance))}
            else if (measure == "richness") {summarise(., val = n_distinct(MF))}
            else if (measure == "biomass"){summarise(., val = sum(adult_mass))}} %>% 
        mutate(val = val/sum(val)) %>% 
        filter(establishmentmeans %in% c("END", "EXO")) %>% 
        ggplot(., aes(x = "", y = val, fill = establishmentmeans)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_fill_manual(values = origin_cols) +
        theme_void() +
        theme(legend.position = "none")
    
    
    main_plot <- ggplot() +
        xlim(-6, 6) + ylim(-6, 6) +
        ggtitle(paste(method, measure)) +
        theme_void()+
        theme(plot.title = element_text(hjust = 0.5))
    
    
    
    main_plot +
        #### Barplots ####
    # B1 barplot an text above       
    annotation_custom(
        ggplotGrob(plot1),
        xmin = -5,xmax = -3, ymin = 1, ymax = 5) +
        annotate("text", label = "B1", x = -4, y = 5.1) +
        
        # B2 barplot an text above 
        annotation_custom(
            ggplotGrob(plot2), xmin = 3, xmax = 5, ymin = 1, ymax = 5) +
        annotate("text", label = "B2", x = 4, y = 5.1) +
        
        # B3 barplot an text above 
        annotation_custom(
            ggplotGrob(plot3),
            xmin = -1, xmax = 1, ymin = -1, ymax = -5) +
        annotate("text", label = "B3", x = 0, y = -0.9) +
        
        #### B1-B2 comparison arrows and points ####     
    
    # origin arrow between B1-B2
    lapply(names(trend_dat), function(x) {
        #x = "END"
        
        origin_dat <- trend_dat[[x]] %>%
            filter(comparison == "B1_B2")
        
        geom_arrow_segment(
            aes(x = -3, y = horiz_arrows[x], xend = 3, yend = horiz_arrows[x]),
            colour = origin_dat$col,
            linewidth = 3,
            alpha =origin_dat$alpha,
            lineend = "butt",
            length = 2
        )}) +
        
        # origin points for B1-B2 start   
        lapply(names(trend_dat), function(x) {geom_point(
            aes(x = -3, y = horiz_arrows[x]),
            colour = origin_cols[x])
        }) +
        
        # origin points for B1-B2 stop
        lapply(names(trend_dat), function(x) {geom_point(
            aes(x = 3,
                y = horiz_arrows[x]),
            colour = origin_cols[x])
        }) +
        
        #### B2-B3 comparison arrows and points ####        
    
    # origin path between B2-B3
    lapply(names(trend_dat), function(x) {
        # x = "END"
        origin_dat <- trend_dat[[x]] %>%
            filter(comparison == "B2_B3")
        
        geom_path(
            aes(x = -side_arrows[[x]]$x,
                y = side_arrows[[x]]$y),
            colour = origin_dat$col,
            linewidth = 3,
            alpha = origin_dat$alpha,
            lineend = "butt"
        )}) +
        
        # tip arrow between B2-B3
        lapply(names(trend_dat), function(x) {
            #x = "END"
            
            origin_dat <- trend_dat[[x]] %>%
                filter(comparison == "B2_B3")
            
            geom_arrow_segment(
                # setting y coords to  minus values without mirroring them 
                # (6 is the height of the plot)
                aes(x = 3, y = horiz_arrows[x]-6, xend = 1, yend = horiz_arrows[x]-6),
                colour = origin_dat$col,
                length = 2,
                linewidth = 3,
                alpha = origin_dat$alpha,
                lineend = "butt"
            )
        }) +
        
        
        # origin points between B2-B3         
        lapply(names(trend_dat), function(x) {
            x_coord = side_arrows[[x]]$x
            x_coord[3] = x_coord[3]+2
            geom_point(
                aes(x = -x_coord,
                    y = side_arrows[[x]]$y),
                colour = origin_cols[x])
        }) +
        
        
        
        #### B1-B3 comparison arrows and points ####                
    
    # origin path between B1-B3
    lapply(names(trend_dat), function(x) {
        
        origin_dat <- trend_dat[[x]] %>%
            filter(comparison == "B1_B3")
        
        
        geom_path(
            aes(x = side_arrows[[x]]$x,
                y = side_arrows[[x]]$y),
            colour = origin_dat$col,
            linewidth = 3,
            alpha = origin_dat$alpha,
            lineend = "butt"
        )
    }) +
        
        # tip arrow between B1-B3
        lapply(names(trend_dat), function(x) {
            #x = "END"
            
            origin_dat <- trend_dat[[x]] %>%
                filter(comparison == "B1_B3")
            
            geom_arrow_segment(
                # setting y coords to  minus values without mirroring them 
                # (6 is the height of the plot)
                aes(x = -3, y = horiz_arrows[x]-6, xend = -1, yend = horiz_arrows[x]-6),
                colour = origin_dat$col,
                length = 2,
                linewidth = 3,
                alpha = origin_dat$alpha,
                lineend = "butt"
            )
        })+
        
        # origin points between B1-B3  
        lapply(names(trend_dat), function(x) {
            x_coord = side_arrows[[x]]$x
            x_coord[3] = x_coord[3]+2
            geom_point(
                aes(x = x_coord, y = side_arrows[[x]]$y),
                colour = origin_cols[x])
        })
    
}
