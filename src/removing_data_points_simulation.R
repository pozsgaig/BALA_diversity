simulation_df<-NULL
i = 0
options(dplyr.summarise.inform = FALSE)

for (measure in c("Abundance", "Biomass", "Richness"))
{
    
    for (method in unique(count_data$sampling_method)){
        
        for(m in 1:10){
            # m = 1
            
            for (n in rev(seq(50, 100, by= 2)))
                
            {
                # n = 10
                
                if (measure == "Biomass"){
                    biomass_row_num<-ceiling(nrow(biomass_data)*n/100)
                    dat<-biomass_data[sample(1:nrow(biomass_data), biomass_row_num, replace = F),]
                } else
                {count_row_num<-ceiling(nrow(count_data)*n/100)
                dat<-count_data[sample(1:nrow(count_data), count_row_num, replace = F),]}
                
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
                
                
                datlist<-lapply(names(all_dat), function(x){
                    out<-trend_calc(all_dat[[x]], measure = measure, magn_alpha = magn_alpha)
                    out$filter<-x
                    out})
                out_df<-do.call("rbind", datlist)
                out_df$measure<-measure
                out_df$data_retained<-n
                out_df$rep<-m
                
                length(seq(50, 100, by= 2))*10*3*2
                simulation_df<-rbind(simulation_df, out_df)
                i = i+1
                cat("\n\n")
                print(Sys.time())
                cat(round(i/(length(seq(50, 100, by= 2))*10*3*2)*100, 2), "% is done")
                cat("\n\n")
            }   
            
        }
    }
}

save(simulation_df, file = "../Calculated_data/simulation_df.RDA")
