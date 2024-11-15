load("../Raw_data/All_BALA_data.RDA")

all_trend_results<-as.data.frame(matrix(data = NA, ncol = 9, nrow = 0,
                               dimnames = list(NULL, c("measure", "spp_set",
                                                       "comparison","sampling_method", 
                                                       "effsize", "magnitude", 
                                                       "conf.low", "conf.high",
                                                       "p_val"))))

all_trend_results$p_val<-round(all_trend_results$p_val, 3)

source("03_abu_yearly_trends.R") 
source("04_richness_yearly_trends.R") 
source("05_biomass_yearly_trends.R")
source("06_proportions_yearly_trends.R")

save(all_trend_results, file="../Calculated_data/all_trend_results.RDA")
write.xlsx(all_trend_results, file="../Calculated_data/all_trend_results.xlsx")

