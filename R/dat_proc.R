library(tidyverse)
library(tbeptools)

# transect data -----------------------------------------------------------

# import entire transect dataset as JSON
transect <- read_transect(training = FALSE) %>% 
  select(-Crew, -MonitoringAgency)

# get transect species occurrence summaries
transectocc <- anlz_transectocc(transect)

save(transect, file = 'data/transect.RData', compress = 'xz')
save(transectocc, file = 'data/transectocc.RData', compress = 'xz')