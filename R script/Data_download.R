library(rgbif)
library(dplyr)
library(purrr)
library(readr)  
library(magrittr) 
library(rgbif) 
library(taxize) 


### Downloading data
### Following example get data from https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/

## access data on gbif
user <- "xxx" # your gbif.org username
pwd <- "xxx" # your gbif.org password
email <- "xxx" # your email


#match names with taxon jeys in gbif
#taxon keys do not capture the 790 species
#I used the online tool, and matched the species with success
gbif_taxon_keys <-
  readr::read_csv('Data\\taxonomy_NA.csv') %>%
  pull("scientific_name") %>% # use fewer names if you want to just test
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
  filter(kingdom == "Animalia") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys


# download all data for species
# occ_download(
#   pred_in("taxonKey", gbif_taxon_keys),
#   format = "SIMPLE_CSV",
#   user=user,pwd=pwd,email=email
# )
# Download key: 0323299-200613084148143
