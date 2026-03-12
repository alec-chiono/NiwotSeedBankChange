# Script to download data directly from EDI
# Alec Chiono; alec.chiono@colorado.edu

# PACKAGES ---------------------------------------------------------------------
librarian::shelf(EDIutils, readr, dplyr)

# DOWNLOAD DATA ----------------------------------------------------------------
## IDs for the data packages we want to download
package_id <- c("knb-lter-nwt.333.1", "knb-lter-nwt.93.10")

## Get entity IDs and names from data packages
entity_df <- rbind(
  cbind(package_id=package_id[1], read_data_entity_names(packageId=package_id[1])), #seed bank and co-located veg data
  cbind(package_id=package_id[2], read_data_entity_names(packageId=package_id[2])) #longterm veg data
) %>%
  filter(entityName!="saddptqd_ancillary_meta.hh.data.txt") #remove metadata file


## Download data into list
data_list <- lapply(
  1:(nrow(entity_df)), #only get first three files
  function(i) read_csv(file=read_data_entity(packageId=entity_df$package_id[i], entityId=entity_df$entityId[i]), show_col_types=FALSE)
)

## Name each dataset in list
names(data_list) <- entity_df$entityName

## Remove unnecessary, intermediate objects
rm(package_id, entity_df)
