download_data <- function(seed = TRUE, veg = TRUE) {
  if (seed == FALSE & veg == FALSE) {
    stop(
      "seed and veg arguments both FALSE. One must be true to download data."
    )
  }

  library(EDIutils)

  ## IDs for the data packages we want to download
  ## [1] seed bank + co-located veg, [2] longterm veg
  package_id <- c("knb-lter-nwt.333.1", "knb-lter-nwt.93.10")

  ## Select packages based on arguments
  selected_ids <- package_id[c(seed, veg)]

  ## Get entity IDs and names from selected data packages
  entity_df <- do.call(
    rbind,
    lapply(selected_ids, function(pid) {
      cbind(package_id = pid, read_data_entity_names(packageId = pid))
    })
  )
  entity_df <- entity_df[
    entity_df$entityName != "saddptqd_ancillary_meta.hh.data.txt",
  ]

  ## Download data into named list
  data_list <- lapply(
    seq_len(nrow(entity_df)),
    function(i) {
      con <- textConnection(
        rawToChar(read_data_entity(
          packageId = entity_df$package_id[i],
          entityId = entity_df$entityId[i]
        ))
      )
      on.exit(close(con))
      read.csv(con)
    }
  )
  names(data_list) <- entity_df$entityName

  return(data_list)
}
