# Function to download data from EDI
# Alec Chiono; alec.chiono@colorado.edu

download_data <- function(
  seed = FALSE,
  coloc_veg = FALSE,
  longterm_veg = FALSE
) {
  if (!any(seed, coloc_veg, longterm_veg)) {
    stop("All arguments are FALSE. At least one must be TRUE.")
  }

  library(EDIutils)

  ## [1] seed bank + co-located veg, [2] longterm veg
  package_id <- c("knb-lter-nwt.333.1", "knb-lter-nwt.93.10")

  ## Fetch entity names only for needed packages
  need_pkg1 <- seed || coloc_veg
  need_pkg2 <- longterm_veg

  entity_df <- do.call(
    rbind,
    lapply(package_id[c(need_pkg1, need_pkg2)], function(pid) {
      cbind(package_id = pid, read_data_entity_names(packageId = pid))
    })
  )
  entity_df <- entity_df[
    entity_df$entityName != "saddptqd_ancillary_meta.hh.data.txt",
  ]

  ## Subset pkg1 entities by requested datasets; keep all pkg2 entities
  pkg1_rows <- which(entity_df$package_id == package_id[1])
  pkg2_rows <- which(entity_df$package_id == package_id[2])

  keep_rows <- c(
    pkg1_rows[c(seed, coloc_veg)], # [1] seed, [2] co-located veg
    pkg2_rows # longterm veg
  )
  entity_df <- entity_df[keep_rows, ]

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