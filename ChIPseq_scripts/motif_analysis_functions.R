library(dplyr)
library(data.table)
library(tibble)


#################################################################################
## this function open the RSAT matrix scan output and returns the dataframe
process_RSAT_matrix_scan_out = function(file, min_weight = NULL){
  
  command = paste("bash -c", '"grep -v \'^;\'', file, '"', sep = " ")
  

  motifs = fread(input = command, sep = "\t", header = T, stringsAsFactors = F, fill = T, data.table = T) %>% 
    dplyr::filter(ft_type != "limit") %>%
    ## optional filter for weight
    { if(is.null(min_weight)) . else dplyr::filter(., weight >= min_weight) }
  
  
  return(motifs)
  
}



#################################################################################








