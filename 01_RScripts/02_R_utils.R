

#' Parse config YAML and add path prefix
#'
#' @param conf a YAML config file parsed as list
#' @param dir starting directory. Default: .
#'
#' @return A list where all the `files` attributes are prefixed with `dir` path
#' @export
#'
#' @examples
prefix_config_paths <- function(conf, dir = "."){
  if(!is.null(conf[["dir"]])){
    ## add prefix
    dir <- paste(dir, "/", conf[["dir"]], sep = "")
  }
  
  ## append dir to files
  if(!is.null(conf[["files"]])){
    conf[["files"]] <- purrr::map(
      .x = conf[["files"]],
      .f = ~ paste(dir, "/", .x, sep = "")
    )
  }
  
  ## recursion on child nodes
  conf <- purrr::map_if(
    .x = conf,
    .p = purrr::map_lgl(conf, is.list) & names(conf) != "files",
    .f = ~prefix_config_paths(.x, dir)
  )
  
  ## update dir for next child nodes
  conf[["dir"]] <- dir
  
  return(conf)
}
