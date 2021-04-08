#' check dependencies and install missing packages, then load them
#' @description Checks which packages are install and which are missing. Missing packages, 
#' including there dependencies will be isntalled. 
#' \nAn exception install routine is created for Rgeostats, which is not avialable on 
#' CRAN or GIT. 
#' \n Once all packages are available on the local machine, all required packages 
#' will be loaded.
#' @param list.of.packages list, list containing all required packages
#' @examples 
#' lop <- list('ggplot2','scales')
#' check_dep(lop)
#' @export
#' 
check_dep <- function(list.of.packages){
  #check for rgeostats as it isn't on git or cran...which makes life a bit harder...
  if ('RGeostats' %in% list.of.packages){
    if("rgeostats" %in% row.names(installed.packages())==FALSE){
      os = Sys.info()['sysname']
      if(os == 'Windows'){
        zip_URL <- "http://cg.ensmp.fr/rgeos/DOWNLOAD/RGeostats_11.2.12.zip"
        zip_filename <- file.path(tempdir(), "rgeostats.zip")
      }else if(os == 'Linux'){
        zip_URL <- "http://cg.ensmp.fr/rgeos/DOWNLOAD/RGeostats_11.2.12_linux64.tar.gz"
        zip_filename <- file.path(tempdir(), "rgeostats.tar.gz")
      }else{
        zip_URL <- "http://cg.ensmp.fr/rgeos/DOWNLOAD/RGeostats_11.2.12_macosx.tgz"
        zip_filename <- file.path(tempdir(), "rgeostats.tgz")
      }
      download.file(zip_URL, destfile = zip_filename, mode = "wb")
      install.packages(pkgs = zip_filename, repos = NULL)
      unlink(zip_filename)
    }
    require('RGeostats')
    list.of.packages = list.of.packages[list.of.packages != 'RGeostats']
  }
  #check which pakacges are installed
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  #install missing ones, with all dependencies
  if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)
  #load all libraries
  tmp = lapply(list.of.packages, require, character.only = TRUE)
  #clean the tmp output
  rm(tmp)
}
