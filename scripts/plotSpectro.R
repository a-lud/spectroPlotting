## Install necessary packages if not already
if(! 'tidyverse' %in% rownames(installed.packages())){install.packages('tidyverse')}
if(! 'magrittr' %in% rownames(installed.packages())){install.packages('magrittr')}
if(! 'devtools' %in% rownames(installed.packages())){install.packages('devtools')}
if(! 'ggpomological' %in% rownames(installed.packages())){devtools::install_github('gadenbuie/ggpomological')}

## Load packages
library(tidyverse)
library(magrittr)
library(assertthat)

## Function to plot spectro data
plotSpectro <- function(vec_samples = NULL, path_data, ext_files = NULL, path_phenotypeCSV,
                        n_commentLines = 14, min_reflectance = 0, max_reflectance = 120, subset_phenotype,
                        save_plot = FALSE, path_plot = NULL, name_plot = NULL, image_type = 'pdf'){
  
  ## Checking inputs - error if conditions aren't met
  assert_that(dir.exists(path_data), msg = "Directory path to spectro data can't be resolved. Check the path exists")
  assert_that(!is.null(ext_files), msg = "File extension hasn't been provided")
  assert_that(file.exists(path_phenotypeCSV), msg = 'No phoenotype data has been provided')
  if(isTRUE(save_plot)){
    assert_that(!is.null(path_plot), msg = "Path to plot directory has not been set")
    assert_that(!is.null(name_plot), msg = "Name of plot file has not been set")
  }
  
  ## Read data
  rawSpectro <- list.files(path = path_data, pattern = ext_files, full.names = TRUE) %>%
    set_names(sub(ext_files, '', basename(.)))
  
  ## Import phoenotype data
  pheno <- read_csv(file = path_phenotypeCSV, col_names = c('id', 'phenotype'), col_types = cols())
  
  ## Subset samples for selection
  if(!is.null(vec_samples)){
    rawSpectro <- rawSpectro[str_detect(string = rawSpectro, pattern = vec_samples)]
  }
  
  ## Read in data
  rawSpectro <- rawSpectro %>%
    map(read_tsv,
        col_names = c("wavelength", "reflectance"),
        col_types = cols(),
        skip = n_commentLines) %>%
    bind_rows(.id = 'sample')
  
  ## Get mean reflectance for each wavelength value grouped by colour + sample
  meanSpectro <- rawSpectro %>%
    separate(col = sample, into = c("id", "colour", "replicate"), extra = 'drop') %>%
    group_by(id, colour, wavelength) %>%
    summarise(mean_reflectance = mean(reflectance))
  
  ## Join with phenotype information
  meanSpectro <- left_join(meanSpectro, pheno, by = 'id')
  
  ## Subset for phenotype
  if(! is.null(subset_phenotype)){
    meanSpectro <- filter(.data = meanSpectro, phenotype %in% subset_phenotype)
  }
  
  ## Plotting data
  p <- meanSpectro %>%
    filter(mean_reflectance > min_reflectance & mean_reflectance <= max_reflectance) %>%
    ggplot(aes(x = wavelength, y = mean_reflectance, colour = colour)) +
    geom_point() +
    geom_smooth() +
    labs(x = 'Wavelength',
         y = 'Reflectance (mean)') +
    scale_colour_manual(values = ggpomological:::pomological_palette) +
    theme_bw() +
    facet_wrap(~id, dir = "v", ncol = 4, scales = "free") +
    xlim(200, 900)
  
  # ## Save plot or just return it
  if(isTRUE(save_plot)){
    
    ## Create output directory if it doesn't exist already
    dir.create(path = normalizePath(path_plot), showWarnings = FALSE, recursive = TRUE)
    
    ## Build file name
    n <- paste(name_plot, image_type, sep = ".")
    
    ggsave(filename = n,
           device = image_type,
           path = normalizePath(path_plot),
           plot = p)
  } else {
    return(p)
  }
}

## Updated - 2020/03/11
plotSpectro(path_data = 'data/jenna_testing', 
            ext_files = '.txt', 
            n_commentLines = 14, 
            path_phenotypeCSV = 'data/jenna_sample.csv', 
            subset_phenotype = c('Stripes', 'Banded'))

## Using function to plot the data
plotSpectro(vec_samples = c("KLS0001", "KLS0002"),
            path_data = "data/testing",
            ext_files = '.txt' ,
            save_plot = TRUE,
            path_plot = "data/plots/juv",
            name_plot = "JUV_2Samples", 
            image_type = "png")

## Plotting juv vs adult
plotSpectro(vec_samples = c("KLS0001", "KLS1030"),
            path_data = "data/testing",
            ext_files = '.txt' ,
            save_plot = TRUE,
            path_plot = "data/plots/juv_adult",
            name_plot = "JUV-ADULT-comparison", 
            image_type = "pdf")
