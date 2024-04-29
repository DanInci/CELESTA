library(Rmixmod)
library(caret)
library(spdep)
library(ggplot2)
library(reshape2)
library(zeallot)
library(RColorBrewer)

run_celesta_for_image <- function(ExperimentConfig, ExperimentDir, dataset, signature_matrix, image_id, no_labels) {
  print(paste0('Applying CELESTA for image ', image_id))
  ### Prepare imaging data
  dataset <- dataset[dataset$image_name == image_id, ]
  
  spatial_columns <- c('X', 'Y')
  other_columns <- c('image_name', 'cell_id', 'cell_type')
  marker_columns <- setdiff(names(dataset), c(spatial_columns, other_columns))
  imaging_data <- dataset[, c(spatial_columns, marker_columns)]
  
  results <- dataset %>%
    select(image_id = image_name, cell_id, label = cell_type)
  
  image_dir <-  paste0(ExperimentDir, '/', image_id)
  if (!dir.exists(image_dir)){
    dir.create(image_dir, recursive = TRUE)
  }
  
  ### Create CELESTA object.
  CelestaObj <- CreateCelestaObject(project_title = "IMCBenchmark", signature_matrix, imaging_data)
  
  ### Filter out questionable cells. 
  CelestaObj <- FilterCells(CelestaObj, high_marker_threshold=ExperimentConfig$high_marker_threshold, low_marker_threshold=ExperimentConfig$low_marker_threshold)
  
  ### Assign cell types. 
  CelestaObj <- AssignCells(
    CelestaObj, max_iteration=ExperimentConfig$max_iteration, cell_change_threshold=ExperimentConfig$cell_change_threshold, save_result = FALSE
  )
  
  results$pred <- CelestaObj@final_cell_type_assignment[, (CelestaObj@total_rounds+2)]
  
  ### Plot cell type assignment for image
  all_cell_numbers <- 0:no_labels
  palette <- brewer.pal(n = no_labels, name = "Set1")
  
  PlotCellsAnyCombination(cell_type_assignment_to_plot=CelestaObj@final_cell_type_assignment[, (CelestaObj@total_rounds+1)],
                          coords = CelestaObj@coords,
                          prior_info = signature_matrix,
                          cell_number_to_use=all_cell_numbers,
                          cell_type_colors=palette,
                          output_dir = image_dir,
                          save_plot = TRUE)
  
  ### Plot marker expression probability
  PlotExpProb(coords=CelestaObj@coords,
              marker_exp_prob=CelestaObj@marker_exp_prob,
              prior_marker_info = signature_matrix,
              output_dir = image_dir,
              save_plot = TRUE)
  
  return(results)
}
