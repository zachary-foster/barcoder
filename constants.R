## ---- general_constants ----
program_name = "Barcoder"
program_version = "1.0"


## ---- filtering_parameters ----
max_tm_difference <- 20
minimum_proportion_detected <- 0.0
minimum_individual_per_species <- 1
minimum_individual_per_taxonomic_level <- 1
target_locus = "ALL"

## ---- file_system_constants ----
primersearch_input_path_name <- 'primersearch_input.txt'
primersearch_database_path_name <- 'primersearch_database.fasta'
primersearch_output_path_name <- 'primersearch_output.txt'
amplicon_proportion_path_name <- 'amplicon_proportion_by_clade.txt'
species_proportion_path_name <- 'species_proportion_by_clade.txt'
genus_proportion_path_name <- 'genus_proportion_by_clade.txt'
taxon_proportion_path_name <- 'taxon_proportion_by_level.txt'
taxon_threshold_path_name <- 'distance_threshold_by_level.txt'
clade_threshold_path_name <- 'distance_threshold_by_clade.txt'
taxon_gap_path_name <- 'barcode_gap_by_level.txt'
clade_gap_path_name <- 'barcode_gap_by_clade.txt'


## ---- ggplot_display_constants ----
blank_theme <- theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none", panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank())


## ---- igraph_display_constants ----
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}
add.vertex.shape("fcircle", 
                 plot=mycircle, 
                 parameters=list(vertex.frame.color=1, vertex.frame.width=1))


## ---- generic_display_constants ----
color_pallette <- c("dodgerblue2","#E31A1C", # red
                    "green4",
                    "#6A3D9A", # purple
                    "#FF7F00", # orange
                    "gold1",
                    "skyblue2","#FB9A99", # lt pink
                    "palegreen2",
                    "#CAB2D6", # lt purple
                    "#FDBF6F", # lt orange
                    "gray70", "khaki2",
                    "maroon","orchid1","deeppink1","blue1","steelblue4",
                    "darkturquoise","green1","yellow4","yellow3",
                    "darkorange4","brown")


## ---- format_constants ----
top_clade = "All"
taxonomy_hierarchy <- c("life", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
#taxonomy_hierarchy <- c("life", "domain", "kingdom", "phylum", "class", "order", "family", "genus")
taxonomy_separator <- "|"