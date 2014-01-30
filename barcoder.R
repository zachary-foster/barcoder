## ---- imports ----
library(png)
library(grid)
library(igraph)
library(ggplot2)
library(spider)
library(tools)
library(phyloch)
library(Hmisc)
library(knitr)
library(xtable)
library(plyr)
library(plotrix)

source("functions.R")
source("constants.R")


## ---- argument_parsing ----
#arguments <- commandArgs(TRUE)
arguments <- c("/home/local/USDA-ARS/fosterz/ITS_analysis/PR2_database/refseq_gb-191_head1000.fasta",
               "/home/local/USDA-ARS/fosterz/ITS_analysis/SSU_primers.fasta",
               "/home/local/USDA-ARS/fosterz/ITS_analysis",
               "0")
names(arguments) <- c("database", "primers", "output", "mismatch")
output_path <- file.path(arguments["output"], basename(file_path_sans_ext(arguments["primers"])))
unlink(output_path, recursive=TRUE)
dir.create(output_path)


## ---- import_reference_database ----
#read sequence database
database <- read.dna(arguments["database"], format="fasta")

#make sequence data frame
sequences <- sapply(as.character(database), paste, collapse='')
split_name <- strsplit(c(labels(sequences)), taxonomy_separator, fixed = TRUE)
sequence_data <- data.frame(row.names=sapply(split_name, function(x) x[1]), 
                            sequence=as.character(unname(sequences)))
sequence_data$taxonomy <- sapply(split_name, function(x) paste(c(top_clade, x[2:length(x)]), collapse=taxonomy_separator))
sequence_data$taxonomy_tips <- unlist(lapply(strsplit(sequence_data$taxonomy, taxonomy_separator, fixed=TRUE), function(x) x[length(x)]))
sequence_data$taxon_depth <- sapply(sequence_data$taxonomy, length)
sequence_data$length <- nchar(as.character(sequence_data$sequence))
sequence_data <- cbind(sequence_data, matrix(unlist(strsplit(sequence_data$taxonomy, taxonomy_separator, fixed = TRUE)), 
                                             ncol=length(taxonomy_hierarchy), 
                                             byrow=TRUE, 
                                             dimnames=list(NULL, taxonomy_hierarchy)))
rm(sequences)
rm(split_name)
rm(database)

#construct taxonomy graph
taxonomy_graph <- graph.edgelist(taxon_edge_list(sequence_data$taxonomy, taxonomy_separator))
taxon_root <- V(taxonomy_graph)[1]

#make taxon data frame (index is veticies of the taxonomy graph)
taxon_data <- data.frame(row.names=V(taxonomy_graph)$name)
taxon_data$index <- c(V(taxonomy_graph))

#store the distance of all taxons from the root
taxon_data$depth <- sapply(get.shortest.paths(taxonomy_graph, from=taxon_root), length)
max_taxon_depth <- max(taxon_data$depth)

#store number of database sequences for each taxon defined at that level explicitly
taxon_data$sequence_count <- tapply(sequence_data$taxonomy, list(sequence_data$taxonomy), length)[rownames(taxon_data)]
taxon_data$sequence_count[is.na(taxon_data$sequence_count)] <- 0

#calculate the number of database sequences at each taxon as the sum of their subgroups
taxon_data$total_count <- get_count(taxon_data, taxonomy_graph, 'sequence_count', max_taxon_depth)

#set default graphing display parameters
taxonomy_graph_layout <- layout.reingold.tilford(taxonomy_graph, root = 1, circular = TRUE)
V(taxonomy_graph)$size <- (log(taxon_data$total_count + .5) / max(log(taxon_data$total_count) + .5)) * 10
V(taxonomy_graph)$label.cex <- V(taxonomy_graph)$size * .05 + .15
V(taxonomy_graph)$label.color <- "black"
V(taxonomy_graph)$color[order(V(taxonomy_graph)$name)] <- rainbow_hcl(length(V(taxonomy_graph)$name))
V(taxonomy_graph)$alpha <- (max_taxon_depth*1.5 - taxon_data$depth) / (max_taxon_depth*1.5)
V(taxonomy_graph)$color <- mapply(add_alpha, V(taxonomy_graph)$color, V(taxonomy_graph)$alpha)
E(taxonomy_graph)$depth <- taxon_data$depth[get_edge_parents(taxonomy_graph)]
max_edge_depth <- max(E(taxonomy_graph)$depth)
E(taxonomy_graph)$width <- V(taxonomy_graph)$size[get_edge_children(taxonomy_graph)] * .4
E(taxonomy_graph)$color <- sapply(((max_edge_depth*4 - E(taxonomy_graph)$depth) / (max_edge_depth*4)) * .3,
                                  function(x) rgb(red=.3,green=.3,blue=.3,alpha=x))


## ---- primersearch ----

#Save modified database file for primersearch
primersearch_database_path <- file.path(output_path, primersearch_database_path_name)
sequences <- sequence_data$sequence
names(sequences) <- apply(matrix(c(rownames(sequence_data), as.character(sequence_data$taxonomy)), ncol = 2),
                          MARGIN=1, FUN=paste, collapse = ' ')
write.dna(sequences, 
          primersearch_database_path, 
          format = "fasta", 
          nbcol = -1, 
          colw =  max(nchar(as.character(sequences))))
rm(sequences)

#load input primers from fasta file
primer_data <- data.frame(sequence=sapply(as.character(read.dna(arguments["primers"], format="fasta")), paste, collapse=''))

#Calculate primer Tm
primer_data$tm <- melting_temperature(primer_data$sequence)

#Create all unique combinations of primers 
primer_pair_data <- data.frame(row.names=apply(t(combn(row.names(primer_data), 2)), 1, paste, collapse='__'), 
                               t(combn(row.names(primer_data), 2)),
                               t(combn(as.character(primer_data$sequence), 2)))
names(primer_pair_data) <- c("name_1", "name_2", "sequence_1", "sequence_2")

#Filter out primer combinations with too large a difference in Tm
primer_pair_data$tm_difference <- abs(primer_data$tm[primer_pair_data$name_1] - primer_data$tm[primer_pair_data$name_2])
primer_pair_data <- droplevels(primer_pair_data[primer_pair_data$tm_difference <= max_tm_difference, ])

#write primersearch input file
primersearch_input_path <- file.path(output_path, primersearch_input_path_name)
write.table(primer_pair_data[,c("sequence_1", "sequence_2")], 
            primersearch_input_path, 
            sep="\t",
            row.names = TRUE,
            col.names = FALSE,
            quote = FALSE) 


#primersearch execution
primersearch_output_path <- file.path(output_path, primersearch_output_path_name)
primersearch_arguments <- unname(c(arguments["database"], primersearch_input_path, arguments["mismatch"], primersearch_output_path))
primersearch_command_line <- gettextf('primersearch -seqall %s -infile %s -mismatchpercent %s -outfile %s',
                                      primersearch_database_path,
                                      primersearch_input_path,
                                      arguments["mismatch"],
                                      primersearch_output_path)
system(primersearch_command_line)


#primersearch output parsing 
primersearch_output <- readLines(primersearch_output_path)
primer_indexes <- grep("Primer name ", primersearch_output, perl=TRUE, value=FALSE)
primer_values <- grep("Primer name ", primersearch_output, perl=TRUE, value=TRUE)
primer_values <- unlist(lapply(strsplit(primer_values, " ", fixed=TRUE), function(x) x[3]))
id_indexes <- grep("\tSequence: ", primersearch_output, perl=TRUE, value=FALSE)
id_values <- grep("\tSequence: ", primersearch_output, perl=TRUE, value=TRUE)
id_values <- sapply(strsplit(id_values, " ", fixed=TRUE), function(x) x[2])

primersearch_data <- data.frame(id_values=id_values)

primersearch_data$taxon_values <- as.factor(sapply(primersearch_output[id_indexes + 1],
                                                   function(x) substr(x, 2, nchar(x))))
primersearch_data$forward_primer <- as.factor(sapply(strsplit(primersearch_output[id_indexes + 2], " ", fixed=TRUE), 
                                                     function(x) substr(x[1], 2, nchar(x))))
primersearch_data$reverse_primer <- as.factor(sapply(strsplit(primersearch_output[id_indexes + 3], " ", fixed=TRUE), 
                                                     function(x) substr(x[1], 2, nchar(x))))
primersearch_data$forward_location <- as.numeric(sapply(strsplit(primersearch_output[id_indexes + 2], " ", fixed=TRUE), 
                                                        function(x) x[6]))
primersearch_data$reverse_location <- as.numeric(sapply(strsplit(primersearch_output[id_indexes + 3], " ", fixed=TRUE), 
                                                        function(x) substr(x[6], 2, nchar(x[6]) - 1)))
primersearch_data$forward_mismatch <- as.numeric(sapply(strsplit(primersearch_output[id_indexes + 2], " ", fixed=TRUE), 
                                                        function(x) x[8]))
primersearch_data$reverse_mismatch <- as.numeric(sapply(strsplit(primersearch_output[id_indexes + 3], " ", fixed=TRUE), 
                                                        function(x) x[8]))              
primersearch_data$amplicon_length <- as.numeric(sapply(strsplit(primersearch_output[id_indexes + 4], " ", fixed=TRUE), 
                                                       function(x) x[3]))              
primersearch_data$primer_pair <- as.factor(rep(primer_values, 
                                               as.numeric(table(cut(id_indexes, c(primer_indexes, length(primersearch_output)))))))
primersearch_data <- cbind(primersearch_data, sequence_data[as.character(primersearch_data$id_values),taxonomy_hierarchy], row.names=NULL)


## ---- sensitivity_statistics ----  

#calculate proportion of sequences detected by each primer pair
primer_pair_data$amplicon_count <- 0
counts <- c(by(primersearch_data, primersearch_data$primer_pair, nrow))
primer_pair_data[labels(counts), c("amplicon_count")] <- counts
primer_pair_data$proportion_detected <- primer_pair_data$amplicon_count/nrow(sequence_data)
rm(counts)

#Filter out primer pairs that have too little amplification
passing_pairs <- row.names(primer_pair_data)[which(primer_pair_data$proportion_detected > minimum_proportion_detected)]
primer_pair_data <- droplevels(primer_pair_data[passing_pairs,])
primersearch_data <- droplevels(primersearch_data[which(primersearch_data$primer_pair %in% passing_pairs), ])
row.names(primersearch_data) = NULL

#store number of amplicons for each taxon defined at that level explicitly
taxon_data$amplicon_count <- 0
unordered_count <- tapply(primersearch_data$taxon_values, 
                          INDEX=list(taxon=primersearch_data$taxon_values), 
                          FUN=length)
taxon_data$amplicon_count <- unordered_count[rownames(taxon_data)]
taxon_data$amplicon_count[is.na(taxon_data$amplicon_count)] <- 0

#calculate the number of amplicons at each taxon as the sum of their subgroups
taxon_data$total_amplicon_count <- get_count(taxon_data, taxonomy_graph, 'amplicon_count', max_taxon_depth)

#Calculate the number of amplicons at each taxon as for each primer set
get_amplicon_count <- function(amplicons, data, graph, max_depth) {
  variable <- 'This_is_arbitrary'
  data[variable] <- 0
  unordered_count <- tapply(amplicons, INDEX=list(taxon=amplicons), FUN=length)
  data[variable] <- unordered_count[rownames(data)]
  data[is.na(data[variable]), variable] <- 0
  return(get_count(data, graph, variable, max_depth))
}
primer_counts <- data.frame(c(by(data=primersearch_data$taxon_values, 
                                 INDICES=list(primer_pair=primersearch_data$primer_pair), 
                                 FUN=function(x) get_amplicon_count(x, taxon_data, taxonomy_graph, max_taxon_depth))))
rownames(primer_counts) <- rownames(taxon_data)

#Calculate the proportion of taxons represented by amplicons 
primer_proportions <- data.frame((apply(primer_counts, MARGIN=2, function(x) x / taxon_data$total_count)))
primer_proportions <- primer_proportions[,rev(order(primer_proportions[1,]))]
if (max(primer_proportions) > 1) {
  warning("Taxon proportions greater than 1 detected. Changing to 1. ")
  primer_proportions[primer_proportions > 1] <- 1
}

### Calculate proportion of taxonomic groups matched for each taxonomic group for each primer pair ###

taxon_data$depth <- factor(taxon_data$depth)
levels(taxon_data$depth) <- taxonomy_hierarchy
taxonomic_rank_data <- data.frame(row.names=taxonomy_hierarchy)
taxonomic_rank_data$count <- c(table(taxon_data$depth))
get_taxonomic_rank_count <- function(proportions, depth, level=0) {
  apply(proportions, MARGIN=2, FUN=tapply, INDEX=depth, function(x) sum(x > level))
}
level_range <- seq(0, .9, by=0.01)
taxonomic_rank_counts <- lapply(level_range,
                                function(x) get_taxonomic_rank_count(primer_proportions, taxon_data$depth, x))
names(taxonomic_rank_counts) <- level_range
taxonomic_rank_counts <- list_to_array(taxonomic_rank_counts)
taxonomic_rank_proportions <- apply(taxonomic_rank_counts, MARGIN=c(2,3), function(x) x / taxonomic_rank_data$count)


## ---- sequence_filtering ----  

primersearch_data <- do.call(rbind, lapply(split(primersearch_data, list(primersearch_data$primer_pair, primersearch_data$taxon_values)),
                                           function(x) if (nrow(x) >= minimum_species_count) x))


filter <- ddply(primersearch_data, .(taxon_values, primer_pair), nrow)
filter <- filter[filter$V1 >= minimum_individual_per_species,]
filter <- (primersearch_data$taxon_values %in% filter$taxon_values) & (primersearch_data$primer_pair %in% filter$primer_pair)
primersearch_data <- primersearch_data[filter, ]


## ---- distance_analysis ----  

#extract amplicon sequences
primersearch_data$amplicon <- sequence_data$sequence[primersearch_data$id_values]
starts <- primersearch_data$forward_location
ends <- starts + primersearch_data$amplicon_length
primersearch_data$amplicon  <- substr(as.character(primersearch_data$amplicon), starts, ends)

primersearch_data <- primersearch_data[nchar(primersearch_data$amplicon) > 0,]



unaligned_sequences <- lapply(split(primersearch_data, list(primersearch_data$primer_pair, primersearch_data$phylum)),
                              function(x) if (nrow(x) > 0) {list_to_DNAbin(x$amplicon, paste(as.character(x$taxon_values), 
                                                                                             as.character(x$id_values), 
                                                                                             sep=taxonomy_separator))})
unaligned_sequences <- unaligned_sequences[sapply(unaligned_sequences, function(x) length(x) >= minimum_individual_per_taxonomic_level)]


#align amplicon sequences
alignments <- lapply(unaligned_sequences[1], mafft, method="retree 1")

#calculate distance matrix
distances <- lapply(alignments, dist.dna, pairwise.deletion = TRUE)