## ---- barcoder_functions ----
taxon_edge_list <- function(taxonomy, separator) {
  get_taxon_edge_list <- function(taxon) {
    apply(matrix(c(1:(length(taxon)-1),2:length(taxon)), ncol = 2), 1, function(x) c(taxon[x[1]], taxon[x[2]]))
  }
  taxons <- strsplit(taxonomy, separator, fixed=TRUE)
  taxons <- lapply(taxons, function(x) sapply(seq(1, length(x)), function(y) paste(x[1:y], collapse=separator)))
  edge_list <- t(do.call(cbind,lapply(taxons, FUN=get_taxon_edge_list)))
  edge_list[!duplicated(edge_list),]
}

get_count <- function(data, graph, variable, max_depth) {
  total_count <- 'This_is_an_arbitrary_name'
  data[total_count] <- rep(NA, nrow(data))
  data[data$depth == max_depth, total_count] <- data[data$depth == max_depth, variable]
  for (depth in seq(max_depth - 1, 1)) {
    data[data$depth == depth, total_count] <- 
      apply(data[data$depth == depth, c(variable, 'index')], MARGIN=1, 
            function(x) x[[variable]] + sum(data[neighbors(graph,x[['index']]), total_count]))
  }
  return(data[[total_count]])
}

list_to_DNAbin <- function(sequences, seq_names) {
  seq_names <- seq_names[nchar(sequences) > 0]
  sequences <- sequences[nchar(sequences) > 0]
  x <- unname(sapply(sequences, strsplit, split=''))
  names(x) <- seq_names
  as.DNAbin(x)
}


## ---- generic_functions ----
melting_temperature <- function(primer_sequence) {
  short_calculation <- function(sequence) {
    4 * sum(sequence == 'G' | sequence == 'C') + 2 * sum(sequence == 'A' | sequence == 'T')
  }
  long_calculation <- function(sequence) {
    64.9 + 41 * (sum(sequence == 'G' | sequence == 'C') - 16.4) / length(sequence)
  }
  split_sequence <- strsplit(as.character(primer_sequence), split="")
  sapply(split_sequence, function(x) if (length(x) < 14) short_calculation(x) else long_calculation(x))
}

children <- function(graph, vertex) {
  which(shortest.paths(graph, V(graph)[vertex], mode="out") != Inf)
}

add_alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

latex_format <- function(string) {
  return(gsub('_', '\\\\_', string))
}

list_to_array <- function(the_list) {
  array(unlist(the_list), 
        dim = c(nrow(the_list[[1]]), 
                ncol(the_list[[1]]), 
                length(the_list)), 
        dimnames=list(row.names(the_list[[1]]), 
                      colnames(the_list[[1]]), 
                      names(the_list)))
}

get_edge_parents <-function(graph) {
  get.edges(taxonomy_graph, 1:ecount(taxonomy_graph))[,1]
}

get_edge_children <- function(graph) {
  get.edges(taxonomy_graph, 1:ecount(taxonomy_graph))[,2]
}