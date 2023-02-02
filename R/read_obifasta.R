#' @import stringr
#' @import readr
#' @import dplyr
#' @import progress
#' @import Matrix
#' @useDynLib ROBIFastread
#'
#' @author Eric Coissac
NULL

#' @title Read an OBI-formatted FASTA file
#'
#' @description
#' This function reads an OBI-formatted FASTA file and returns a tibble containing
#' the sequence information. The FASTA definition line is expected to contain
#' key-value pairs expressed in the JSON format. If the keys argument is provided,
#' the specified keys will be extracted from the JSON annotations and added as columns to the tibble.
#'
#' @param file a character string specifying the path to the FASTA file.
#' @param keys an optional character vector of keys to extract from the definition line.
#' @param verbose a logical flag indicating whether to print messages during the processing.
#'
#' @return a tibble with columns "id", "features", "definition", and "sequence".
#' If keys is provided, the additional columns for the specified keys will be added.
#'
#' @examples
#' file <- system.file("extdata", "sample.fasta", package="ROBIFastread")
#' read_obifasta(file)
#'
#' @export
read_obifasta <- function(file,
                          keys = NULL,
                          verbose = is_robi_verbose()) {
  # keys_pattern <- .build_key_value_pattern(keys)

  robimessage("Loading the {file} file in memory...")

  mfasta <- read_file_raw(file)

  robimessage("Looking for the beginning of the sequences...")

  mfasta <- .Call("R_parse_fasta", mfasta)
  names(mfasta) <- c("id", "features", "definition", "sequence")

  mfasta <- as_tibble(mfasta)

  robimessage("{length(mfasta[[1]])} sequences found")

  if (!is.null(keys)) {
    extract_features(mfasta,keys,verbose)
  } else {
    mfasta
  }
}

#' Extract Features from JSON Strings
#'
#' This function extracts features from JSON strings within a tibble and binds the features to the tibble as additional columns.
#'
#' @param sequences A tibble with a column named "features" containing JSON strings.
#' @param ... A character vector or several characters values specifying the keys to extract from the JSON strings.
#' @param verbose A logical indicating whether to show the progress bar.
#' @return A tibble with the extracted features bound as additional columns.
#' @export
#' @examples
#' \dontrun{
#' sequences <- tibble(
#' features = c('{"a": 1, "b": 2}', '{"a": 3, "b": 4}')
#' )
#' extract_features(sequences, c("a", "b"))
#' }
extract_features <- function(sequences,
                             ...,
                             verbose = is_robi_verbose()) {
  keys = unlist(list(...))
  sequences <- as_tibble(sequences)
  pb <- progress_bar$new(total = length(sequences[[1]]))

  features <- lapply(
    sequences$features,
    function(x) {
      pb$tick()
      val <- list()
      json <- suppressWarnings(jsonify::from_json(x))
      data <- json[keys,drop=FALSE]
      names(data) <- keys
      data
    }
  )
  features <- do.call(bind_rows, features)
  colnames(features) <- keys
  sequences <- bind_cols(sequences,features)
  sequences
}

#' Extract Read Counts
#'
#' This function extracts the read count data from the input sequences.
#'
#' @param sequences A data frame containing the sequences information, as provided by the read_obifasta function.
#' @param key A character string specifying the key name in the JSON data.
#'
#' @return A sparse matrix representing the read counts for each sample and each MOTU.
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "sample.fasta", package="ROBIFastread")
#' seq <- read_obifasta(file)
#' m <- extract_readcount(seq)
#'
extract_readcount <- function(sequences, key = "merged_sample") {
  pb <- progress_bar$new(total = length(sequences[[1]]))
  counts <- lapply(
    seq_along(sequences$features),
    function(x) {
      pb$tick()
      val <- list()
      json <- suppressWarnings(jsonify::from_json(sequences$features[[x]]))
      data <- json[[key]]
      tibble::tibble(
        motu = sequences$id[[x]],
        sample = names(data),
        count = unlist(data)
      )
    }
  )
  data  <- do.call(bind_rows,counts)
  data$motu <- factor(data$motu,levels = sequences$id)
  samples <- unique(data$sample)
  data$sample <- factor(data$sample,levels = samples)
  sparseMatrix(i = as.integer(data$sample),
              j = as.integer(data$motu),
              x = data$count,
              dims = c(length(samples),nrow(sequences)),
              dimnames= list(samples,sequences$id))

}
