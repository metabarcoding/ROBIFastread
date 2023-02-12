#' @import stringr
#' @import readr
#' @import dplyr
#' @import progress
#' @import Matrix
#' @useDynLib ROBIFastread
#'
#' @author Eric Coissac
NULL

#' @title Read an OBI-formatted FASTQ file
#'
#' @description
#' This function reads an OBI-formatted FASTQ file and returns a tibble containing
#' the sequence information. The FASTQ definition line is expected to contain
#' key-value pairs expressed in the JSON format. If the keys argument is provided,
#' the specified keys will be extracted from the JSON annotations and added as columns to the tibble.
#'
#' @param file a character string specifying the path to the FASTQ file.
#' @param keys an optional character vector of keys to extract from the definition line.
#' @param verbose a logical flag indicating whether to print messages during the processing.
#'
#' @return a tibble with columns "id", "features", "definition", and "sequence".
#' If keys is provided, the additional columns for the specified keys will be added.
#'
#' @examples
#' file <- system.file("extdata", "sample.fasta", package="ROBIFastread")
#' read_obifastq(file)
#'
#' @export
read_obifastq <- function(file,
                          keys = NULL,
                          verbose = is_robi_verbose()) {
  # keys_pattern <- .build_key_value_pattern(keys)

  robimessage("Loading the {file} file in memory...")

  lines <- read_lines(file)

  mfastq <- .Call("R_parse_obiJSONheader",
                  lines[c(TRUE,FALSE,FALSE,FALSE)])
  mfastq <- tibble(id = mfastq[[1]],
                   features = mfastq[[2]],
                   definition = mfastq[[3]],
                   sequence = lines[c(FALSE,TRUE,FALSE,FALSE)],
                   quality = lines[c(FALSE,FALSE,FALSE,TRUE)])

  robimessage("{length(mfasta[[1]])} sequences found")

  if (!is.null(keys)) {
    extract_features(mfastq,keys,verbose)
  } else {
    mfastq
  }
}

