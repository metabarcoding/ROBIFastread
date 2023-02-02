# ROBIReadfasta

A preliminary R package allowing to read files produced by OBITools4

## Installation

You first need to install the `devtools` package.

```{r}
install.packages("devtools")
```

Then you can install the `ROBIReadfasta` package itself.

```{r}
devtools::install_git("https://git.metabarcoding.org/obitools/obitools4/robireadfasta.git")
```

## Provided functions

The package provides three main functions:

-   `read_obifasta(file, keys = NULL, verbose = is_robi_verbose())`:

    Reads a FASTA file including annotations inserted in the sequence header by OBITools in JSON format.

-   `extract_features(sequences, ..., verbose = is_robi_verbose())`:

    Extracts some annotations from the JSON, and adds them as supplentary columns to the returned tibble

-   `extract_readcount(sequences, key = "merged_sample")`:

    Extracts from then JSON annotation information about the MOTUs abundances per PCR. The function returns a $PCRs \times MOTUs$ matrix.

## Usage

A simple example relying on a small sample FASTA file provided with the package.

```{r}
library(ROBIFastread)
library(tidyverse)
library(magrittr)

filename <- system.file("extdata", 
                        "sample.fasta", 
                        package="ROBIFastread")

sequences <- read_obifasta(filename,
                           keys = c("forward_match","reverse_match"))

sequences %<>% extract_features("forward_score","reverse_score")

reads <- sequences %>% extract_readcount()
```
