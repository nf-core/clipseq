#!/usr/bin/env Rscript

library(tidyverse)

# Keep backslashes and dollar signs out of this file: Nextflow templates treat them as escapes/placeholders
all_files <- list.files(path = ".")
subtype_files <- all_files[endsWith(all_files, "subtype_premapadjusted.tsv")]
type_files <- all_files[endsWith(all_files, "_type_premapadjusted.tsv")]

files <- list(subtype_files, type_files)

# Subtype-----------------------------------------------------------------------
i = 1

for (file in subtype_files){

  if (i == 1) {
    df <- read.delim(file, header = TRUE) %>%
      select(-Length, -`cDNA..`) %>%
      rename(!!paste0("cDNA_percent_", str_extract(paste(file[1]), "^[^.]+")) := `cDNA...1`)
  }

  else if (i > 1) {

    df_2 <- read.delim(file, header = TRUE) %>%
      select(-Length, -`cDNA..`) %>%
      rename(!!paste0("cDNA_percent_", str_extract(paste(file[1]), "^[^.]+")) := `cDNA...1`)

    df <- full_join(df, df_2, by = "Subtype")

  }

    i <- i + 1

  }

# reformat to mqc requirements
subtype_reformatted <- df %>%
  pivot_longer(starts_with("cDNA"), names_to = "sample", values_to = "value") %>%
  pivot_wider(names_from = "Subtype", values_from = "value")

write_tsv(subtype_reformatted, "summary_subtype_qc.tsv")



# Type -------------------------------------------------------------------------
i = 1

for (file in type_files){

  if (i == 1) {
    df <- read.delim(file, header = TRUE) %>%
      select(-Length, -`cDNA..`) %>%
      rename(!!paste0("cDNA_percent_", str_extract(paste(file[1]), "^[^.]+")) := `cDNA...1`)
  }

  else if (i > 1) {

    df_2 <- read.delim(file, header = TRUE) %>%
      select(-Length, -`cDNA..`) %>%
      rename(!!paste0("cDNA_percent_", str_extract(paste(file[1]), "^[^.]+")) := `cDNA...1`)

    df <- full_join(df, df_2, by = "Type")

  }

  i <- i + 1

}

subtype_reformatted <- df %>%
  pivot_longer(starts_with("cDNA"), names_to = "sample", values_to = "value") %>%
  pivot_wider(names_from = "Type", values_from = "value")

write_tsv(subtype_reformatted, "summary_type_qc.tsv")

writeLines(c(
  '"${task.process}":',
  paste0("    r-tidyverse: ", as.character(packageVersion("tidyverse")))
), "versions.yml")
