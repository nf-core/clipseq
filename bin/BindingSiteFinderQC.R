args <- commandArgs(trailingOnly = TRUE)

print(args)
bds = args[1]


rmarkdown::render(xxreport_tmp_path,
                  #output_dir = paste0(snakemake@params[[1]], "/results/"),
                  args = args,
                  output_format = "html_document"
)