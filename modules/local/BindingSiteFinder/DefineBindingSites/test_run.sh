# without nextflow
docker run --rm -it -v /Users/melinaklostermann/Documents/projects:/mnt/ melinak/bindingsitefinder:1.0 bin/bash

Rscript /mnt/nf-core-clipseq/modules/local/BindingSiteFinder/DefineBindingSites.R \
	--bw_files_folder "/mnt/nf-core-clipseq/devel_folder/example_inputs/AGO_iCLIP" \
	--peaks "/mnt/nf-core-clipseq/devel_folder/example_inputs/AGO_iCLIP/IP_WT_pureclip_sites.bed" \
	--anno_genes "/mnt/nf-core-clipseq/devel_folder/outputs/gns.rds" \
	--anno_regions "/mnt/nf-core-clipseq/devel_folder/outputs/regions.rds" \
	--output_path "/mnt/nf-core-clipseq/devel_folder/outputs"
