process sortAnnotationForBindingSiteFinder {
	container = 'melinak/bindingsitefinder:latest'
    
    input:
        path gtf_file
    output:
        path "gns.rds"
		path "regions.rds"

	script:
    """
    Rscript /home/mek24iv/nfcore-clipseq/devel_BindingSiteFinder/modules/local/sortAnnotationForBindingSiteFinder/sortAnnotationForBindingSiteFinder.R \\
		$gtf_file \\
		gns.rds \\
		regions.rds
    """
}