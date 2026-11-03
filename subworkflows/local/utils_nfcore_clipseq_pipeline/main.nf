// CLIPseq initialization and reporting around the shared nf-core utilities.
include { UTILS_NEXTFLOW_PIPELINE; getWorkflowVersion } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFCORE_PIPELINE; completionEmail; completionSummary } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFSCHEMA_PLUGIN } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap } from 'plugin/nf-schema'

workflow PIPELINE_INITIALISATION {
    take:
    nextflow_cli_args

    main:
    UTILS_NEXTFLOW_PIPELINE(
        params.version,
        false,
        params.outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    def before_text = "${workflow.manifest.name} ${getWorkflowVersion()}\n"
    def after_text = """If you use ${workflow.manifest.name} for your analysis please cite:

* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md
"""
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        params.validate_params,
        null,
        params.help,
        params.help_full,
        params.show_hidden || params.show_hidden_params,
        before_text,
        after_text,
        "nextflow run ${workflow.manifest.name} -profile <docker/singularity/apptainer/institute> --input samplesheet.csv --outdir <OUTDIR>",
        null
    )
    UTILS_NFCORE_PIPELINE(nextflow_cli_args)
    validateInputParameters()

    emit:
    dummy_emit = true
}

workflow PIPELINE_COMPLETION {
    take:
    email
    email_on_fail
    plaintext_email
    outdir
    monochrome_logs
    hook_url
    multiqc_report

    main:
    def summary_params = paramsSummaryMap(workflow)
    def multiqc_reports = multiqc_report.toList()

    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : []
            )
        }
        completionSummary(monochrome_logs)
        if (hook_url) {
            completionNotification(summary_params, hook_url)
        }
    }
}

// Keep pipeline-specific checks here; samplesheet content is checked by INPUT_CHECK.
def validateInputParameters() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        error("Genome '${params.genome}' not found in the configured genomes. Available keys: ${params.genomes.keySet().join(', ')}")
    }
    if (!params.input) {
        error("Please provide an input samplesheet using --input samplesheet.csv")
    }
    if (workflow.profile.contains('awsbatch')) {
        if (!params.awsqueue || !params.awsregion) {
            error('Specify --awsqueue and --awsregion when running on AWS Batch.')
        }
        if (!params.outdir.startsWith('s3:')) {
            error('Specify an S3 bucket for --outdir when running on AWS Batch.')
        }
    }
}

def methodsDescriptionText(run_workflow, mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = run_workflow.toMap()
    meta["manifest_map"] = run_workflow.manifest.toMap()

    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}


// Preserve the existing optional Teams/Slack notification payloads.
def completionNotification(summary_params, hook_url) {

    def summary = [:]
    for (group in summary_params.keySet()) {
        summary << summary_params[group]
    }

    def misc_fields = [:]
    misc_fields['start']                                = workflow.start
    misc_fields['complete']                             = workflow.complete
    misc_fields['scriptfile']                           = workflow.scriptFile
    misc_fields['scriptid']                             = workflow.scriptId
    if (workflow.repository) misc_fields['repository']  = workflow.repository
    if (workflow.commitId)   misc_fields['commitid']    = workflow.commitId
    if (workflow.revision)   misc_fields['revision']    = workflow.revision
    misc_fields['nxf_version']                          = workflow.nextflow.version
    misc_fields['nxf_build']                            = workflow.nextflow.build
    misc_fields['nxf_timestamp']                        = workflow.nextflow.timestamp

    def msg_fields = [:]
    msg_fields['version']      = getWorkflowVersion()
    msg_fields['runName']      = workflow.runName
    msg_fields['success']      = workflow.success
    msg_fields['dateComplete'] = workflow.complete
    msg_fields['duration']     = workflow.duration
    msg_fields['exitStatus']   = workflow.exitStatus
    msg_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    msg_fields['errorReport']  = (workflow.errorReport ?: 'None')
    msg_fields['commandLine']  = workflow.commandLine.replaceFirst(/ +--hook_url +[^ ]+/, "")
    msg_fields['projectDir']   = workflow.projectDir
    msg_fields['summary']      = summary << misc_fields

    // Render the JSON template
    def engine       = new groovy.text.GStringTemplateEngine()
    // Different JSON depending on the service provider
    // Defaults to "Adaptive Cards" (https://adaptivecards.io), except Slack which has its own format
    def json_path     = hook_url.contains("hooks.slack.com") ? "slackreport.json" : "adaptivecard.json"
    def hf            = new File("${workflow.projectDir}/assets/${json_path}")
    def json_template = engine.createTemplate(hf).make(msg_fields)
    def json_message  = json_template.toString()

    // POST
    def post = new URL(hook_url).openConnection();
    post.setRequestMethod("POST")
    post.setDoOutput(true)
    post.setRequestProperty("Content-Type", "application/json")
    post.getOutputStream().write(json_message.getBytes("UTF-8"));
    def postRC = post.getResponseCode();
    if (! postRC.equals(200)) {
        log.warn(post.getErrorStream().getText());
    }
}
