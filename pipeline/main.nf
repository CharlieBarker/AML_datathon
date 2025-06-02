// main.nf
nextflow.enable.dsl = 2

// Define parameters with default values
params.data_input_dir = "${launchDir}/Data"    // Path to the main data directory
params.results_output_dir = "${launchDir}/results" // Where all results will be published
params.project_base_dir = "${launchDir}"       // Base for scripts and other fixed assets


// Process for Feature Selection
process FEATURE_SELECTION {
    tag 'feature_selection'
    publishDir "${params.results_output_dir}", mode: 'copy', pattern: '*.RData'

    input:
    path script: "${params.project_base_dir}/feature_selection.R"
    path data_dir: "${params.data_input_dir}" // Stage the entire Data directory
    path initial_results_dir: "${params.project_base_dir}/Results" // Stage initial Results dir for patient_pathway_activities etc.

    output:
    path "features_data.RData"
    path "partition.RData"
    path "features.RData"

    script:
    """
    # Create the 'Data' and 'Results' symlinks/copies within the task's work dir
    # so R can find them as 'Data/file.RData' and 'Results/file.rds'
    cp -r ${data_dir} ./Data
    cp -r ${initial_results_dir} ./Results

    Rscript ${script}
    """
}


// Workflow definition
workflow {
    // Stage initial data and R scripts
    // All R scripts implicitly use relative paths now.

    // 1. Run Feature Selection
    FEATURE_SELECTION()

}

