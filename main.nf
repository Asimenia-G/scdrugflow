nextflow.enable.dsl = 2

params.input              = null
params.outdir             = 'results'
params.mode               = 'inspect'  // inspect | analyse | ra
params.condition_column   = 'perturbation'
params.control_label      = 'control'
params.matrix             = 'X'        // X or a layer name
params.input_kind         = 'counts'   // counts | log1p
params.subset_column      = ''
params.subset_value       = ''
params.min_genes          = 100
params.min_counts         = 300
params.min_cells_per_gene = 3
params.min_group_cells    = 20

// Gehring-Pachter example: retinoic acid contrast with other doses fixed.
params.ra_bmp4            = 8
params.ra_egf_bfgf        = 200
params.ra_mix             = 0
params.ra_dose            = 10
params.ra_batches         = '0,1'

process INSPECT {
    tag { dataset.simpleName }
    conda "${projectDir}/environment.yml"
    publishDir params.outdir, mode: 'copy', saveAs: { name ->
        name == 'inspection.json' ? 'inspection.json' : null
    }

    input:
    path dataset

    output:
    tuple path(dataset), path('inspection.json')

    script:
    """
    python ${projectDir}/scripts/inspect_dataset.py \\
        --input '${dataset}' \\
        --output inspection.json
    """
}

process ANALYSE {
    tag { dataset.simpleName }
    conda "${projectDir}/environment.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(dataset), path(inspection)

    output:
    path 'analysis'

    script:
    """
    python ${projectDir}/scripts/analyse.py \\
        --input '${dataset}' \\
        --inspection '${inspection}' \\
        --outdir analysis \\
        --condition-column '${params.condition_column}' \\
        --control-label '${params.control_label}' \\
        --matrix '${params.matrix}' \\
        --input-kind '${params.input_kind}' \\
        --subset-column '${params.subset_column}' \\
        --subset-value '${params.subset_value}' \\
        --min-genes ${params.min_genes} \\
        --min-counts ${params.min_counts} \\
        --min-cells-per-gene ${params.min_cells_per_gene} \\
        --min-group-cells ${params.min_group_cells}
    """
}

// Real-data, dataset-specific descriptive analysis. The input X is already
// log1p-normalised; this process does not perform count QC or re-normalisation.
process RA_CONTRAST {
    tag { dataset.simpleName }
    conda "${projectDir}/environment.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path dataset

    output:
    path 'ra_analysis'

    script:
    """
    python ${projectDir}/scripts/compare_ra.py \\
        --input '${dataset}' \\
        --outdir ra_analysis \\
        --bmp4 ${params.ra_bmp4} \\
        --egf-bfgf ${params.ra_egf_bfgf} \\
        --mix ${params.ra_mix} \\
        --ra ${params.ra_dose} \\
        --batches '${params.ra_batches}' \\
        --min-cells ${params.min_group_cells}
    """
}

// Verify that the RA summary reflects the input dataset and report
// gene detection prevalence and descriptive batch concordance.
process RA_QUALITY {
    tag { dataset.simpleName }
    conda "${projectDir}/environment.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path dataset
    path ra_results

    output:
    path 'ra_quality'

    script:
    """
    python ${projectDir}/scripts/ra_quality.py \\
        --input '${dataset}' \\
        --summary '${ra_results}/contrast_summary.json' \\
        --genes '${ra_results}/gene_differences.tsv' \\
        --outdir ra_quality
    """
}

workflow {
    if (!params.input) error 'Please supply --input path/to/file.h5ad'
    if (!(params.mode in ['inspect', 'analyse', 'ra'])) error '--mode must be inspect, analyse, or ra'

    // This starter handles one dataset per run to avoid ambiguous sample pairing.
    def datasets = Channel.fromPath(params.input, checkIfExists: true)
        .collect()
        .map { files ->
            if (files.size() != 1) error 'Supply exactly one .h5ad input file'
            files[0]
        }

    inspected = INSPECT(datasets)
    if (params.mode == 'analyse') {
        ANALYSE(inspected)
    }
    else if (params.mode == 'ra') {
        ra_results = RA_CONTRAST(datasets)
        RA_QUALITY(datasets, ra_results)
    }
}
