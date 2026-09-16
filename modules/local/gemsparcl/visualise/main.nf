process GEMSPARCL_VISUALISE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
   
    container "quay.io/microbiome-informatics/gemsparcl@sha256:a5b6296d9d4ce3f7d03845574f75baac72fbe4ea0737762824aef9a5c8d7741a"

    input:
    path dists
    path clusters_csv

    output:
    path("wdp_vis_part*.graphml")        , emit: graphml
    path("wdp_vis_annotations_part*.csv"), emit: annotations
    tuple val("${task.process}"), val('gemsparcl'), eval("gemsparcl --version 2>&1 | tail -n1"), emit: versions_gemsparcl, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    gemsparcl visualise \\
        --existing-distances ${dists} \\
        --clusters-file ${clusters_csv} \\
        -o wdp_vis \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch wdp_vis_part1.graphml wdp_vis_annotations_part1.csv
    """
}
