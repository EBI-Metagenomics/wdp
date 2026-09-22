process GEMSPARCL_CLUSTER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "quay.io/microbiome-informatics/gemsparcl@sha256:a5b6296d9d4ce3f7d03845574f75baac72fbe4ea0737762824aef9a5c8d7741a"

    input:
    path genomes_file
    path completeness_file

    output:
    path("wdp_clusters.csv")                       , emit: clusters
    path("wdp_cluster_stats.txt")                  , emit: stats
    path("wdp.log")                                , emit: log
    path("wdp_representatives.txt"), optional: true, emit: representatives
    path("wdp.skm")                                , emit: sketch_skm
    path("wdp.skd")                                , emit: sketch_skd
    path("wdp.dists")                              , emit: dists
    tuple val("${task.process}"), val('gemsparcl'), eval("gemsparcl --version 2>&1 | tail -n1"), emit: versions_gemsparcl, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    gemsparcl cluster \\
        -i ${genomes_file} \\
        -o wdp \\
        --completeness-file ${completeness_file} \\
        --threads ${task.cpus} \\
        --representatives \\
        ${args}\\
        2>&1 | tee wdp.log
    """

    stub:
    """
    touch wdp_clusters.csv wdp_cluster_stats.txt wdp.log wdp_representatives.txt wdp.skm wdp.skd wdp.dists
    """
}
