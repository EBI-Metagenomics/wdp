process GTRANSLATE_DETECTTABLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/microbiome-informatics/gtranslate@sha256:549f007342ce8a7e19bc411c60c172610ce6feab6e29dd88026dd1c3487b2147"
   
    input:
    tuple val(meta), path(fastas, stageAs: "batch_input/*")
    path model_dir

    output:
    tuple val(meta), path("results/${meta.id}.translation_table_summary.tsv"), emit: translation_table_summary
    tuple val("${task.process}"), val('gtranslate'), eval("gtranslate --version | sed -E 's/.*version ([0-9.]+).*/\\1/'"), emit: versions_gtranslate, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    for f in batch_input/*; do
        acc=\$(basename "\$f" | sed -E 's/\\.fasta\\.gz\$//')
        printf "%s\\t%s\\n" "\$f" "\$acc"
    done > batchfile.tsv

    GTRANSLATE_MODEL_PATH=${model_dir} gtranslate detect_table \\
        --batchfile batchfile.tsv \\
        --out_dir results \\
        --cpus ${task.cpus} \\
        --prefix ${meta.id} \\
        ${args}
    """

    stub:
    """
    mkdir -p results
    touch results/${meta.id}.translation_table_summary.tsv
    """
}
