process GTRANSLATE_DOWNLOADMODELS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    output:
    path("gtranslate_models")                                  , emit: model_dir
    tuple val("${task.process}"), val('wget'), eval("wget --version | head -n1"), emit: versions_wget, topic: versions

    script:
    """
    wget https://data.gtdb.ecogenomic.org/tools/gtranslate/gtranslate_r232_classifiers.tar.gz

    mkdir -p gtranslate_models
    tar -xzf gtranslate_r232_classifiers.tar.gz -C gtranslate_models --strip-components=1
    rm gtranslate_r232_classifiers.tar.gz
    """

    stub:
    """
    mkdir -p gtranslate_models
    touch gtranslate_models/.stub
    """
}
