process GTRANSLATE_DOWNLOADMODELS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/microbiome-informatics/gtranslate@sha256:549f007342ce8a7e19bc411c60c172610ce6feab6e29dd88026dd1c3487b2147"
   
    output:
    path("gtranslate_models")                                  , emit: model_dir
    tuple val("${task.process}"), val('gtranslate'), eval("gtranslate --version | sed -E 's/.*version ([0-9.]+).*/\\1/'"), emit: versions_gtranslate, topic: versions

    script:
    """
    python3 -c "
    import requests

    url = 'https://data.gtdb.ecogenomic.org/tools/gtranslate/gtranslate_r232_classifiers.tar.gz'
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open('gtranslate_r232_classifiers.tar.gz', 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    "
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
