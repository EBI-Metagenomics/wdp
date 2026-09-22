process GUNC_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(maxcss_tsv)
    path completeness_tsv

    output:
    tuple val(meta), path("*.gunc_filter.tsv"), emit: result
    tuple val("${task.process}"), val('python'), eval("python3 --version"), emit: versions_python, topic: versions

    script:
    // Reproduces EBI-Metagenomics/genomes-catalogue-pipeline's modules/gunc.nf filtering logic.
    """
    python3 <<'PYEOF'
    import csv
    import math

    def is_num(x):
        try:
            return not math.isnan(float(x))
        except ValueError:
            return False

    completeness = {}
    with open("${completeness_tsv}") as f:
        for row in csv.DictReader(f, delimiter='\\t'):
            completeness[row['genome']] = row['completeness']

    with open("${maxcss_tsv}") as f, open("${meta.id}.gunc_filter.tsv", 'w', newline='') as out:
        reader = csv.DictReader(f, delimiter='\\t')
        writer = csv.writer(out, delimiter='\\t')
        writer.writerow(['genome', 'gunc_contaminated', 'gunc_excluded'])
        for row in reader:
            css = row['clade_separation_score']
            cp = row['contamination_portion']
            rrs = row['reference_representation_score']
            # GUNC can write "nan" for small genomes with few genes -- is_num() guards against
            # treating that as a valid, thresholdable score.
            contaminated = (
                is_num(css) and is_num(cp) and is_num(rrs)
                and float(css) > 0.45 and float(cp) > 0.05 and float(rrs) > 0.5
            )
            genome = row['genome']
            # re-visit the filtering logic here
            excluded = contaminated and genome in completeness and float(completeness[genome]) < 90
            writer.writerow([genome, str(contaminated).lower(), str(excluded).lower()])
    PYEOF
    """

    stub:
    """
    printf "genome\\tgunc_contaminated\\tgunc_excluded\\n" > ${meta.id}.gunc_filter.tsv
    """
}
