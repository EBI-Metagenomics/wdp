process GUNC_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'biocontainers/gawk:5.3.1' }"

    input:
    tuple val(meta), path(maxcss_tsv)
    path completeness_tsv

    output:
    tuple val(meta), path("*.gunc_filter.tsv"), emit: result
    tuple val("${task.process}"), val('gawk'), eval("awk --version | head -n1"), emit: versions_gawk, topic: versions

    script:
    // Reproduces EBI-Metagenomics/genomes-catalogue-pipeline's modules/gunc.nf filtering logic. 
    
    """
    awk -F'\\t' '
        NR==FNR {
            if (FNR==1) { for(i=1;i<=NF;i++) ch[\$i]=i; next }
            comp[\$(ch["genome"])] = \$(ch["completeness"])
            next
        }
        FNR==1 { for(i=1;i<=NF;i++) gh[\$i]=i; print "genome\\tgunc_contaminated\\tgunc_excluded"; next }
        {
            g = \$(gh["genome"])
            contaminated = (\$(gh["clade_separation_score"]) > 0.45 && \$(gh["contamination_portion"]) > 0.05 && \$(gh["reference_representation_score"]) > 0.5) ? "true" : "false"
            excluded = (contaminated == "true" && (g in comp) && comp[g] < 90) ? "true" : "false"
            print g "\\t" contaminated "\\t" excluded
        }
    ' ${completeness_tsv} ${maxcss_tsv} > ${meta.id}.gunc_filter.tsv
    """

    stub:
    """
    printf "genome\\tgunc_contaminated\\tgunc_excluded\\n" > ${meta.id}.gunc_filter.tsv
    """
}
