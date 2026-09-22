//
// Cluster QC-passing genomes with gemsparcl, then export the clustering as a Cytoscape-viewable
// network. GEMSPARCL_CLUSTER does the actual sketching/distance/clustering work;
// GEMSPARCL_VISUALISE is a separate, re-runnable step that rebuilds a network view
// from CLUSTER's .dists + clusters.csv without re-sketching.
//

include { GEMSPARCL_CLUSTER    } from '../../../modules/local/gemsparcl/cluster/main'
include { GEMSPARCL_VISUALISE  } from '../../../modules/local/gemsparcl/visualise/main'

workflow GEMSPARCL_CLUSTERING {

    take:
    ch_passing_genomes    // channel: [ meta (+completeness [0-100 scale]), fasta ] --
                          //          GENOME_QC.out.genomes_passing_qc, already filtered to
                          //          passes_qc==true

    main:
    
    // gemsparcl input: genome_id<TAB>path list 
    ch_genome_paths = ch_passing_genomes
        .map { meta, fasta -> "${meta.id}\t${fasta}" }
        .collectFile(name: 'gemsparcl_rfile.tsv', newLine: true, sort: true)

    // Completeness file: genome_id<TAB>completeness. Gemsparcl expects completeness on a 0-1 scale, 
    // which is not CheckM2's native 0-100 scale, hence the /100.0 conversion below.
    ch_completeness_file = ch_passing_genomes
        .map { meta, fasta -> "${meta.id}\t${meta.completeness / 100.0}" }
        .collectFile(name: 'gemsparcl_completeness.tsv', newLine: true, sort: true)

    GEMSPARCL_CLUSTER(ch_genome_paths, ch_completeness_file)

    GEMSPARCL_VISUALISE(GEMSPARCL_CLUSTER.out.dists, GEMSPARCL_CLUSTER.out.clusters)

    emit:
    clusters        = GEMSPARCL_CLUSTER.out.clusters         // channel: path(wdp_clusters.csv)
    cluster_stats   = GEMSPARCL_CLUSTER.out.stats             // channel: path(wdp_cluster_stats.txt)
    representatives = GEMSPARCL_CLUSTER.out.representatives   // channel: path(wdp_representatives.txt)
    sketch_skm      = GEMSPARCL_CLUSTER.out.sketch_skm        // channel: path(wdp.skm)
    sketch_skd      = GEMSPARCL_CLUSTER.out.sketch_skd        // channel: path(wdp.skd)
    dists           = GEMSPARCL_CLUSTER.out.dists             // channel: path(wdp.dists)
    graphml         = GEMSPARCL_VISUALISE.out.graphml         // channel: path(wdp_vis_part*.graphml)
    annotations     = GEMSPARCL_VISUALISE.out.annotations     // channel: path(wdp_vis_annotations_part*.csv)
}
