//
// Cluster QC-passing genomes with gemsparcl, then export the clustering as a Cytoscape-viewable
// network. Two steps, not one -- GEMSPARCL_CLUSTER does the actual sketching/distance/clustering
// work; GEMSPARCL_VISUALISE is a separate, cheap, re-runnable step that rebuilds a network view
// from CLUSTER's .dists + clusters.csv without re-sketching (per gemsparcl's own real docs,
// https://johannahelene.github.io/gemsparcl/guides/visualise.html -- "you can experiment freely
// with different thresholds and outputs without touching your cluster assignments").
//
// Deliberately its own subworkflow (not folded into GENOME_QC): clustering is a genuinely
// different concern from QC (fact 17-adjacent reasoning -- separate subworkflows where the
// pieces are independently meaningful, not just build-order artifacts).
//

include { GEMSPARCL_CLUSTER    } from '../../../modules/local/gemsparcl/cluster/main'
include { GEMSPARCL_VISUALISE  } from '../../../modules/local/gemsparcl/visualise/main'

workflow GEMSPARCL_CLUSTERING {

    take:
    ch_genomes_with_qc   // channel: [ meta (+passes_qc, +completeness [0-100 scale]), fasta ] --
                          //          GENOME_QC.out.genomes_with_qc, unfiltered (filtering to
                          //          passes_qc==true happens here, not by the caller)

    main:

    // Only genomes that passed the combined QC gate go into clustering.
    ch_passing = ch_genomes_with_qc.filter { meta, fasta -> meta.passes_qc }

    // rfile: genome_id<TAB>path, no header, per gemsparcl's real --input format (docs/source/
    // guides/cluster.rst). Built from the ORIGINAL fasta paths (this collectFile runs at the
    // channel level, before any process staging, so ${fasta} here is the real shared-filesystem
    // path -- e.g. /hps/nobackup/.../ENA_all_fetch/data/<accession>.fasta.gz -- not a
    // work-dir-staged copy. Matches how the samplesheet's own `assembly` column already works;
    // gemsparcl's own subprocess needs filesystem access to these paths at runtime, same
    // assumption the whole pipeline already makes.
    ch_rfile = ch_passing
        .map { meta, fasta -> "${meta.id}\t${fasta}" }
        .collectFile(name: 'gemsparcl_rfile.tsv', newLine: true, sort: true)

    // completeness file: genome_id<TAB>completeness on a 0-1 scale (gemsparcl's real format,
    // confirmed in docs/source/guides/cluster.rst) -- CheckM2 reports 0-100 (fact 1/12), so this
    // is a real conversion, not a passthrough.
    ch_completeness_file = ch_passing
        .map { meta, fasta -> "${meta.id}\t${meta.completeness / 100.0}" }
        .collectFile(name: 'gemsparcl_completeness.tsv', newLine: true, sort: true)

    GEMSPARCL_CLUSTER(ch_rfile, ch_completeness_file)

    GEMSPARCL_VISUALISE(GEMSPARCL_CLUSTER.out.dists, GEMSPARCL_CLUSTER.out.clusters)

    // Both modules self-report via the global versions topic (topic: versions) -- no explicit
    // ch_versions channel needed in this subworkflow at all, unlike GENOME_QC (which has to
    // explicitly mix CheckM2's classic-pattern versions.yml output).

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
