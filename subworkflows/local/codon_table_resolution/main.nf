//
// Resolve each genome's codon (translation) table: per-row known_ttable override ->
// gTranslate fallback for whatever's left unresolved.
//
// Deliberately a standalone subworkflow (not inlined into the larger genome_qc_and_clustering
// subworkflow) so it can be wired and tested on its own before CheckM2/GUNC/gemsparcl modules
// exist -- genome_qc_and_clustering will include and extend this once those land.
//

include { GTRANSLATE_DOWNLOADMODELS } from '../../../modules/local/gtranslate/downloadmodels/main'
include { GTRANSLATE_DETECTTABLE    } from '../../../modules/local/gtranslate/detecttable/main'

workflow CODON_TABLE_RESOLUTION {

    take:
    ch_samplesheet          //    channel: [ meta, fasta ] -- meta carries .id, .taxid, .known_ttable
    gtranslate_model_path   //    path: pre-staged gTranslate classifier model dir, or null to
                            //          trigger GTRANSLATE_DOWNLOADMODELS
    gtranslate_chunk_size   //    int: max genomes per GTRANSLATE_DETECTTABLE call

    main:

    // known_ttable override -> gTranslate fallback for everything else
    ch_branched = ch_samplesheet
        .map { meta, fasta ->
            def resolved = meta.known_ttable ? meta.known_ttable as Integer : null
            [ meta + [known_table: resolved], fasta ]
        }
        .branch {
            override_hit:  it[0].known_table != null
            override_miss: true
        }

    // gTranslate fallback for whatever wasn't resolved by the per-row override
    ch_model_dir = gtranslate_model_path
        ? Channel.value(file(gtranslate_model_path))
        : GTRANSLATE_DOWNLOADMODELS().model_dir

    ch_gtranslate_in = ch_branched.override_miss
        .map { meta, fasta -> fasta }
        .toList()
        .flatMap { fastas ->
            fastas.collate(gtranslate_chunk_size).withIndex().collect { chunk, idx ->
                [ [id: "gtranslate_chunk${idx}"], chunk ]
            }
        }
    GTRANSLATE_DETECTTABLE(ch_gtranslate_in, ch_model_dir)

    ch_detected = GTRANSLATE_DETECTTABLE.out.translation_table_summary
        .flatMap { meta, tsv -> tsv.splitCsv(header: true, sep: '\t') }
        .map { row -> [row.user_genome, row.best_tln_table as Integer] }

    ch_resolved_miss = ch_branched.override_miss
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_detected, failOnMismatch: true)
        .map { id, meta, fasta, table -> [meta + [known_table: table], fasta] }

    ch_genomes_with_table = ch_branched.override_hit.mix(ch_resolved_miss)

    emit:
    genomes_with_table = ch_genomes_with_table   // channel: [ meta (+known_table), fasta ]
}
