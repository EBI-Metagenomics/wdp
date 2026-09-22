//
// Resolve each genome's codon (translation) table: per-row known_ttable override ->
// gTranslate fallback for whatever's left unresolved.
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
            [ meta + [ttable: resolved], fasta ]
        }
        .branch { meta, fasta ->
            override_hit:  meta.ttable != null
            override_miss: true
        }

    // gTranslate fallback for whatever wasn't resolved by the per-row override
    ch_model_dir = gtranslate_model_path
        ? channel.value(file(gtranslate_model_path))
        : GTRANSLATE_DOWNLOADMODELS().model_dir

    ch_gtranslate_in = ch_branched.override_miss
        // Collect all [meta, fasta] tuples before chunking
        // This will hurt parallelization as gTranslate will only run once all the genomes have reach this step
        .toList()
        .flatMap { entries ->
            entries
                // Ensure deterministic ordering regardless of channel arrival order.
                .sort { a, b -> a[0].id <=> b[0].id }
                .collect { meta, fasta -> fasta }
                // Split the sorted entries into fixed-size chunks.
                .collate(gtranslate_chunk_size)
                // Add a stable index to each chunk.
                .withIndex()
                .collect { chunk, idx ->
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
        .map { id, meta, fasta, table -> [meta + [ttable: table], fasta] }

    ch_genomes_with_table = ch_branched.override_hit.mix(ch_resolved_miss)

    emit:
    genomes_with_table = ch_genomes_with_table   // channel: [ meta (+ttable), fasta ]
}
