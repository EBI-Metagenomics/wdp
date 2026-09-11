//
// Genome QC: CheckM2 completeness/contamination, chunked and grouped by resolved codon table
// (CheckM2's --ttable is whole-invocation, so genomes sharing a table are batched together).


include { CODON_TABLE_RESOLUTION  } from '../codon_table_resolution/main'
include { CHECKM2                 } from '../../../modules/ebi-metagenomics/checkm2/checkm2/main'
include { CHECKM2_DOWNLOAD_DB     } from '../../../modules/ebi-metagenomics/checkm2/download_db/main'  

workflow GENOME_QC {

    take:
    ch_samplesheet          // channel: [ meta, fasta ] -- meta carries .id, .taxid, .known_ttable
    gtranslate_model_path   //    path: pre-staged gTranslate classifier model dir, or null to trigger GTRANSLATE_DOWNLOADMODELS
    gtranslate_chunk_size   //     int: max genomes per GTRANSLATE_DETECTTABLE call
    checkm2_db              //    path: pre-staged CheckM2 .dmnd database, or null to trigger CHECKM2_DOWNLOAD_DB
    checkm2_chunk_size      //     int: max genomes per CHECKM2 call (within one table-group)

    main:

    def ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Resolve each genome's codon (translation) table
    //
    CODON_TABLE_RESOLUTION(
        ch_samplesheet,
        gtranslate_model_path,
        gtranslate_chunk_size
    )

    ch_all_genomes = CODON_TABLE_RESOLUTION.out.genomes_with_table

    // CheckM2: group by resolved table (mandatory -- --ttable is whole-invocation), then chunk
    ch_checkm2_in = ch_all_genomes
        .map { meta, fasta -> [meta.known_table, meta, fasta] }
        .groupTuple()
        .flatMap { table, metas, fastas ->
            [metas, fastas].transpose().collate(checkm2_chunk_size).withIndex().collect { chunk, idx ->
                [ [id: "ttable_${table}_chunk${idx}", ttable: table], chunk.collect { it[1] } ]
            }
        }

    if (checkm2_db) {
        ch_checkm2_db = Channel.value(file(checkm2_db))
    } else {
        CHECKM2_DOWNLOAD_DB()
        ch_checkm2_db = CHECKM2_DOWNLOAD_DB.out.checkm2_db
        ch_versions = ch_versions.mix(CHECKM2_DOWNLOAD_DB.out.versions)
    }

    CHECKM2(ch_checkm2_in, ch_checkm2_db)
    ch_versions = ch_versions.mix(CHECKM2.out.versions)

    ch_checkm2_per_genome = CHECKM2.out.checkm2_stats
        .flatMap { meta, tsv -> tsv.splitCsv(header: true, sep: '\t') }
        .map { row -> [row.Name, row.Completeness as Double, row.Contamination as Double] }

    ch_genomes_with_qc = ch_all_genomes
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_checkm2_per_genome)
        .map { id, meta, fasta, completeness, contamination ->
            // Quality score (QS) = completeness - 5*contamination, matching
            // EBI-Metagenomics/genomes-catalogue-pipeline's real formula (bin/filter_qs50.py:
            // `qs50()`). qs50/qs80 apply that same function at two thresholds -- both also
            // require contamination <= 5.0, per that function's own combined check (not just
            // the raw QS score alone).
            def quality_score = completeness - (5.0 * contamination)
            [ meta + [
                completeness: completeness,
                contamination: contamination,
                passes_qc_80_5: (completeness >= 80.0 && contamination < 5.0),
                quality_score: quality_score,
                qs50: (contamination <= 5.0 && quality_score >= 50.0),
                qs80: (contamination <= 5.0 && quality_score >= 80.0)
            ], fasta ]
        }

    emit:
    genomes_with_qc = ch_genomes_with_qc   // channel: [ meta (+completeness/contamination/
                                           //           passes_qc_80_5), fasta ]
    versions        = ch_versions
}
