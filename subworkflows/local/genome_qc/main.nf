//
// Genome QC: CheckM2 completeness/contamination, chunked and grouped by resolved codon table
// (CheckM2's --ttable is whole-invocation, so genomes sharing a table are batched together).
//

include { CODON_TABLE_RESOLUTION  } from '../codon_table_resolution/main'
include { CHECKM2                 } from '../../../modules/ebi-metagenomics/checkm2/checkm2/main'
include { GUNC_RUN                } from '../../../modules/nf-core/gunc/run/main'
include { GUNC_FILTER             } from '../../../modules/local/gunc/filter/main'

workflow GENOME_QC {

    take:
    ch_samplesheet          //    channel: [ meta, fasta ] -- meta carries .id, .taxid, .known_ttable
    gtranslate_model_path   //    path: pre-staged gTranslate classifier model dir. Required --
                            //       pre-trained models are expected to exist already; downloading
                            //          them is a separate, standalone workflow's job.
    gtranslate_chunk_size   //    int: max genomes per GTRANSLATE_DETECTTABLE call
    checkm2_db              //    path: pre-staged CheckM2 .dmnd database, required.
    checkm2_chunk_size      //    int: max genomes per CHECKM2 call (within one table-group)
    gunc_db                 //    path: pre-staged GUNC .dmnd database, required.
    gunc_chunk_size         //    int: max genomes per GUNC_RUN call

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
        .map { meta, fasta -> [meta.ttable, meta, fasta] }   // [meta, fasta] -> [table, meta, fasta]
        .groupTuple()                                        // group by table: [table, [meta...], [fasta...]]
        .flatMap { table, metas, fastas ->
            [metas, fastas]
                .transpose()                                 // re-pair each meta with its own fasta
                // Sort deterministically before chunking -- channel arrival order isn't
                // reproducible across runs, and unsorted collate() would shuffle genomes
                // between chunks, busting -resume's cache for no reason.
                .sort { a, b -> a[0].id <=> b[0].id }
                .collate(checkm2_chunk_size)                  // split the sorted group into fixed-size chunks
                .withIndex()                                  // pair each chunk with a stable index
                .collect { chunk, idx ->
                    [ [id: "ttable_${table}_chunk${idx}", ttable: table], chunk.collect { genome -> genome[1] } ]
                }
        }

    ch_checkm2_db = channel.value(file(checkm2_db))

    CHECKM2(ch_checkm2_in, ch_checkm2_db)
    ch_versions = ch_versions.mix(CHECKM2.out.versions)

    ch_checkm2_per_genome = CHECKM2.out.checkm2_stats
        .flatMap { meta, tsv -> tsv.splitCsv(header: true, sep: '\t') }
        .map { row -> [row.Name, row.Completeness as Double, row.Contamination as Double] }

    // GUNC: flat chunking (no codon-table grouping needed)
    ch_gunc_in = ch_all_genomes
        .toList()
        .flatMap { entries ->
            entries
                // Sort deterministically before chunking, same reproducibility reason as
                // CheckM2's chunking above.
                .sort { a, b -> a[0].id <=> b[0].id }
                .collect { meta, fasta -> fasta }
                .collate(gunc_chunk_size)
                .withIndex()
                .collect { chunk, idx ->
                    [ [id: "gunc_chunk${idx}"], chunk ]
                }
        }

    ch_gunc_db = channel.value(file(gunc_db))

    GUNC_RUN(ch_gunc_in, ch_gunc_db)

    ch_completeness_tsv = ch_checkm2_per_genome
        .map { name, completeness, contamination -> "${name}\t${completeness}\t${contamination}" }
        .collectFile(name: 'checkm2_completeness.tsv', newLine: true, sort: true,
            seed: "genome\tcompleteness\tcontamination")

    GUNC_FILTER(GUNC_RUN.out.maxcss_level_tsv, ch_completeness_tsv.first())
   
    ch_gunc_per_genome = GUNC_FILTER.out.result
        .flatMap { meta, tsv -> tsv.splitCsv(header: true, sep: '\t') }
        .map { row -> [row.genome, row.gunc_contaminated == 'true', row.gunc_excluded == 'true'] }

    ch_genomes_with_qc = ch_all_genomes
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(ch_checkm2_per_genome, failOnMismatch: true)
        .join(ch_gunc_per_genome, failOnMismatch: true)
        .map { id, meta, fasta, completeness, contamination, gunc_contaminated, gunc_excluded ->
            // Quality score (QS) = completeness - 5*contamination, matching
            // EBI-Metagenomics/genomes-catalogue-pipeline's formula (bin/filter_qs50.py:
            // `qs50()`). qs50/qs80 apply that same function at two thresholds -- both also
            // require contamination <= 5.0, per that function's own combined check (not just
            // the raw QS score alone).
            def quality_score = Math.round((completeness - (5.0 * contamination)) * 100) / 100
            def passes_qc_80_5 = (completeness >= 80.0 && contamination < 5.0)
            [ meta + [
                completeness: completeness,
                contamination: contamination,
                passes_qc_80_5: passes_qc_80_5,
                quality_score: quality_score,
                qs50: (contamination <= 5.0 && quality_score >= 50.0),
                qs80: (contamination <= 5.0 && quality_score >= 80.0),
                gunc_contaminated: gunc_contaminated,
                // Final combined QC gate -- gunc_excluded is GUNC_FILTER's own awk-computed
                // intersection (gunc_contaminated AND completeness<90), matching
                // genomes-catalogue-pipeline's real bad.txt logic exactly.
                passes_qc: passes_qc_80_5 && !gunc_excluded
            ], fasta ]
        }

    ch_genomes_passing_qc = ch_genomes_with_qc.filter { meta, fasta -> meta.passes_qc }

    emit:
    genomes_with_qc    = ch_genomes_with_qc    // channel: [ meta (+completeness/contamination/
                                               //           passes_qc_80_5/quality_score/qs50/qs80/
                                               //           gunc_contaminated/passes_qc), fasta ]
    genomes_passing_qc = ch_genomes_passing_qc // channel: same shape, filtered to passes_qc==true
    versions           = ch_versions
}
