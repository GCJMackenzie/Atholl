//
// Performs GATK best practice pre-processing of reads using markduplicates, sortsam, samtools index and BQSR.
//

include { GATK4_APPLYBQSR                                   } from '../../../modules/gatk4/applybqsr/main.nf'
include { GATK4_BASERECALIBRATOR                            } from '../../../modules/gatk4/baserecalibrator/main.nf'
include { PICARD_MARKDUPLICATES                             } from '../../../modules/picard/markduplicates/main.nf'
include { PICARD_SORTSAM as PICARD_SORTSAM_DUPLICATESMARKED } from '../../../modules/picard/sortsam/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RECAL            } from '../../../modules/samtools/index/main.nf'

workflow GATK_PREPROCESS {
    take:
    input          // channel: [ val(meta), [ bam ], [ bai ], intervals ]
    fasta          // channel: /path/to/reference/fasta
    fai            // channel: /path/to/reference/fasta/index
    dict           // channel: /path/to/reference/fasta/dictionary
    sort_order     // channel: which sort order to use for PICARD_SORTSAM_DUPLICATESMARKED
    knownsites     // channel: /path/to/known/sites/vcf
    knownsites_tbi // channel: /path/to/known/sites/tbi

    main:
    ch_versions = Channel.empty()

    ch_input = input.map {
        meta, reads, indexes, intervals ->
        [meta, reads]
    }

    ch_bais = input.map {
        meta, reads, indexes, intervals ->
        [meta, indexes]
    }
    ch_intervals = input.map {
        meta, reads, indexes, intervals ->
        [meta, intervals]
    }

    //
    //use picard markduplicates to mark duplicates in the alignment bams
    //
    PICARD_MARKDUPLICATES ( ch_input )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    ch_sortsam_in = PICARD_MARKDUPLICATES.out.bam

    //
    //Bam files sorted using picard sortsam.
    //
    PICARD_SORTSAM_DUPLICATESMARKED ( ch_sortsam_in, sort_order )
    ch_versions = ch_versions.mix(PICARD_SORTSAM_DUPLICATESMARKED.out.versions)
    ch_samindex_in = PICARD_SORTSAM_DUPLICATESMARKED.out.bam

    //
    //Index for sorted bam file made using samtools index
    //
    SAMTOOLS_INDEX_RECAL (ch_samindex_in)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_RECAL.out.versions)
    ch_bai = SAMTOOLS_INDEX_RECAL.out.bai

    //
    //Perform first pass of BQSR using gatk baserecalibrator.
    //
    ch_baserecal_in = ch_samindex_in.combine(ch_bai, by: 0).combine(ch_intervals, by: 0)
    GATK4_BASERECALIBRATOR( ch_baserecal_in, fasta, fai, dict, knownsites, knownsites_tbi )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)
    ch_bqsrtable = GATK4_BASERECALIBRATOR.out.table

    //
    //Perform second pass of BQSR using gatk applybqsr.
    //
    ch_baserecal_out = ch_baserecal_in.combine(ch_bqsrtable, by: 0)
    ch_bqsr_in = ch_baserecal_out.map {
       meta, input, input_index, intervals, bqsrtable ->
       [meta, input, input_index, bqsrtable, intervals]
    }
    GATK4_APPLYBQSR( ch_bqsr_in, fasta, fai, dict )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    emit:
    versions                = ch_versions                                       // channel: [ versions.yml ]
    // markdup_out             = PICARD_MARKDUPLICATES.out.bam                     // channel: [ val(meta), [ bam ] ]
    metrics_out             = PICARD_MARKDUPLICATES.out.metrics                 // channel: [ val(meta), [ metrics ] ]
    // samtools_index_out      = SAMTOOLS_INDEX_RECAL.out.bai                      // channel: [ val(meta), [ bai ] ]
    // baserecalibrator_out    = GATK4_BASERECALIBRATOR.out.table                  // channel: [ val(meta), [ table ] ]
    applybqsr_index_out     = GATK4_APPLYBQSR.out.bai
    applybqsr_out           = GATK4_APPLYBQSR.out.bam                           // channel: [ val(meta), [ bam ] ]
    ch_intervals_out        = ch_intervals
    // sortsam_dupesmarked_out = PICARD_SORTSAM_DUPLICATESMARKED.out.bam           // channel: [ val(meta), [ bam ] ]
}
