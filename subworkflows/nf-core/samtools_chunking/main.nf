//
// Performs GATK best practice alignment and pre-processing of reads using BWA, GATK mergebamalignments (where necessary), markduplicates, sortsam, samtools index and BQSR.
// BWA index created from fasta file if not already provided
//

include { SAMTOOLS_VIEW } from '../../../modules/samtools/view/main.nf'

workflow SAMTOOLS_CHUNK {
    take:
    input           // channel: [ val(meta), [ input ], indexes, intervals ]
    joint_intervals // channel: /path/to/joint/intervals/file

    main:
    ch_versions = Channel.empty()
    println("meep")
    ch_input = input.map {
        meta, reads, indexes, intervals ->
        [meta, reads, indexes]
    }

    ch_bams = input.map {
        meta, reads, indexes, intervals ->
        [meta, indexes]
    }
    ch_intervals = input.map {
        meta, reads, indexes, intervals ->
        [meta, intervals]
    }

    //
    //Index for sorted bam file made using samtools index
    //
    SAMTOOLS_VIEW (ch_input, [])
    // ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    // ch_bai = SAMTOOLS_VIEW.out.bai

    //
    //Perform first pass of BQSR using gatk baserecalibrator.
    //
    //ch_baserecal_in = ch_samindex_in.combine(ch_bai, by: 0).combine(ch_intervals, by: 0)


    emit:
    versions                = ch_versions                                       // channel: [ versions.yml ]

}
