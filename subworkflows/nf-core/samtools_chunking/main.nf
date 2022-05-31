//
// Splits bam files into chunks using joint intervals file containing genomic intervals.
//

include { SAMTOOLS_VIEW } from '../../../modules/samtools/view/main.nf'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/local/picard/addorreplacereadgroups/main.nf'

workflow SAMTOOLS_CHUNK {
    take:
    input           // channel: [ val(meta), [ input ], indexes, intervals ]
    joint_intervals // channel: /path/to/joint/intervals/file

    main:
    ch_versions = Channel.empty()
    // joint intervals file read in and rows combined with the input channel.
    intervals_file = Channel.from(joint_intervals).splitCsv(sep: "\t")
    ch_added_ints = input.combine(intervals_file)
    // rows from joint intervals file mapped to crate intervals values
    ch_input = ch_added_ints.map {
        meta, reads, indexes, chr, start, end, seconds, time ->
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_${chr}_${start}_${end}"
        new_meta.splits = "${chr}:${start}-${end}"
        [new_meta, reads, indexes, "${chr}:${start}-${end}" ]
    }
    ch_indexes = ch_input.map {
       meta, reads, indexes, intervals ->
       [meta, indexes]
    }
    ch_intervals = ch_input.map {
       meta, reads, indexes, intervals ->
       [meta, intervals]
    }

    //
    //Bam files split into chunks using samtools view
    //
    SAMTOOLS_VIEW (ch_input, [])
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    //chunked bams and their read group info sent to addorreplacereadgroups
    ch_picard_in = SAMTOOLS_VIEW.out.bam.map{meta, bam ->
    ID = "$meta.rgID"
    LB = "$meta.rgLB"
    PL = "$meta.rgPL"
    PU = "$meta.rgPU"
    SM = "$meta.id"
    [meta, bam, ID, LB, PL, PU, SM]}
    PICARD_ADDORREPLACEREADGROUPS(ch_picard_in)
    ch_reformmatted = PICARD_ADDORREPLACEREADGROUPS.out.bam.combine(ch_indexes, by: 0).combine(ch_intervals, by: 0)

    emit:
    ch_format_out   = ch_reformmatted
    // ch_bam_out      = SAMTOOLS_VIEW.out.bam // channel: [ file.bam ]
    ch_rg_bam_out = PICARD_ADDORREPLACEREADGROUPS.out.bam
    // ch_index_out    = ch_indexes            // channel: [ file.bai ]
    ch_interval_out = ch_intervals          // channel: [ interval_string ]
    versions        = ch_versions           // channel: [ versions.yml ]

}
