//
// Performs GATK best practice alignment and pre-processing of reads using BWA, GATK mergebamalignments (where necessary), markduplicates, sortsam, samtools index and BQSR.
// BWA index created from fasta file if not already provided
//

include { SAMTOOLS_VIEW } from '../../../modules/samtools/view/main.nf'

workflow SAMTOOLS_CHUNK {
    take:
    input           // channel: [ val(meta), [ input ], indexes, intervals ]
    joint_intervals

    main:
    ch_versions = Channel.empty()
    intervals_to_split = Channel.from(joint_intervals).splitCsv(sep: "\t").map{row ->
        "${row[0]}:${row[1]}-${row[2]}"}
    ch_added_ints = input.combine(intervals_to_split)
    ch_input = ch_added_ints.map {
        meta, reads, indexes, intervals ->
        new_meta = meta.clone()
        new_meta.id = "${meta.id}_$intervals"
        new_meta.splits = intervals
        [new_meta, reads, indexes, intervals]
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
    //Index for sorted bam file made using samtools index
    //
    SAMTOOLS_VIEW (ch_input, [])
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_reformmatted = SAMTOOLS_VIEW.out.bam.combine(ch_indexes, by: 0).combine(ch_intervals, by: 0)

    emit:
    ch_format_out   = ch_reformmatted
    ch_bam_out      = SAMTOOLS_VIEW.out.bam // channel: [ file.bam ]
    ch_index_out    = ch_indexes            // channel: [ file.bai ]
    ch_interval_out = ch_indexes            // channel: [ interval_string ]
    versions        = ch_versions           // channel: [ versions.yml ]

}
