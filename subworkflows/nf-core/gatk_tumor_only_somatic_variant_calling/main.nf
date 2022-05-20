//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

include { GATK4_GETPILEUPSUMMARIES     } from '../../../modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION } from '../../../modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS      } from '../../../modules/gatk4/filtermutectcalls/main'
include { PICARD_RENAMESAMPLEINVCF as PICARD_RENAMESOMATICSAMPLEINVCF} from '../../../modules/picard/renamesampleinvcf/main'

workflow GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING {
    take:
    input                 // channel: [ val(meta), bam, bai, [ vcf ], [ tbi ], [stats] , intervals ]
    fasta                 // channel: /path/to/reference/fasta
    fai                   // channel: /path/to/reference/fasta/index
    dict                  // channel: /path/to/reference/fasta/dictionary
    germline_resource     // channel: /path/to/germline/resource
    germline_resource_tbi // channel: /path/to/germline/index


    main:
    ch_versions = Channel.empty()

    //
    // Generate pileup summary table using getepileupsummaries.
    //
    pileup_input = input.map {
         meta, bam, bai, input_file, input_index, stats, intervals ->
         [meta, bam, bai, intervals]
     }
    GATK4_GETPILEUPSUMMARIES ( pileup_input , germline_resource , germline_resource_tbi )
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions)

    //
    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup = GATK4_GETPILEUPSUMMARIES.out.table.map{meta, table -> [meta, table, []]}
    //[] is a placeholder for the optional input where the matched normal sample would be passed in for tumor-normal samples, which is not necessary for this workflow.
    GATK4_CALCULATECONTAMINATION ( ch_pileup, true )
    ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions)

    //
    // Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    ch_vcf = input
    ch_segment       = GATK4_CALCULATECONTAMINATION.out.segmentation
    ch_contamination = GATK4_CALCULATECONTAMINATION.out.contamination
    ch_filter_results = ch_vcf.combine(ch_segment, by: 0).combine(ch_contamination, by: 0)
    ch_filtermutect_in = ch_filter_results.map{meta, bam, bai, vcf, tbi, stats, intervals, segment, contamination ->
         [meta, vcf, tbi, stats, intervals, [], segment, contamination, []]}
     ch_filtermutect_in.view()
    GATK4_FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)

    PICARD_RENAMESOMATICSAMPLEINVCF ( GATK4_FILTERMUTECTCALLS.out.vcf )
    ch_versions = ch_versions.mix(PICARD_RENAMESOMATICSAMPLEINVCF.out.versions)

    emit:

    // pileup_table        = GATK4_GETPILEUPSUMMARIES.out.table             // channel: [ val(meta), [ table ] ]

    // contamination_table = GATK4_CALCULATECONTAMINATION.out.contamination // channel: [ val(meta), [ contamination ] ]
    // segmentation_table  = GATK4_CALCULATECONTAMINATION.out.segmentation  // channel: [ val(meta), [ segmentation ] ]

    // filtered_vcf        = GATK4_FILTERMUTECTCALLS.out.vcf                // channel: [ val(meta), [ vcf ] ]
    // filtered_index      = GATK4_FILTERMUTECTCALLS.out.tbi                // channel: [ val(meta), [ tbi ] ]
    filtered_stats      = GATK4_FILTERMUTECTCALLS.out.stats              // channel: [ val(meta), [ stats ] ]

    renamed_vcf    = PICARD_RENAMESOMATICSAMPLEINVCF.out.vcf    // channel: [ val(meta), [ vcf ] ]

    versions            = ch_versions                                              // channel: [ versions.yml ]
}
