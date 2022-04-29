//
// Run GATK haplotypecaller.
//

include { GATK4_HAPLOTYPECALLER } from '../../../modules/gatk4/haplotypecaller/main'
include { PICARD_RENAMESAMPLEINVCF } from '../../../modules/picard/renamesampleinvcf/main'

workflow GATK_HAPLOTYPECALLING {
    take:
    input            // channel: [ val(meta), [ input ], [ input_index ], [ intervals ] ]
    fasta            // channel: /path/to/reference/fasta
    fai              // channel: /path/to/reference/fasta/index
    dict             // channel: /path/to/reference/fasta/dictionary
    sites            // channel: /path/to/known/sites/file
    sites_index      // channel: /path/to/known/sites/index

    main:
    ch_versions = Channel.empty()

    haplotc_input = input
    haplotc_intervals = input.map{meta, bam, bai, intervals -> [meta, intervals] }
    //
    //Perform variant calling using haplotypecaller module. Additional argument "-ERC GVCF" used to run in gvcf mode.
    //
    GATK4_HAPLOTYPECALLER ( haplotc_input, fasta, fai, dict, sites, sites_index)

    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    PICARD_RENAMESAMPLEINVCF ( GATK4_HAPLOTYPECALLER.out.vcf )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
    emit:
    versions       = ch_versions                      // channel: [ versions.yml ]
    haplotc_vcf    = GATK4_HAPLOTYPECALLER.out.vcf    // channel: [ val(meta), [ vcf ] ]
    haplotc_index  = GATK4_HAPLOTYPECALLER.out.tbi    // channel: [ val(meta), [ tbi ] ]
    renamed_vcf    = PICARD_RENAMESAMPLEINVCF.out.vcf    // channel: [ val(meta), [ vcf ] ]
    renamed_index  = PICARD_RENAMESAMPLEINVCF.out.tbi    // channel: [ val(meta), [ tbi ] ]
    haplotc_interval_out = haplotc_intervals
}
