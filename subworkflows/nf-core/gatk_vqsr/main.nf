//
// Recalibrate with variantrecalibrator & applyvqsr.
//

include { GATK4_VARIANTRECALIBRATOR } from '../../../modules/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR } from '../../../modules/gatk4/applyvqsr/main'

workflow GATK_VQSR {
    take:
    input            // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta            // channel: /path/to/reference/fasta
    fai              // channel: /path/to/reference/fasta/index
    dict             // channel: /path/to/reference/fasta/dictionary
    allelespecific   // channel: true/false run allelespecific mode of vqsr modules
    resources        // channel: [[resource, vcfs, forvariantrecal], [resource, tbis, forvariantrecal], [resource, labels, forvariantrecal]]
    annotation       // channel: [annotations, to, use, for, variantrecal, filtering]
    mode             // channel: which mode to run variantrecal: SNP/INDEL/BOTH
    create_rscript   // channel: true/false whether to generate rscript plots in variantrecal
    truthsensitivity // channel: 0-100.0 truthsensitivity cutoff for applyvqsr

    main:
    ch_versions = Channel.empty()

    ch_vrecal_in  = input

    GATK4_VARIANTRECALIBRATOR ( ch_vrecal_in, fasta, fai, dict, allelespecific, resources, annotation, mode, create_rscript )

    ch_versions   = ch_versions.mix(GATK4_VARIANTRECALIBRATOR.out.versions)

    //
    //Perform second step in VQSR using ApplyVQSR
    //
    ch_recal      = GATK4_VARIANTRECALIBRATOR.out.recal
    ch_idx        = GATK4_VARIANTRECALIBRATOR.out.idx
    ch_tranches   = GATK4_VARIANTRECALIBRATOR.out.tranches
    ch_vqsr_in    = ch_vrecal_in.combine(ch_recal, by: 0).combine(ch_idx, by: 0).combine(ch_tranches, by: 0)

    GATK4_APPLYVQSR ( ch_vqsr_in, fasta, fai, dict, allelespecific, truthsensitivity, mode )

    ch_versions   = ch_versions.mix(GATK4_APPLYVQSR.out.versions)


    emit:
    versions       = ch_versions                                     // channel: [ versions.yml ]
    recal_file     = GATK4_VARIANTRECALIBRATOR.out.recal // channel: [ val(meta), [ recal ] ]
    recal_index    = GATK4_VARIANTRECALIBRATOR.out.idx   // channel: [ val(meta), [ idx ] ]
    recal_tranches = GATK4_VARIANTRECALIBRATOR.out.tranches // channel: [ val(meta), [ tranches ] ]
    vqsr_vcf       = GATK4_APPLYVQSR.out.vcf             // channel: [ val(meta), [ vcf ] ]
    vqsr_index     = GATK4_APPLYVQSR.out.tbi             // channel: [ val(meta), [ tbi ] ]
}
