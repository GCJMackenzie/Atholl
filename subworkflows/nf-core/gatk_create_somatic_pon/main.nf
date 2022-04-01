//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//

include { GATK4_MUTECT2                     } from '../../../modules/gatk4/mutect2/main'
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/gatk4/genomicsdbimport/main'
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/gatk4/createsomaticpanelofnormals/main'

workflow GATK_CREATE_SOMATIC_PON {
    take:
    gendb_input         // channel: [ val(meta), [ vcf ], [ tbi ], intervals, [] ]
    fasta               // channel: /path/to/reference/fasta
    fai                 // channel: /path/to/reference/fasta/index
    dict                // channel: /path/to/reference/fasta/dictionary
    pon_name            // channel: name for panel of normals
    joint_interval_file       // channel: /path/to/interval/file

    main:
    ch_versions      = Channel.empty()

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    //
    //Panel of normals made from genomicsdb workspace using createsomaticpanelofnormals.
    //
    GATK4_GENOMICSDBIMPORT.out.genomicsdb.view()
    GATK4_CREATESOMATICPANELOFNORMALS ( GATK4_GENOMICSDBIMPORT.out.genomicsdb, fasta, fai, dict )
    ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS.out.versions.first())

    emit:

    genomicsdb       = GATK4_GENOMICSDBIMPORT.out.genomicsdb               // channel: [ val(meta), [ genomicsdb ] ]

    pon_vcf          = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf           // channel: [ val(meta), [ vcf.gz ] ]
    pon_index        = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi           // channel: [ val(meta), [ tbi ] ]

    versions         = ch_versions                                         // channel: [ versions.yml ]
}
