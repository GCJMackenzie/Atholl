//
// Run GATK mutect2
//

include { GATK4_MUTECT2                         } from '../../../modules/gatk4/mutect2/main'
include { GATK4_MUTECT2 as GATK4_MUTECT2_SOMPON } from '../../../modules/gatk4/mutect2/main'

workflow GATK_MUTECT2_CALLING {
    take:
    input                 // channel: [ val(meta), [ input ], [ input_index ], [intervals], [which_norm] or [] ]
    tumor_normal          // channel: true/false
    tumor_only            // channel: true/false
    create_sompon        // channel: true/false
    fasta                 // channel: /path/to/reference/fasta
    fai                   // channel: /path/to/reference/fasta/index
    dict                  // channel: /path/to/reference/fasta/dictionary
    germline_resource     // channel: /path/to/germline/resource
    germline_resource_tbi // channel: /path/to/germline/index
    panel_of_normals      // channel: /path/to/panel/of/normals
    panel_of_normals_tbi  // channel: /path/to/panel/of/normals/index

    main:
    ch_intervals = input.map{meta, bam, bai, intervals, which_norm -> [meta, intervals]}
    ch_versions = Channel.empty()
    //
    //Perform variant calling for each input bam using mutect2 module in the selected mode.
    //
    if (tumor_normal){
        ch_input = input.map{old_meta, sample_id, bam, bai, intervals, which_norm ->
            meta = [:]
            meta.id = "${sample_id}_${intervals}"
            meta.intervals = intervals
            [meta, bam, bai, intervals, which_norm[0]]
        }
        ch_input.view()
    } else {
        ch_input = input
    }
    if (create_sompon){
        GATK4_MUTECT2_SOMPON ( ch_input, tumor_only, create_sompon, false, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
        ch_versions = ch_versions.mix(GATK4_MUTECT2_SOMPON.out.versions.first())
    } else {
        GATK4_MUTECT2 ( ch_input, tumor_only, create_sompon, false, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
        ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())
    }

    emit:
    mutect2_vcf   = create_sompon ? GATK4_MUTECT2_SOMPON.out.vcf : GATK4_MUTECT2.out.vcf                      // channel: [ val(meta), [ vcf ] ]
    mutect2_tbi   = create_sompon ? GATK4_MUTECT2_SOMPON.out.tbi : GATK4_MUTECT2.out.tbi                      // channel: [ val(meta), [ tbi ] ]
    mutect2_stats = create_sompon ? GATK4_MUTECT2_SOMPON.out.stats : GATK4_MUTECT2.out.stats                    // channel: [ val(meta), [ stats ] ]
    mutect2_intervals = ch_intervals
    mutect2_f1r2  = tumor_normal ? GATK4_MUTECT2.out.f1r2 : [] // channel: [ val(meta), [ f1r2 ] ]
    versions      = ch_versions                                // channel: [ versions.yml ]
}
