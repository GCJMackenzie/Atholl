#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_ALIGN_AND_PREPROCESS } from './subworkflows/nf-core/gatk_align_and_preprocess/main'
include { GATK_CREATE_SOMATIC_PON } from './subworkflows/nf-core/gatk_create_somatic_pon/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING } from './subworkflows/nf-core/gatk_joint_germline_variant_calling/main'
include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING } from './subworkflows/nf-core/gatk_tumor_only_somatic_variant_calling/main'
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from './subworkflows/nf-core/gatk_tumor_normal_somatic_variant_calling/main'

workflow ATHOLL {

    take:

    //universal args: these args are used by every subworkflow, input arg is used for passing in sample specific data, some of the array entries included are
    //not necessary for everysubworkflow, e.g which_norm, these are only passed in to the subworkflows that require them.
    input
    fasta
    fai
    dict

    //entry params: used to control which subworkflow(s) are run
    alignment
    create_som_pon
    joint_germline
    tumor_somatic
    tumor_normal_somatic
    paired

    //shared args: these args are used by two or more subworkflows, but are not universal, which subworkflows use them are noted.
    sites                 // channel: /path/to/known/sites/file       align and preprocess + joint germline
    sites_index           // channel: /path/to/known/sites/index      align and preprocess + joint germline

    joint_id              // channel: joint id for gendbs and pons    joint germline + create_som_pon
    joint_intervals       // channel: joint intervals file            joint germline + create_som_pon

    germline_resource     // channel: /path/to/germline/resource      tumor_only + tumor_normal
    germline_resource_tbi // channel: /path/to/germline/index         tumor_only + tumor_normal
    panel_of_normals      // channel: /path/to/panel/of/normals       tumor_only + tumor_normal
    panel_of_normals_tbi  // channel: /path/to/panel/of/normals/index tumor_only + tumor_normal

    // aligner args: args exclusive to align and preprocess subworkflow
    bwaindex              // channel: /path/to/bwa/index/directory
    is_ubam               // channel: true/false whether input is in ubam format or not
    sort_order            // channel: which sort order to use for PICARD_SORTSAM_DUPLICATESMARKED

    // joint germline args: args exclusive to joint germline subworkflow
    run_haplotc
    run_vqsr
    allelespecific        // channel: true/false run allelespecific mode of vqsr modules
    resources             // channel: [[resource, vcfs, forvariantrecal], [resource, tbis, forvariantrecal], [resource, labels, forvariantrecal]]
    annotation            // channel: [annotations, to, use, for, variantrecal, filtering]
    mode                  // channel: which mode to run variantrecal: SNP/INDEL/BOTH
    truthsensitivity      // channel: 0-100.0 truthsensitivity cutoff for applyvqsr

    main:
    filetest = extract_samples(input, alignment, paired, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic)
    println(filetest)
    filetest.view()

    if (alignment) {
        println("The aligner is running")
        GATK_ALIGN_AND_PREPROCESS( filetest , fasta , fai , dict , bwaindex , is_ubam , sort_order , sites , sites_index )
    }

    if (create_som_pon) {
        GATK_CREATE_SOMATIC_PON(  filetest , fasta , fai , dict , joint_id, joint_intervals )
    }

    if (joint_germline) {
        println("Performing joint germline variant calling")
        GATK_JOINT_GERMLINE_VARIANT_CALLING(  filetest , run_haplotc , run_vqsr , fasta , fai , dict, sites , sites_index , joint_id , joint_intervals , allelespecific , resources , annotation , mode , false , truthsensitivity )
    }

    if (tumor_somatic) {
        println("Performing tumor-only somatic variant calling")
        GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING(  filetest , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    }

    if (tumor_normal_somatic) {
        println("Performing tumor-normal somatic variant calling")
        GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(  filetest , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    }

}

def extract_samples(csv_file, alignment, paired, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic) {
    firststep = Channel.from(csv_file).splitCsv(header: true).map{ row ->
        def meta = [:]
        meta.id = row.SampleID
        if (alignment) {
            meta.single_end = "$row.single_end"
            meta.read_group = "$row.readgroup"
            if (paired) {
                println("paired end data")
                [meta, [ file(row.input_file , checkIfExists : true), file(row.paired_file_2 , checkIfExists : true) ], row.input_index, file(row.intervals , checkIfExists : true), row.which_norm ]

            } else {
                println("interleaved or ubam")
                [meta, file(row.input_file , checkIfExists : true), row.input_index, file(row.intervals , checkIfExists : true), row.which_norm ]
            }
        } else if (create_som_pon) {
            println("pon")
            [meta, file(row.input_file , checkIfExists : true), file(row.input_index , checkIfExists : true), row.intervals, row.which_norm ]
        } else {
            println("variant calling")
            [meta, file(row.input_file , checkIfExists : true), file(row.input_index , checkIfExists : true), file(row.intervals , checkIfExists : true), row.which_norm ]
        }
    }.groupTuple().map { meta, input_file, input_index, intervals_file, which_norm ->
        def the_meta = meta
        def input_files = input_file
        def input_indexes = input_index
        def intervals = intervals_file
        def which_norms = which_norm
        if (alignment){
            if (paired) {
                return [the_meta, input_files[0], intervals]
            } else {
                return [the_meta, input_files, intervals]
            }
        } else if (joint_germline) {
            return [the_meta, input_files, input_indexes, intervals]
        } else if (tumor_normal_somatic) {
            filtered_normal = which_norms.unique().toList()
            return [the_meta, input_files, input_indexes, intervals, filtered_normal]
        } else {
            return [the_meta, input_files, input_indexes, intervals, which_norms]
        }
    }
}
