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
    filetest = Channel.from(input).splitCsv(header: true).map{ row ->
        metaid = "id:$row.SampleID"
        if (alignment) {
            metasingle = "single_end:$row.single_end"
            metaread = "read_group:$row.readgroup"
            if (paired) {
                println("paired end data")
                [[metaid, metasingle, metaread], [ file(row.input_file , checkIfExists : true), file(row.paired_file_2 , checkIfExists : true) ], row.input_index, file(row.intervals , checkIfExists : true), row.which_norm ]
            } else {
                println("interleaved or ubam")
                [[metaid, metasingle, metaread], file(row.input_file , checkIfExists : true), row.input_index, file(row.intervals , checkIfExists : true), row.which_norm ]
            }
        } else {
            println("variant calling or pon")
            [[metaid], file(row.input_file , checkIfExists : true), file(row.input_index , checkIfExists : true), file(row.intervals , checkIfExists : true), row.which_norm ]
        }
    }.groupTuple()

    if (alignment) {
        ch_align_in = filetest.map {
            meta, input_file, input_index, intervals_file, which_normal ->
            if (paired) {
                [meta, input_file[0], intervals_file]
            } else {
                [meta, input_file, intervals_file]
            }
        }.collect()
        ch_align_in.view()
        println("The aligner is running")
        GATK_ALIGN_AND_PREPROCESS( ch_align_in , fasta , fai , dict , bwaindex , is_ubam , sort_order , sites , sites_index )
    }

    if (create_som_pon) {
        ch_sompon_in = filetest.toList()
        println("Panel of normals is being made")
        ch_sompon_in.view()
        // GATK_CREATE_SOMATIC_PON(  ch_sompon_in , fasta , fai , dict , joint_id, joint_intervals )
    }

    if (joint_germline) {
        ch_joint_in = filetest.map {
            meta, reads, index, intervals, which_norm ->
            [meta, reads, index, intervals]
        }.toList()
        println("Performing joint germline variant calling")
        ch_joint_in.view()
        // GATK_JOINT_GERMLINE_VARIANT_CALLING(  ch_joint_in , run_haplotc , run_vqsr , fasta , fai , dict, sites , sites_index , joint_id , joint_intervals , allelespecific , resources , annotation , mode , false , truthsensitivity )
    }

    if (tumor_somatic) {
        ch_tumor_in = filetest.toList()
        println("Performing tumor-only somatic variant calling")
        ch_tumor_in.view()
        // GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING(  ch_tumor_in , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    }

    if (tumor_normal_somatic) {
        ch_tumor_normal_in = filetest.map{ meta, input_file, input_index, intervals_file, which_normal ->
            filtered_normal = which_normal.unique().toList()
            [meta, input_file, input_index, intervals_file, filtered_normal]

        }.collect()
        println("Performing tumor-normal somatic variant calling")
        ch_tumor_normal_in.view()
        GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(  ch_tumor_normal_in , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    }

}
