#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../modules/bwamem2/index/main.nf'
include { ATHOLL        } from '../../Atholl'

workflow test_align_reads {

    input                 = [
                              [[ id:'test', single_end:false, read_group:'"@RG\\tID:Seq01p\\tSM:Seq01\\tPL:ILLUMINA"' ], // meta map
                              [
                                  file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
                                  file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)
                              ],
                              [],
                              file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true),
                              [] ]
                            ]

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    alignment             = true
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = false
    tumor_normal_somatic  = false

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = []
    joint_intervals       = []

    germline_resource     = []
    germline_resource_tbi = []
    panel_of_normals      = []
    panel_of_normals_tbi  = []

    temp_intervals        = []

    BWAMEM2_INDEX( fasta )
    bwaindex              = BWAMEM2_INDEX.out.index
    is_ubam               = false
    sort_order            = "coordinate"

    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, temp_intervals, bwaindex, is_ubam, sort_order, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_create_somatic_pon {
    input                 = [
                              [[ id:'test1' ], // meta map
                              [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
                              [file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
                              [],
                              [] ],
                              [[ id:'test2' ], // meta map
                              [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
                              [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
                              [],
                              [] ]
                            ]

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = false
    create_som_pon        = true
    joint_germline        = false
    tumor_somatic         = false
    tumor_normal_somatic  = false

    sites                 = []
    sites_index           = []

    joint_id              = "test_panel"
    joint_intervals       = file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true)

    germline_resource     = []
    germline_resource_tbi = []
    panel_of_normals      = []
    panel_of_normals_tbi  = []

    temp_intervals        = []

    bwaindex              = []
    is_ubam               = []
    sort_order            = []

    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, temp_intervals, bwaindex, is_ubam, sort_order, allelespecific, resources, annotation, mode, truthsensitivity )
}
