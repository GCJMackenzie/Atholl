#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../modules/bwamem2/index/main.nf'
include { ATHOLL        } from '../../Atholl'

workflow test_align_reads {

    input                 = [
                              [[ id:'test', single_end:false, read_group:'"@RG\\tID:Seq01p\\tSM:Seq01\\tPL:ILLUMINA"' ], // meta map
                              [file(params.test_data['homo_sapiens']['illumina']['test_umi_1_fastq_gz'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_umi_2_fastq_gz'], checkIfExists: true)],
                              [],
                              [file(params.test_data['homo_sapiens']['genome']['genome_interval_list'], checkIfExists: true)],
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

    BWAMEM2_INDEX( fasta )
    bwaindex              = BWAMEM2_INDEX.out.index
    is_ubam               = false
    sort_order            = "coordinate"

    run_haplotc           = true
    run_vqsr              = true
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
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

    bwaindex              = []
    is_ubam               = []
    sort_order            = []

    run_haplotc           = true
    run_vqsr              = true
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_tumor_only_somatic {
input = [
    [   [ id:'test' ], // meta map
        [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
        [file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true)],
        []
    ]
]

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = false
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = true
    tumor_normal_somatic  = false

    sites                 = []
    sites_index           = []

    joint_id              = []
    joint_intervals       = []

    germline_resource     = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    bwaindex              = []
    is_ubam               = []
    sort_order            = []

    run_haplotc           = true
    run_vqsr              = true
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_tumor_normal_somatic {
input = [
    [
        [ id:'test'], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true) , file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam'], checkIfExists: true)],
        [ file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true) , file(params.test_data['homo_sapiens']['illumina']['test_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true)],
        [ file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true) , file('/home/AD/gmackenz/Atholl/genome2.interval_list' , checkIfExists: true)],
        ["normal"]
    ]
]
    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = false
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = false
    tumor_normal_somatic  = true

    sites                 = []
    sites_index           = []

    joint_id              = []
    joint_intervals       = []

    germline_resource     = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    bwaindex              = []
    is_ubam               = []
    sort_order            = []

    run_haplotc           = true
    run_vqsr              = true
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_joint_skip_haplotc {
    input = [
        [
            [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_g_vcf_gz'],  checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_g_vcf_gz_tbi'],  checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true),
            []
        ],
        [
            [ id:'test2' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_g_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_g_vcf_gz_tbi'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true),
            []
        ]
    ]
    fasta           = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai             = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict            = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = false
    create_som_pon        = false
    joint_germline        = true
    tumor_somatic         = false
    tumor_normal_somatic  = false

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = "test_joint"
    joint_intervals       = file(params.test_data['homo_sapiens']['genome']['genome_21_interval_list'], checkIfExists: true)

    germline_resource     = []
    germline_resource_tbi = []
    panel_of_normals      = []
    panel_of_normals_tbi  = []

    bwaindex              = []
    is_ubam               = []
    sort_order            = []

    run_haplotc           = false
    run_vqsr              = true
    allelespecific        = false
    resources             = [
        [
            file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
        ],
        [
            file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz_tbi'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz_tbi'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz_tbi'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)
        ],
        [
            'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
            'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
            '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.hg38.vcf.gz',
            'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf.gz'
        ]
    ]
    annotation            = ['QD', 'FS', 'SOR']
    mode                  = 'SNP'
    truthsensitivity      = '99.0'

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_joint_skip_vqsr {
    input = [
        [
            [ id:'test' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
            []
        ],
        [
            [ id:'test2' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
            []
        ]
    ]
    fasta            = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai              = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict             = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    alignment             = false
    create_som_pon        = false
    joint_germline        = true
    tumor_somatic         = false
    tumor_normal_somatic  = false

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = "test_joint"
    joint_intervals       = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    germline_resource     = []
    germline_resource_tbi = []
    panel_of_normals      = []
    panel_of_normals_tbi  = []

    bwaindex              = []
    is_ubam               = []
    sort_order            = []

    run_haplotc           = true
    run_vqsr              = false
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}
