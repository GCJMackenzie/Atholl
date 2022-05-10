#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../modules/bwamem2/index/main.nf'
include { ATHOLL        } from '../../Atholl'

workflow test_fastq_to_recalbam {

    input                 = file('/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/test_long_chunk_germline.csv', checkIfExists : true)

    alignment             = true
    checkpoint            = false
    chunking              = false
    start_calling         = false
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = true
    tumor_normal_somatic  = false
    paired                = true

    joint_id              = []
    joint_intervals       = file('/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/test_shortened_intervals.bed', checkIfExists: true)

    is_ubam               = false
    sort_order            = "coordinate"

    allelespecific        = []
    truthsensitivity      = []

    ATHOLL ( input, alignment, checkpoint, chunking, start_calling, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, joint_id, joint_intervals, is_ubam, sort_order, allelespecific, truthsensitivity )
}

workflow test_fastq_to_recalbam_multi {

    input                 = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_align_multi_infile.csv', checkIfExists : true)

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)

    alignment             = true
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = false
    tumor_normal_somatic  = false
    paired                = true

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = []
    joint_intervals       = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_shortened_intervals.bed', checkIfExists: true)

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

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_mutect_sompon {

    input                 = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_align_to_mutect.csv', checkIfExists : true)

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = true
    create_som_pon        = true
    joint_germline        = false
    tumor_somatic         = false
    tumor_normal_somatic  = false
    paired                = true

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = 'test_joint'
    joint_intervals       = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_shortened_intervals.bed', checkIfExists: true)

    germline_resource     = []
    germline_resource_tbi = []
    panel_of_normals      = []
    panel_of_normals_tbi  = []

    BWAMEM2_INDEX( fasta )
    bwaindex              = BWAMEM2_INDEX.out.index
    is_ubam               = false
    sort_order            = "coordinate"

    run_haplotc           = false
    run_vqsr              = false
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, joint_id, joint_intervals, panel_of_normals, panel_of_normals_tbi, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_mutect_to {

    input                 = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_align_to_mutect.csv', checkIfExists : true)

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = true
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = true
    tumor_normal_somatic  = false
    paired                = true

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = []
    joint_intervals       = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_shortened_intervals.bed', checkIfExists: true)

    germline_resource     = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    // BWAMEM2_INDEX( fasta )
    // bwaindex              = BWAMEM2_INDEX.out.index
    is_ubam               = false
    sort_order            = "coordinate"

    run_haplotc           = false
    run_vqsr              = false
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, joint_id, joint_intervals, is_ubam, sort_order, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_mutect_tn {

    input                 = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_align_to_mutect_TN.csv', checkIfExists : true)

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = true
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = false
    tumor_normal_somatic  = true
    paired                = true

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = []
    joint_intervals       = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_shortened_intervals.bed', checkIfExists: true)

    germline_resource     = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    BWAMEM2_INDEX( fasta )
    bwaindex              = BWAMEM2_INDEX.out.index
    is_ubam               = false
    sort_order            = "coordinate"

    run_haplotc           = false
    run_vqsr              = false
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_empty {

    input                 = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_align_to_mutect_TN.csv', checkIfExists : true)

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    alignment             = false
    create_som_pon        = false
    joint_germline        = false
    tumor_somatic         = false
    tumor_normal_somatic  = false
    paired                = false

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    joint_id              = []
    joint_intervals       = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_shortened_intervals.bed', checkIfExists: true)

    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)
    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)


    bwaindex              = []
    is_ubam               = false
    sort_order            = "coordinate"

    run_haplotc           = false
    run_vqsr              = false
    allelespecific        = []
    resources             = []
    annotation            = []
    mode                  = []
    truthsensitivity      = []

    ATHOLL ( input, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, joint_id, joint_intervals, panel_of_normals, panel_of_normals_tbi, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}

workflow test_joint_germ {
    input           = file('/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/test_long_chunk_germline.csv', checkIfExists : true)

    // input           = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_long_germline.csv', checkIfExists : true)

    alignment             = false
    chunking              = true
    create_som_pon        = false
    joint_germline        = true
    tumor_somatic         = false
    tumor_normal_somatic  = false
    paired                = false

    joint_id              = "joint_germline"
    joint_intervals       = file('/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/wgs_calling_regions_hg38_latest.bed', checkIfExists: true)

    // joint_intervals       = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/wgs_calling_regions_hg38_latest.bed', checkIfExists: true)

    is_ubam               = false
    sort_order            = "coordinate"

    allelespecific        = false

    // resources             = [
    //     [
    //         file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz'], checkIfExists: true),
    //         file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz'], checkIfExists: true),
    //         file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz'], checkIfExists: true),
    //         file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    //     ],
    //     [
    //         file(params.test_data['homo_sapiens']['genome']['hapmap_3_3_hg38_21_vcf_gz_tbi'], checkIfExists: true),
    //         file(params.test_data['homo_sapiens']['genome']['res_1000g_omni2_5_hg38_21_vcf_gz_tbi'], checkIfExists: true),
    //         file(params.test_data['homo_sapiens']['genome']['res_1000g_phase1_snps_hg38_21_vcf_gz_tbi'], checkIfExists: true),
    //         file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)
    //     ],
    //     [
    //         'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
    //         'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
    //         '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.hg38.vcf.gz',
    //         'dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf'
    //     ]
    // ]
    mode                  = 'SNP'
    truthsensitivity      = '99.0'

    ATHOLL ( input, alignment, chunking, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, joint_id, joint_intervals, is_ubam, sort_order, allelespecific, truthsensitivity )
}
