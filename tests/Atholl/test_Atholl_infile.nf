#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX } from '../../modules/bwamem2/index/main.nf'
include { ATHOLL        } from '../../Atholl'

workflow test_file_in {

    input                 = file('/home/AD/gmackenz/Atholl/Atholl/tests/Atholl/test_infile.csv', checkIfExists: true)

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

    ATHOLL ( input, fasta, fai, dict, alignment, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic, paired, sites, sites_index, joint_id, joint_intervals, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi, bwaindex, is_ubam, sort_order, run_haplotc, run_vqsr, allelespecific, resources, annotation, mode, truthsensitivity )
}
