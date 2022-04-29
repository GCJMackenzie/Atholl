#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_JOINT_GERMLINE_VARIANT_CALLING } from '../../../../subworkflows/nf-core/gatk_joint_germline_variant_calling/main'

workflow test_atholl {
    input = [
        [id:'joint_chr21:1-46709983'],
        [
            '/home/AD/gmackenz/Atholl/testing/work/ca/712344c28614d98bae7e5f1438893b/test2_chr21_1_46709983.g.vcf.gz',
            '/home/AD/gmackenz/Atholl/testing/work/d9/fb23ffd906e36b18b264a051392fcf/test_chr21_1_46709983.g.vcf.gz'
        ],
        [
            '/home/AD/gmackenz/Atholl/testing/work/ca/712344c28614d98bae7e5f1438893b/test2_chr21_1_46709983.g.vcf.gz.tbi',
            '/home/AD/gmackenz/Atholl/testing/work/d9/fb23ffd906e36b18b264a051392fcf/test_chr21_1_46709983.g.vcf.gz.tbi'
        ],
        [],
        'chr21:1-46709983',
        '/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/sequence/genome.dict'
    ]
    fasta       = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai         = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict        = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)
    sites       = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    sites_index = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    GATK_JOINT_GERMLINE_VARIANT_CALLING ( input, fasta, fai, dict, sites, sites_index )
}
