#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_VQSR } from '../../../../subworkflows/nf-core/gatk_vqsr/main'

workflow test_gatk_vqsr {
    input = [
        [ id:'joint_germline' ], // meta map
        file('gs://mhra-ngs-dev-0yzc-nextflow/testing-work/4a/c777520ac621d23f6818d5be2cb8f1/joint_germline.vcf.gz',  checkIfExists: true),
        file('gs://mhra-ngs-dev-0yzc-nextflow/testing-work/4a/c777520ac621d23f6818d5be2cb8f1/joint_germline.vcf.gz.tbi',  checkIfExists: true)
      ]

    fasta                 = file('gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta', checkIfExists: true)
    fai                   = file('gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai', checkIfExists: true)
    dict                  = file('gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict', checkIfExists: true)
    allelespecific  = false
    resources_SNP             = [
       [
           file('gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz', checkIfExists: true),
           file('gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz', checkIfExists: true),
           file('gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz', checkIfExists: true),
           file('gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf', checkIfExists: true)
        ],
        [
            file('gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi', checkIfExists: true),
            file('gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi', checkIfExists: true),
            file('gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi', checkIfExists: true),
            file('gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx', checkIfExists: true)
        ],
        [
            'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
            'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
            '1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz',
            'dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf'
        ]
    ]
    
        resources_INDEL             = [
       [
           file('gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz', checkIfExists: true),
           file('gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz', checkIfExists: true),
           file('gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz', checkIfExists: true),
           file('gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz', checkIfExists: true)
        ],
        [
            file('gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi', checkIfExists: true),
            file('gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi', checkIfExists: true),
            file('gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi', checkIfExists: true),
            file('gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi', checkIfExists: true)
        ],
        [
            'hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
            'omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
            '1000G,known=false,training=true,truth=false,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
            'dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.known_indels.vcf.gz'
        ]
    ]
    
    annotation_SNP       = ['QD', 'MQ', 'FS', 'SOR']
    annotation_INDEL       = ['QD', 'FS', 'SOR']
    truthsensitivity = '99.0'
    GATK_VQSR(input, fasta, fai, dict, allelespecific , resources_SNP , resources_INDEL , annotation_SNP , annotation_INDEL , false , truthsensitivity)
}
