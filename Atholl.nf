#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK_ALIGN } from './subworkflows/nf-core/gatk_align/main'
include { SAMTOOLS_CHUNK } from './subworkflows/nf-core/samtools_chunking/main'
include { GATK_PREPROCESS } from './subworkflows/nf-core/gatk_preprocess/main'
include { GATK_MUTECT2_CALLING } from './subworkflows/nf-core/gatk_mutect2_calling/main'
include { GATK_CREATE_SOMATIC_PON } from './subworkflows/nf-core/gatk_create_somatic_pon/main'
include { GATK_HAPLOTYPECALLING} from './subworkflows/nf-core/gatk_haplotypecalling/main'
include { GATK_JOINT_GERMLINE_VARIANT_CALLING } from './subworkflows/nf-core/gatk_joint_germline_variant_calling/main'
include { GATK_VQSR} from './subworkflows/nf-core/gatk_vqsr/main'
include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING } from './subworkflows/nf-core/gatk_tumor_only_somatic_variant_calling/main'
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from './subworkflows/nf-core/gatk_tumor_normal_somatic_variant_calling/main'
include { GATK4_MERGEVCFS } from './modules/gatk4/mergevcfs/main'
include { BWAMEM2_INDEX } from './modules/bwamem2/index/main'

workflow ATHOLL {

//    fasta                 = params.genomes.'GATK.GRCh38'.fasta
//    fai                   = params.genomes.'GATK.GRCh38'.fasta_fai
//    dict                  = params.genomes.'GATK.GRCh38'.dict

//    bwaindex              = params.genomes.'GATK.GRCh38'.bwa

//    germline_resource     = params.genomes.'GATK.GRCh38'.germline_resource
//    germline_resource_tbi = params.genomes.'GATK.GRCh38'.germline_resource_tbi

//    sites                 = params.genomes.'GATK.GRCh38'.dbsnp
//    sites_index           = params.genomes.'GATK.GRCh38'.dbsnp_tbi

    fasta                 = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)
    fai                   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)
    dict                  = file(params.test_data['homo_sapiens']['genome']['genome_21_dict'], checkIfExists: true)

    BWAMEM2_INDEX( fasta )
    bwaindex = BWAMEM2_INDEX.out.index
    //bwaindex              = Channel.fromPath('/home/AD/gmackenz/Atholl/bwamem2/genome.fasta.{amb,ann,bwt.2bit.64,pac,0123}').collect()

    germline_resource     = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz'], checkIfExists: true)
    germline_resource_tbi = file(params.test_data['homo_sapiens']['genome']['gnomad_r2_1_1_21_vcf_gz_tbi'], checkIfExists: true)

    sites                 = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz'], checkIfExists: true)
    sites_index           = file(params.test_data['homo_sapiens']['genome']['dbsnp_138_hg38_21_vcf_gz_tbi'], checkIfExists: true)

    panel_of_normals      = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz'], checkIfExists: true)
    panel_of_normals_tbi  = file(params.test_data['homo_sapiens']['genome']['mills_and_1000g_indels_21_vcf_gz_tbi'], checkIfExists: true)

    take:

    //universal args: these args are used by every subworkflow, input arg is used for passing in sample specific data, some of the array entries included are
    //not necessary for everysubworkflow, e.g which_norm, these are only passed in to the subworkflows that require them.
    input

    //entry params: used to control which subworkflow(s) are run
    alignment
    create_som_pon
    joint_germline
    tumor_somatic
    tumor_normal_somatic
    paired

    //shared args: these args are used by two or more subworkflows, but are not universal, which subworkflows use them are noted.

    joint_id              // channel: joint id for gendbs and pons    joint germline + create_som_pon
    joint_intervals       // channel: joint intervals file            joint germline + create_som_pon

    // aligner args: args exclusive to align and preprocess subworkflow
    is_ubam               // channel: true/false whether input is in ubam format or not
    sort_order            // channel: which sort order to use for PICARD_SORTSAM_DUPLICATESMARKED

    // joint germline args: args exclusive to joint germline subworkflow
    allelespecific        // channel: true/false run allelespecific mode of vqsr modules
    resources             // channel: [[resource, vcfs, forvariantrecal], [resource, tbis, forvariantrecal], [resource, labels, forvariantrecal]]
    annotation            // channel: [annotations, to, use, for, variantrecal, filtering]
    mode                  // channel: which mode to run variantrecal: SNP/INDEL/BOTH
    truthsensitivity      // channel: 0-100.0 truthsensitivity cutoff for applyvqsr

    main:
        filetest = extract_samples(input, alignment, paired, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic)

    if (alignment) {
        println("The aligner is running")
        GATK_ALIGN( filetest , fasta , fai , dict , bwaindex , is_ubam , sort_order )
        ch_chunk_in = GATK_ALIGN.out.sortsam_out.combine(GATK_ALIGN.out.samtools_index_out, by: 0)
        SAMTOOLS_CHUNK(ch_chunk_in, joint_intervals)
        GATK_PREPROCESS( SAMTOOLS_CHUNK.out.ch_format_out , fasta , fai , dict , sort_order, sites, sites_index )
    }

    if (create_som_pon) {
        ch_mutect2_sub_in = GATK_PREPROCESS.out.applybqsr_out.combine(GATK_PREPROCESS.out.samtools_index_out, by: 0).combine(GATK_PREPROCESS.out.ch_intervals_out, by: 0).map{meta, bam, bai, intervals -> [meta, bam, bai, intervals, [] ]}
        GATK_MUTECT2_CALLING(ch_mutect2_sub_in, false, false, true, fasta, fai, dict, [], [], [], [])
        ch_som_pon_vcf =  GATK_MUTECT2_CALLING.out.mutect2_vcf.collect{it[1]}.toList()
        ch__som_pon_index =  GATK_MUTECT2_CALLING.out.mutect2_tbi.collect{it[1]}.toList()
        ch_som_pon_in = Channel.of([[ id:joint_id ]]).combine(ch_som_pon_vcf).combine(ch__som_pon_index).combine([joint_intervals]).combine(['']).combine([dict])
        GATK_CREATE_SOMATIC_PON(  ch_som_pon_in , fasta , fai , dict , joint_id, joint_intervals )
    }

    if (tumor_somatic) {
        ch_mutect2_sub_in = GATK_PREPROCESS.out.applybqsr_out.combine(GATK_PREPROCESS.out.samtools_index_out, by: 0).combine(GATK_PREPROCESS.out.ch_intervals_out, by: 0).map{meta, bam, bai, intervals -> [meta, bam, bai, intervals, [] ]}
        GATK_MUTECT2_CALLING(ch_mutect2_sub_in, false, true, false, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
        ch_tumor_only_in = GATK_PREPROCESS.out.applybqsr_out.combine(GATK_PREPROCESS.out.samtools_index_out, by: 0).combine(GATK_MUTECT2_CALLING.out.mutect2_vcf, by: 0).combine(GATK_MUTECT2_CALLING.out.mutect2_tbi, by: 0).combine(GATK_MUTECT2_CALLING.out.mutect2_stats, by: 0).combine(GATK_MUTECT2_CALLING.out.mutect2_intervals, by: 0)
        println("Performing tumor-only somatic variant calling")
        GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING(  ch_tumor_only_in , fasta , fai , dict , germline_resource , germline_resource_tbi)
    }

    if (tumor_normal_somatic) {
        if (alignment){ ch_mutect2_sub_in = GATK_PREPROCESS.out.applybqsr_out.combine(GATK_PREPROCESS.out.samtools_index_out, by: 0).combine(GATK_PREPROCESS.out.ch_intervals_out, by: 0)
        .map{meta, bam, bai, intervals ->
            [meta, meta.sample, bam, bai, intervals, meta.norms ]}.groupTuple(by:[1,4])
        .map{meta, sample_id, bam, bai, intervals, which_norms ->
            which_norm = which_norms.unique().toList()
            [meta, sample_id, bam, bai, intervals, which_norm]
            }
        GATK_MUTECT2_CALLING(ch_mutect2_sub_in, true, false, false, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
        // GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(  filetest , fasta , fai , dict , germline_resource , germline_resource_tbi , panel_of_normals , panel_of_normals_tbi )
    }}

    if (joint_germline) {
        println("Performing joint germline variant calling")

        ch_haplotc_sub_in = GATK_PREPROCESS.out.applybqsr_out.combine(GATK_PREPROCESS.out.applybqsr_index_out, by: 0).combine(GATK_PREPROCESS.out.ch_intervals_out, by: 0)
        GATK_HAPLOTYPECALLING(ch_haplotc_sub_in, fasta, fai, dict, sites, sites_index)

        ch_haplo_out = GATK_HAPLOTYPECALLING.out.haplotc_vcf.combine(GATK_HAPLOTYPECALLING.out.haplotc_index, by: 0).combine(GATK_HAPLOTYPECALLING.out.haplotc_interval_out, by: 0).groupTuple(by: 3).map{meta, vcf, tbi, intervals ->
            def inter_meta = [:]
            inter_meta.id = "joint_$intervals"
            [ meta, vcf, tbi, [], intervals ]}
        ch_joint_germ_in = ch_haplo_out.combine([dict])
        ch_joint_germ_in.view()
        GATK_JOINT_GERMLINE_VARIANT_CALLING(  ch_joint_germ_in, fasta, fai, dict, sites, sites_index )

        ch_merge_vcf =  GATK_JOINT_GERMLINE_VARIANT_CALLING.out.genotype_vcf.collect{it[1]}.toList()
        ch_merge_in = Channel.of([[ id:joint_id ]]).combine(ch_merge_vcf)
        ch_merge_in.view()

        GATK4_MERGEVCFS(ch_merge_in, dict, true)

        ch_vqsr_in = GATK4_MERGEVCFS.out.vcf.combine(GATK4_MERGEVCFS.out.tbi, by: 0)
        ch_vqsr_in.view()
        GATK_VQSR(ch_vqsr_in, fasta, fai, dict, allelespecific , resources , annotation , mode , false , truthsensitivity)

    }

    else {
        println(fasta)
        println(fai)
        println(dict)
        println(sites)
        println(sites_index)
        println(bwaindex)
    }

}

def extract_samples(csv_file, alignment, paired, create_som_pon, joint_germline, tumor_somatic, tumor_normal_somatic) {
    firststep = Channel.from(csv_file).splitCsv(header: true).map{ row ->
        def meta = [:]
        meta.id = row.EntryID
        meta.sample = row.SampleID
        if (alignment) {
            meta.single_end = "$row.single_end"
            meta.rgID = "$row.rgID"
            meta.rgLB = "$row.rgLB"
            meta.rgPL = "$row.rgPL"
            meta.rgPU = "$row.rgPU"
            meta.rgSM = "$row.rgSM"
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
        def the_meta = meta.clone()
        def input_files = input_file
        def input_indexes = input_index
        def intervals = intervals_file
        the_meta.norms = which_norm
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
