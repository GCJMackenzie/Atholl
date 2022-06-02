#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// params for the input csv, intervals file used for splitting into chunks and the ID to be used for when samples are joined together
params.input                 = '/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/test_align_to_mutect.csv'
params.joint_intervals       = '/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/test_shortened_intervals.bed'
params.joint_id              = "results_joint"

// params deciding what Atholl runs and where to start the workflow.
//turns on alignment step, aligning fastqs (or u_bams, see below) into aligned bam files.
params.run_alignment         = false
//if this is true, aligned, preprocessed split bam files are merged by sample so they can be saved as a checkpoint
params.save_checkpoint       = true
//when true, input files are assumed to be aligned, but not preprocessed bam files, process skips alignment step and splits input into chunks for preprocessing
params.start_at_chunking     = false
//when true, input is assumed to be aligned, preprocessed bam files, files are split into chunks, ready for variant calling or create sompon steps
params.start_at_calling      = false
//when true, run the create somatic panel of normals steps.
params.run_create_som_pon        = false
//when true, run the joint germline calling steps.
params.run_joint_germline        = false
//when true, assume input files are a vcf and tbi file, containing joint germline calls and skip straight to vqsr step.
params.start_at_vqsr                  = false
//when true, run somatic variant calling with just tumor samples.
params.run_tumor_somatic         = false
//when true, run somatic variant calling with tumor normal pairs.
params.run_tumor_normal_somatic  = false

// params determining options related to input csv. eg. Is the file fastq or an unaligned bam? If its a fastq is it paired or single end?
params.is_ubam               = false
params.paired                = true

//Params related to other settings. The sort order for sorting bam files by, whether variantrecalibrator is allelespecific or not, the truth sensitivity for variant recalibrator.
params.bwamem2_index         = ''
params.use_f1r2              = false
params.sort_order            = "coordinate"
params.allelespecific        = false
params.truthsensitivity      = 99.0

include { ATHOLL } from './workflows/Atholl'
include { DUMMY } from './workflows/dummy'

workflow {
    ATHOLL ()
    // DUMMY ()
}
