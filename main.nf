#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// params for the input csv, intervals file used for splitting into chunks and the ID to be used for when samples are joined together
params.input                 = '/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/test_align_to_mutect.csv'
params.joint_intervals       = '/home/gavin_mackenzie_nibsc_org/code/Atholl/tests/Atholl/test_shortened_intervals.bed'
params.joint_id              = "results_joint"

// params deciding what Atholl runs and where to start the workflow.
params.alignment             = false
params.checkpoint            = true
params.chunking              = false
params.start_calling         = false
params.create_som_pon        = false
params.joint_germline        = false
params.vqsr                  = false
params.tumor_somatic         = false
params.tumor_normal_somatic  = false

// params determining options related to input csv. eg. Is the file fastq or an unaligned bam? If its a fastq is it paired or single end?
params.is_ubam               = false
params.paired                = true

//Params related to other settings. The sort order for sorting bam files by, whether variantrecalibrator is allelespecific or not, the truth sensitivity for variant recalibrator.
params.bwamem2_index         = ''
params.sort_order            = "coordinate"
params.allelespecific        = false
params.truthsensitivity      = 99.0

include { ATHOLL } from './workflows/Atholl'
include { DUMMY } from './workflows/dummy'

workflow {
    ATHOLL ()
    // DUMMY ()
}
