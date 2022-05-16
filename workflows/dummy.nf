#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow DUMMY {
  main:
  println("input file: $params.input")
  println("intervals file: $params.joint_intervals")
  println("joint_id: $params.joint_id")
  println("run alignment: $params.alignment")            
  println("re-merge aligned bams to save as a checkpoint: $params.checkpoint")           
  println("input files are aligned bams, split by intervals then reun preprocessing: $params.chunking")             
  println("input files are aligned and preprocessed bams, start one of the variant calling processes: $params.start_calling")        
  println("Run create somatic panel of normals subworkflow: $params.create_som_pon")       
  println("Run joint germline subworkflows: $params.joint_germline")       
  println("Run tumor only somatic variant calling: $params.tumor_somatic")        
  println("Run tumor normal paired somatic variant calling: $params.tumor_normal_somatic") 
  println("Input file is an unaligned bam: $params.is_ubam")           
  println("Input file is a paired fastq file: $params.paired")             
  println("How Samtools sort will sort bam files: $params.sort_order")          
  println("Variant_Recalibrator files are allele specific: $params.allelespecific")      
  println("Variant_Recalibrator will run with truth sensitivity set to: $params.truthsensitivity")    
}
