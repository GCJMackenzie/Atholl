#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow DUMMY {
  main:
  println("input file: $params.input")
  println("intervals file: params.joint_intervals")
  println("joint_id: params.joint_id ")
  println("run alignment: ${params.alignment}")            
  println(params.checkpoint)           
  println(params.chunking)             
  println(params.start_calling)        
  println(params.create_som_pon)       
  println(params.joint_germline)       
  println(params.tumor_somatic)        
  println(params.tumor_normal_somatic) 
  println(params.is_ubam)           
  println(params.paired)             
  println(params.sort_order)          
  println(params.allelespecific)      
  println(params.truthsensitivity)    
}
