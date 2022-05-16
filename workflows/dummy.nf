#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow DUMMY {
  main:
  println(params.chunking)
}
