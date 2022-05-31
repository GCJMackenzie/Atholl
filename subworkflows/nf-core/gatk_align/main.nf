//
// Performs GATK best practice alignment  BWA, Picard sortsam, samtools index. For unaligned bam files: Samtools fastq and GATK mergebamalignments.
//

include { BWAMEM2_MEM                                       } from '../../../modules/bwamem2/mem/main.nf'
include { GATK4_MERGEBAMALIGNMENT                           } from '../../../modules/gatk4/mergebamalignment/main.nf'
include { PICARD_SORTSAM as PICARD_SORTSAM_ALIGNED          } from '../../../modules/picard/sortsam/main.nf'
include { PICARD_SORTSAM as PICARD_SORTSAM_UNMAPPED         } from '../../../modules/picard/sortsam/main.nf'
include { SAMTOOLS_FASTQ                                    } from '../../../modules/samtools/fastq/main.nf'
include { SAMTOOLS_INDEX                                    } from '../../../modules/samtools/index/main.nf'

workflow GATK_ALIGN {
    take:
    input          // channel: [ val(meta), [ input ], intervals ]
    fasta          // channel: /path/to/reference/fasta
    fai            // channel: /path/to/reference/fasta/index
    dict           // channel: /path/to/reference/fasta/dictionary
    bwaindex       // channel: /path/to/bwa/index/directory
    is_ubam        // channel: true/false whether input is in ubam format or not
    sort_order     // channel: which sort order to use for PICARD_SORTSAM_ALIGNED

    main:
    ch_versions = Channel.empty()

    ch_input = input.map {
        meta, reads, intervals ->
        [meta, reads]
    }

    ch_intervals = input.map {
        meta, reads, intervals ->
        [meta, intervals]
    }

    //
    //If the user inputs a unaligned bam file, it must be converted to fastq format, then aligned with bwamem2. Mergebamalignment is then used to add useful info from unaligned to the aligned bam.
    //If the user inputs a fastq file(s), the they are aligned with bwamem2 and mergebamalignment is skipped, since there is no unaligned bam file to extract info from.
    //
    if (is_ubam) {
        //
        //if input file is a ubam, convert unaligned bam to fastq format
        //
        SAMTOOLS_FASTQ( ch_input )
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
        ch_bwa_in = SAMTOOLS_FASTQ.out.fastq

        //
        //Align reads using bwamem2 mem
        //
        BWAMEM2_MEM ( ch_bwa_in, bwaindex, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        ch_mem_out = BWAMEM2_MEM.out.bam
        //

        //Bam files sorted using picard sortsam.
        //
        PICARD_SORTSAM_UNMAPPED ( ch_input, "queryname" )
        ch_versions = ch_versions.mix(PICARD_SORTSAM_UNMAPPED.out.versions)
        ch_unmapped = PICARD_SORTSAM_UNMAPPED.out.bam

        //
        //Use GATK4 Mergebamalignment to add additional info from the ubam, that was dropped by samtools fastq, back into the aligned bam file
        //
        ch_mergebam_in = ch_mem_out.combine(ch_unmapped, by: 0)
        GATK4_MERGEBAMALIGNMENT ( ch_mergebam_in, fasta, dict )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        ch_sortsam_in = BWAMEM2_MEM.out.bam
    } else {
        //
        //If input is a fastq file then align reads using bwamem2 mem
        //
        BWAMEM2_MEM ( ch_input, bwaindex, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
        ch_sortsam_in = BWAMEM2_MEM.out.bam
    }

    //
    //Bam files sorted using picard sortsam.
    //
    PICARD_SORTSAM_ALIGNED ( ch_sortsam_in, sort_order )
    ch_versions = ch_versions.mix(PICARD_SORTSAM_ALIGNED.out.versions)
    ch_samindex_in = PICARD_SORTSAM_ALIGNED.out.bam

    //
    //Index for sorted bam file made using samtools index
    //
    SAMTOOLS_INDEX (ch_samindex_in)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_bai = SAMTOOLS_INDEX.out.bai

    emit:
    versions                = ch_versions                                       // channel: [ versions.yml ]
    // intervals_out           = ch_intervals                                      // channel: [ val(meta), [intervals]]
    // bwa_mem_out             = BWAMEM2_MEM.out.bam                               // channel: [ val(meta), [ bam ] ]
    samtools_index_out      = SAMTOOLS_INDEX.out.bai                            // channel: [ val(meta), [ bai ] ]
    sortsam_out = PICARD_SORTSAM_ALIGNED.out.bam                                // channel: [ val(meta), [ bam ] ]
    // sortsam_unmapped_out    = is_ubam ? PICARD_SORTSAM_UNMAPPED.out.bam  : []   // channel: [ val(meta), [ bam ] ]
    // merge_bam_out           = is_ubam ? GATK4_MERGEBAMALIGNMENT.out.bam  : []   // channel: [ val(meta), [ bam ] ]
    // fastq_out               = is_ubam ? SAMTOOLS_FASTQ.out.fastq         : []   // channel: [ val(meta), [ fastq ] ]
}
