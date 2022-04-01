process PICARD_ADDORREPLACEREADGROUPS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), val(ID), val(LB), val(PL), val(PU), val(SM)

    output:
    tuple val(meta), path("*.bam") , emit: bam
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension()
    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    id_command = ID ? "RGID=${ID}" : ''
    lb_command = LB ? "RGLB=${LB}" : ''
    pl_command = PL ? "RGPL=${PL}" : ''
    pu_command = PU ? "RGPU=${PU}" : ''
    sm_command = SM ? "RGSM=${SM}" : ''
    """
    picard \\
        -Xmx${avail_mem}g \\
        AddOrReplaceReadGroups \
        I=$input \
        O=${prefix}.bam \
        $id_command \
        $lb_command \
        $pl_command \
        $pu_command \
        $sm_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectMultipleMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
