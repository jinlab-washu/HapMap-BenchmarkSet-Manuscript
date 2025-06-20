version 1.0
# Remove PolyG artifacts from reads. These show up due to PCR so theoreticvlly they will not appear in the simulated data but to match the pipeline for aligning read data we still perfrom this step

import "Structs.wdl"

task RemovePolyG{
    input {
        Array[File] Fastq
        String FastpDocker
        Int read_length
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 1.7,
      disk_gb: 10,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }

    runtime {
      cpu: default_attr.cpu_cores
      memory: default_attr.mem_gb + " GiB"
      disks: "local-disk " + default_attr.disk_gb + " HDD"
      bootDiskSizeGb: default_attr.boot_disk_gb
      preemptible: default_attr.preemptible_tries
      maxRetries: default_attr.max_retries
      docker: FastpDocker
    }
    String basenameFoward = basename(Fastq[0], ".fq")
    String basenameReverse = basename(Fastq[1], ".fq")
    command <<<
        fastp \
            --dont_eval_duplication \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --trim_poly_g \
            --length_required ~{read_length} \
            -i ~{Fastq[0]} -I ~{Fastq[1]} \
            -o ~{basenameFoward}.trimmed.fastq -O ~{basenameReverse}.trimmed.fastq
    >>>

    output {
        Array[File] out = ["${basenameFoward}.trimmed.fastq", "${basenameReverse}.trimmed.fastq"]
    }
}