version 1.0

# Use BWA to align the simulated reads to a reference genome

import "Structs.wdl"

task BwaAlignment{
    input {
        Array[File] Fastq
        File Reference
        File RefFAI
        File RefAMB
        File RefANN
        File RefBWT
        File RefPAC
        File RefSA
        String BwaDocker
        String sample
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
      docker: BwaDocker
    }
    String basename = basename(Fastq[0], "_combined_Foward.fq")
    String method="wgsim_simulated"
    command <<<
        bwa mem -K 10000000 -t 30 -R '@RG\tID:~{sample}_~{method}\tSM:~{sample}\tPL:~{method}' ~{Reference} ~{Fastq[0]} ~{Fastq[1]} > ~{basename}_temp.bam
        samtools sort --no-PG -@ 30 -o ~{basename}_sorted.bam ~{basename}_temp.bam 
    >>>

    output {
        File AlignedFile = "~{basename}_sorted.bam"
    }
}