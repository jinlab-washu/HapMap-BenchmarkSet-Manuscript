version 1.0

import "Structs.wdl"

# Idenify and remove duplicate reads in the aligned data

task MarkDuplicates{
  
    input {
        File Bam
        String PicardDocker
        Int OPTICAL_DUPLICATE_PIXEL_DISTANCE
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
      docker: PicardDocker
    }
    
    String basename = basename(Bam, ".bam")

    command <<<
        java -jar /miniconda3/envs/picard_env/share/picard-2.9.0-0/picard.jar MarkDuplicates \
        I=~{Bam} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{OPTICAL_DUPLICATE_PIXEL_DISTANCE} \
        O=~{basename}_marked_duplicates.bam \
        M=~{basename}_marked_dup_metrics.txt \
        TMP_DIR=./tmp

        samtools index -@ 20 ~{basename}_marked_duplicates.bam 
    >>>

    output {
        File Out = "~{basename}_marked_duplicates.bam"
        File Index = "~{basename}_marked_duplicates.bam.bai"
        File Metrics = "~{basename}_marked_dup_metrics.txt"
    }
}