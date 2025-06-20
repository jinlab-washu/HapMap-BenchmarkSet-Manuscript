version 1.0

import "Structs.wdl"

task combineFASTQ{
    input {
        Array[File] Fastqs
        String Output
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
      docker: "anman1227/rutt_basic:vs3"
    }

    command <<<
      cat ~{sep=' ' Fastqs} > ~{Output}
    >>>

    output {
      File combinedFile = "~{Output}"
    }
}

task combineBams {
  input {
    Array[File] Bams
    String OutName
    String Docker
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
      docker: Docker
    }

    command <<<
      samtools merge -@ 30 -o ~{OutName} ~{sep=' ' Bams}
    >>>

    output {
      File combinedFile = "~{OutName}"
    }
}

task combineHG005 {
  input {
    File SimulatedBam
    File HG005Bam
    String OutName
    String Docker
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
      docker: Docker
    }

    command <<<
      samtools merge -@ 30 -o ~{OutName} ~{SimulatedBam} ~{HG005Bam}
    >>>

    output {
      File combinedFile = "~{OutName}"
    }
}

task SortAndIndex {
  input {
    File Bam
    String Docker
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
      docker: Docker
    }

    String base = basename(Bam, ".bam")
    String sorted = base + "_sorted.bam"
    String indexBam = base + "_sorted.bam.bai"

    command <<<
      samtools sort -@ 32 -o ~{sorted} ~{Bam}
      samtools index -@ 32 ~{sorted}
    >>>

    output {
      File OutBam = "~{sorted}"
      File Index = "~{indexBam}"
    }
}