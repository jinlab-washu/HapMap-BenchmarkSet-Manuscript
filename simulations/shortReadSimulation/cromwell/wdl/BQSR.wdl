version 1.0

import "Structs.wdl"

workflow BQSR {
  meta {
    allowNestedInputs: true
  }

  input {
    File Bam
    File Reference
    File ReferenceIndex
    File ReferenceDict
    File KnownSNVs
    File SNVIndex
    File KnownIndels
    File IndelIndex
    String Docker
  }


  call BaseRecalibrator as BaseRecalibrator {
    input:
      Bam=Bam,
      Reference=Reference,
      ReferenceIndex=ReferenceIndex,
      ReferenceDict=ReferenceDict,
      KnownSNVs=KnownSNVs,
      SNVIndex=SNVIndex,
      KnownIndels=KnownIndels,
      IndelIndex=IndelIndex,
      Docker=Docker
  }

  call ApplyBQSR as ApplyBQSR {
    input:
      Bam=Bam,
      Reference=Reference,
      ReferenceIndex=ReferenceIndex,
      ReferenceDict=ReferenceDict,
      BqrsFile=BaseRecalibrator.Out,
      Docker=Docker
  }


  output {
    File Out = ApplyBQSR.Out
  }
}

task BaseRecalibrator {
  input {
    File Bam
    File Reference
    File ReferenceIndex
    File ReferenceDict
    File KnownSNVs
    File SNVIndex
    File KnownIndels
    File IndelIndex
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
    
  String basename = basename(Bam, ".bam")

  command <<<
    gatk BaseRecalibrator \
      -R ~{Reference} \
      -I ~{Bam} \
      --enable-baq \
      --known-sites ~{KnownSNVs} \
      --known-sites ~{KnownIndels} \
      -O ~{basename}_BQSR.table
  >>>

  output {
    File Out = "~{basename}_BQSR.table"
  }
}


task ApplyBQSR {
  input {
    File Bam
    File Reference
    File ReferenceIndex
    File ReferenceDict
    File BqrsFile
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
    
  String basename = basename(Bam, ".bam")

  command <<<
    gatk ApplyBQSR \
      -R ~{Reference} \
      -I ~{Bam} \
      -bqsr ~{BqrsFile} \
      -O ~{basename}_recalibrated.bam
  >>>

  output {
    File Out = "~{basename}_recalibrated.bam"
  }
}