version 1.0

# perform local realignment to standardize location of insertions and deltions in alinged reads

import "Structs.wdl"

workflow LocalRealignment {
  meta {
    allowNestedInputs: true
  }

  input {
        File Bam
        File BamIndex
        File Reference
        File ReferenceIndex
        File ReferenceDict
        File KnownIndels
        String Docker
    }

  # Idneify nlocations that need local realignment

  call RealignerTargetCreator as RealignerTargetCreator {
    input:
      Bam=Bam,
      BamIndex=BamIndex,
      Reference=Reference,
      ReferenceIndex=ReferenceIndex,
      ReferenceDict=ReferenceDict,
      KnownSites=KnownIndels,
      Docker=Docker
  }

  # perform local realigned on target regions prom prior step

  call IndelRealigner as IndelRealigner {
    input:
      Bam=Bam,
      BamIndex=BamIndex,
      Reference=Reference,
      ReferenceIndex=ReferenceIndex,
      ReferenceDict=ReferenceDict,
      TargetList=RealignerTargetCreator.Out,
      KnownSites=KnownIndels,
      Docker=Docker
  }


  output {
    File Out = IndelRealigner.Out
  }
}

task RealignerTargetCreator{
    input {
        File Bam
        File BamIndex
        File Reference
        File ReferenceIndex
        File ReferenceDict
        File KnownSites
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
        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
            -nt 30 \
            -known ~{KnownSites} \
            -R ~{Reference} \
            -I ~{Bam} \
            -o ~{basename}_target.list 
    >>>

    output {
        File Out = "~{basename}_target.list"
    }
}

task IndelRealigner{
    input {
        File Bam
        File BamIndex
        File Reference
        File ReferenceIndex
        File ReferenceDict
        File TargetList
        File KnownSites
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
        java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
            -R ~{Reference} \
            -I ~{Bam} \
            -targetIntervals ~{TargetList} \
            -known ~{KnownSites} \
            -o ~{basename}_hg38_realigned.bam
    >>>

    output {
        File Out = "~{basename}_hg38_realigned.bam"
    }
}