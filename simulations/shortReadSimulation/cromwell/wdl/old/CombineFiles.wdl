version 1.0

import "Structs.wdl"

workflow CombineFiles {
    meta {
        allowNestedInputs: true
    }

    input {

        ## 12 main inputs are fq files for maternal and paternal sequncing for each samples
        Array[File] HG002_Foward
        Array[File] HG002_Reverse
        Array[File] HG00438_Foward
        Array[File] HG00438_Reverse
        Array[File] HG02257_Foward
        Array[File] HG02257_Reverse
        Array[File] HG02486_Foward
        Array[File] HG02486_Reverse
        Array[File] HG02622_Foward
        Array[File] HG02622_Reverse

        ## HG005 can be optional since it doesn't need multiple trials, so the final BAM can just be inputed instead
        Array[File]? HG005_Foward
        Array[File]? HG005_Reverse
    }

    call combineFASTQ as combineHG002_Foward {
        input:
            Fastqs=HG002_Foward,
            Output="HG002_combined_Foward.fq"
    }

    call combineFASTQ as combineHG002_Reverse {
        input:
            Fastqs=HG002_Reverse,
            Output="HG002_combined_Reverse.fq"
    }

    call combineFASTQ as combineHG00438_Foward {
        input:
            Fastqs=HG00438_Foward,
            Output="HG00438_combined_Foward.fq"
    }

    call combineFASTQ as combineHG00438_Reverse {
        input:
            Fastqs=HG00438_Reverse,
            Output="HG00438_combined_Reverse.fq"
    }

    call combineFASTQ as combineHG02257_Foward {
        input:
            Fastqs=HG02257_Foward,
            Output="HG02257_combined_Foward.fq"
    }

    call combineFASTQ as combineHG02257_Reverse {
        input:
            Fastqs=HG02257_Reverse,
            Output="HG02257_combined_Reverse.fq"
    }

    call combineFASTQ as combineHG02486_Foward {
        input:
            Fastqs=HG02486_Foward,
            Output="HG02486_combined_Foward.fq"
    }

    call combineFASTQ as combineHG02486_Reverse {
        input:
            Fastqs=HG02486_Reverse,
            Output="HG02486_combined_Reverse.fq"
    }

    call combineFASTQ as combineHG02622_Foward {
        input:
            Fastqs=HG02622_Foward,
            Output="HG02622_combined_Foward.fq"
    }

    call combineFASTQ as combineHG02622_Reverse {
        input:
            Fastqs=HG02622_Reverse,
            Output="HG02622_combined_Reverse.fq"
    }

    output {
        #File out = combineHG002_Foward.combinedFile
        Array[File] HG002 = [combineHG002_Foward.combinedFile, combineHG002_Reverse.combinedFile]
        Array[File] HG00438 = [combineHG00438_Foward.combinedFile, combineHG00438_Reverse.combinedFile]
        Array[File] HG02257 = [combineHG02257_Foward.combinedFile, combineHG02257_Reverse.combinedFile]
        Array[File] HG02486 = [combineHG02486_Foward.combinedFile, combineHG02486_Reverse.combinedFile]
        Array[File] HG02622 = [combineHG02622_Foward.combinedFile, combineHG02622_Reverse.combinedFile]
    }
}

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

