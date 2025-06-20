version 1.0

import "Structs.wdl"

task CollectAlignmentSummaryMetrics {
    input {
        File Ref
        File input_bam
        String GatkDocker

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
        docker: GatkDocker
    }

    command <<<
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" CollectAlignmentSummaryMetrics \
            --REFERENCE_SEQUENCE ~{Ref} \
            --INPUT ~{input_bam} \
            --OUTPUT "CollectAlignmentSummaryMetrics.txt" \
            --TMP_DIR "temp"
    >>>

    output {
        File Metrics = "CollectAlignmentSummaryMetrics.txt"
    }
}

task CollectBaseDistributionByCycle {
    input {
        File input_bam
        String GatkDocker

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
        docker: GatkDocker
    }

    command <<<
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" CollectBaseDistributionByCycle \
            --CHART_OUTPUT "CollectBaseDistributionByCycle.pdf" \
            --INPUT ~{input_bam} \
            --OUTPUT "CollectBaseDistributionByCycle.txt" \
            --TMP_DIR "temp"
    >>>

    output {
        File Metrics = "CollectBaseDistributionByCycle.txt"
        File PDF = "CollectBaseDistributionByCycle.pdf"
    }
}

task CollectGcBiasMetrics {
    input {
        File Ref
        File input_bam
        String GatkDocker

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
        docker: GatkDocker
    }

    command <<<
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" CollectGcBiasMetrics \
            --REFERENCE_SEQUENCE ~{Ref} \
            --CHART_OUTPUT "CollectGcBiasMetrics.pdf" \
            --SUMMARY_OUTPUT "CollectGcBiasMetrics.summary.txt" \
            --INPUT ~{input_bam} \
            --OUTPUT "CollectGcBiasMetrics.txt" \
            --TMP_DIR "temp"
    >>>

    output {
        File PDF = "CollectGcBiasMetrics.pdf"
        File Summary = "CollectGcBiasMetrics.summary.txt"
        File Metrics = "CollectGcBiasMetrics.txt"
    }
}

task CollectInsertSizeMetrics {
    input {
        File input_bam
        String GatkDocker

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
        docker: GatkDocker
    }

    command <<<
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" CollectInsertSizeMetrics \
            --Histogram_FILE "CollectInsertSizeMetrics.pdf" \
            --INPUT ~{input_bam} \
            --OUTPUT "CollectInsertSizeMetrics.txt" \
            --TMP_DIR "temp"
    >>>

    output {
        File Hist = "CollectInsertSizeMetrics.pdf"
        File Metrics = "CollectInsertSizeMetrics.txt"
    }
}

task CollectWgsMetrics {
    input {
        File Ref
        File input_bam
        String GatkDocker

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
        docker: GatkDocker
    }

    command <<<
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" CollectWgsMetrics \
            --REFERENCE_SEQUENCE ~{Ref} \
            --COVERAGE_CAP 1000 \
            --INPUT ~{input_bam} \
            --OUTPUT "CollectWgsMetrics.txt" \
            --TMP_DIR "temp"
    >>>

    output {
        File Metrics = "CollectWgsMetrics.txt"
    }
}

task MeanQualityByCycle {
    input {
        File input_bam
        String GatkDocker

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
        docker: GatkDocker
    }

    command <<<
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" MeanQualityByCycle \
            --CHART_OUTPUT "MeanQualityByCycle.pdf" \
            --INPUT ~{input_bam} \
            --OUTPUT "MeanQualityByCycle.txt" \
            --TMP_DIR "temp"
    >>>

    output {
        File Metrics = "MeanQualityByCycle.txt"
        File PDF = "MeanQualityByCycle.pdf"
    }
}

task Stats {
    input {
        File input_bam
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
       samtools stats -@ 4 ~{input_bam} > "samtools.stats.txt"
    >>>

    output {
        File Metrics = "samtools.stats.txt"
    }
}

task Flagstat {
    input {
        File input_bam
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
        samtools flagstat -@ 4 ~{input_bam} > "samtools.flagstat.txt"
    >>>

    output {
        File Metrics = "samtools.flagstat.txt"
    }
}

task Idxstats {
    input {
        File input_bam
        File input_bam_index
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

    String index="~{input_bam}.bai"

    command <<<

        if [ ! -f ~{index} ]; then
            samtools index ~{input_bam}
        fi

        samtools idxstats ~{input_bam} > "samtools.idxstats.txt"
    >>>

    output {
        File Metrics = "samtools.idxstats.txt"
    }
}

