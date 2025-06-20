version 1.0

import "QC_tasks.wdl" as tasks

workflow QC {
  meta {
    allowNestedInputs: true
  }

  input {
    File Ref
    File input_bam
    File input_bam_index
    String GATKDocker
    String SamtoolsDocker

  }

  call tasks.CollectAlignmentSummaryMetrics as CollectAlignmentSummaryMetrics {
    input:
      input_bam=input_bam,
      Ref=Ref,
      GatkDocker=GATKDocker
  }

  call tasks.CollectBaseDistributionByCycle as CollectBaseDistributionByCycle {
    input:
      input_bam=input_bam,
      GatkDocker=GATKDocker
  }

  call tasks.CollectGcBiasMetrics as CollectGcBiasMetrics {
    input:
      input_bam=input_bam,
      Ref=Ref,
      GatkDocker=GATKDocker
  }

  call tasks.CollectInsertSizeMetrics as CollectInsertSizeMetrics {
    input:
      input_bam=input_bam,
      GatkDocker=GATKDocker
  }

  call tasks.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam=input_bam,
      Ref=Ref,
      GatkDocker=GATKDocker
  }

  call tasks.MeanQualityByCycle as MeanQualityByCycle {
    input:
      input_bam=input_bam,
      GatkDocker=GATKDocker
  }

  call tasks.Stats as Stats {
    input:
      input_bam=input_bam,
      Docker=SamtoolsDocker
  }

  call tasks.Flagstat as Flagstat {
    input:
      input_bam=input_bam,
      Docker=SamtoolsDocker
  }

  call tasks.Idxstats as Idxstats {
    input:
      input_bam=input_bam,
      input_bam_index=input_bam_index,
      Docker=SamtoolsDocker
  }

  output {
    File CollectAlignmentSummaryMetricsMetrics = CollectAlignmentSummaryMetrics.Metrics

    File CollectBaseDistributionByCycleMetrics = CollectBaseDistributionByCycle.Metrics
    File CollectBaseDistributionByCyclePDF = CollectBaseDistributionByCycle.PDF

    File CollectGcBiasMetricsPDF = CollectGcBiasMetrics.PDF
    File CollectGcBiasMetricsSummary = CollectGcBiasMetrics.Summary
    File CollectGcBiasMetricsMetrics = CollectGcBiasMetrics.Metrics

    File CollectInsertSizeMetricsHist = CollectInsertSizeMetrics.Hist
    File CollectInsertSizeMetricsMetrics = CollectInsertSizeMetrics.Metrics

    File CollectWgsMetricsMetrics = CollectWgsMetrics.Metrics

    File MeanQualityByCycleMetrics = MeanQualityByCycle.Metrics
    File MeanQualityByCyclePDF = MeanQualityByCycle.PDF

    File StatsMetrics = Stats.Metrics
    File FlagstatMetrics = Flagstat.Metrics
    File IdxstatsMetrics = Idxstats.Metrics
  }

}