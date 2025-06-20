version 1.0

# Align a fastq file to an inputed reference.
# This consists of
# 1. Remove PolyG Artifacts: Remove Repeats from the end of a read. Not an issue with the simulated data but still done to match the process of alignment for read data 
# 2. BWA alignment: the actual alignment of the reads to the reference genome 
# 3. Mark Duplicates: mark and remove duplicate reads 
# 4. Local Alignment: Realign gaps to be consistent in where the gap is located (left aligned) 
# 5. Base Quality Score Recalibration: identify the accuracy of the bases called to provide to the user


import "RemovePolyGArtifacts.wdl" as polyg
import "BwaAlignment.wdl" as bwa
import "MarkDuplicates.wdl" as md
import "LocalRealignment.wdl" as lr
import "BQSR.wdl" as bqsr

workflow Alignment {
  meta {
    allowNestedInputs: true
  }

  input {
    String sample
    File FowardFastq
    File ReverseFastq

    File Reference
    File RefFAI
    File RefAMB
    File RefANN
    File RefBWT
    File RefPAC
    File RefSA
    File ReferenceDict

    File KnownSNVs
    File SNVIndex
    File KnownIndels
    File IndelIndex

    Int read_length
    Int OPTICAL_DUPLICATE_PIXEL_DISTANCE

    String FastpDocker
    String BwaDocker
    String GATK3Docker
    String GATK4Docker
    String PicardDocker
  }

  call polyg.RemovePolyG as RemovePolyG {
    input:
      Fastq=[FowardFastq, ReverseFastq],
      read_length=read_length,
      FastpDocker=FastpDocker
  }

  call bwa.BwaAlignment as BwaAlignment {
    input:
      sample=sample,
      Fastq=RemovePolyG.out,
      Reference=Reference,
      RefFAI=RefFAI,
      RefAMB=RefAMB,
      RefANN=RefANN,
      RefBWT=RefBWT,
      RefPAC=RefPAC,
      RefSA=RefSA,
      BwaDocker=BwaDocker
  }

  call md.MarkDuplicates as MarkDuplicates {
    input:
      Bam=BwaAlignment.AlignedFile,
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=OPTICAL_DUPLICATE_PIXEL_DISTANCE,
      PicardDocker=PicardDocker
  }

  call lr.LocalRealignment as LocalRealignment {
    input:
      Bam=MarkDuplicates.Out,
      BamIndex=MarkDuplicates.Index,
      Reference=Reference,
      ReferenceIndex=RefFAI,
      ReferenceDict=ReferenceDict,
      KnownIndels=KnownIndels,
      Docker=GATK3Docker
  }

  call bqsr.BQSR as BQSR {
    input:
      Bam=LocalRealignment.Out,
      Reference=Reference,
      ReferenceIndex=RefFAI,
      ReferenceDict=ReferenceDict,
      KnownSNVs=KnownSNVs,
      SNVIndex=SNVIndex,
      KnownIndels=KnownIndels,
      IndelIndex=IndelIndex,
      Docker=GATK4Docker
  }

  output {
    File Bam = BQSR.Out
  }
}
