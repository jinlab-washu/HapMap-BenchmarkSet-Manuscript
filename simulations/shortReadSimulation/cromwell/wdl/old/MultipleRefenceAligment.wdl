version 1.0

import "Structs.wdl"
import "Alignment.wdl" as align
import "AlignmentTasks.wdl" as tasks

workflow HapMapFastqSim {
    meta {
        allowNestedInputs: true
    }

    input {
        Int TotalReadDepth
        File HG00438_Mother
        File HG00438_Father
        File HG002_Mother
        File HG002_Father
        File HG02257_Mother
        File HG02257_Father
        File HG02486_Mother
        File HG02486_Father
        File HG02622_Mother
        File HG02622_Father
        File HG005_Mother
        File HG005_Father

        Int ReadPairDistance
        Int ReadPairDistanceSD
        Int ReadLength
        Int OPTICAL_DUPLICATE_PIXEL_DISTANCE
        
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

        String WGSimDocker
        String FastpDocker
        String BwaDocker
        String GATK3Docker
        String GATK4Docker
        String PicardDocker
        String SamtoolsDocker
        
    }

    Sample HG00438 = object {
        Depth: 0.25,
        MotherFasta: HG00438_Mother,
        FatherFasta: HG00438_Father,
        SampleName: "HG00438"
    }

    Sample HG002 = object {
        Depth: 1,
        MotherFasta: HG002_Mother,
        FatherFasta: HG002_Father,
        SampleName: "HG002"
    }

    Sample HG02257 = object {
        Depth: 1,
        MotherFasta: HG02257_Mother,
        FatherFasta: HG02257_Father,
        SampleName: "HG02257"
    }

    Sample HG02486 = object {
        Depth: 1,
        MotherFasta: HG02486_Mother,
        FatherFasta: HG02486_Father,
        SampleName: "HG02486"
    }

    Sample HG02622 = object {
        Depth: 5,
        MotherFasta: HG02622_Mother,
        FatherFasta: HG02622_Father,
        SampleName: "HG02622"
    }

    Sample HG005 = object {
        Depth: 41.75,
        MotherFasta: HG005_Mother,
        FatherFasta: HG005_Father,
        SampleName: "HG005"
    }

    Int Runs=TotalReadDepth/50
    scatter (i in range(Runs)) {

        Array[Sample] MySamples=[HG00438, HG002, HG02257, HG02486, HG02622, HG005]
        scatter (samples in MySamples) {
            call RunWGSim as RunWGSimMother {
                input:
                    reference=samples.MotherFasta,
                    ReadDepth=samples.Depth,
                    ReadPairDistance=ReadPairDistance,
                    ReadPairDistanceSD=ReadPairDistanceSD,
                    ReadLength=ReadLength,
                    Docker=WGSimDocker
            }

            call RunWGSim as RunWGSimFather {
                input:
                    reference=samples.FatherFasta,
                    ReadDepth=samples.Depth,
                    ReadPairDistance=ReadPairDistance,
                    ReadPairDistanceSD=ReadPairDistanceSD,
                    ReadLength=ReadLength,
                    Docker=WGSimDocker
            }

            String FowardOut="${samples.SampleName}_foward.fa"
            String ReverseOut="${samples.SampleName}_reverse.fa"

            call tasks.combineFASTQ as combineFoward {
                input: 
                    Fastqs=[RunWGSimMother.fowardFasta, RunWGSimFather.fowardFasta],
                    Output=FowardOut
            }

            call tasks.combineFASTQ as combineReverse {
                input: 
                    Fastqs=[RunWGSimMother.reverseFasta, RunWGSimFather.reverseFasta],
                    Output=ReverseOut
            }
        }
        String FowardAllOut="HapMap_50x_${i}_foward.fa"
        String ReverseAllOut="HapMap_50x_${i}_reverse.fa"

        call tasks.combineFASTQ as combineAllFoward {
            input:
                Fastqs=combineFoward.combinedFile,
                Output=FowardAllOut
        }

        call tasks.combineFASTQ as combineAllReverse {
            input:
                Fastqs=combineReverse.combinedFile,
                Output=ReverseAllOut
        }

        call align.Alignment as Alignment {
            input:
                sample="HapMapMixture",
                FowardFastq=combineAllFoward.combinedFile,
                ReverseFastq=combineAllReverse.combinedFile,
                Reference=Reference,
                RefFAI=RefFAI,
                RefAMB=RefAMB,
                RefANN=RefANN,
                RefBWT=RefBWT,
                RefPAC=RefPAC,
                RefSA=RefSA,
                ReferenceDict=ReferenceDict,
                KnownSNVs=KnownSNVs,
                SNVIndex=SNVIndex,
                KnownIndels=KnownIndels,
                IndelIndex=IndelIndex,
                read_length=ReadLength,
                OPTICAL_DUPLICATE_PIXEL_DISTANCE=OPTICAL_DUPLICATE_PIXEL_DISTANCE,
                FastpDocker=FastpDocker,
                BwaDocker=BwaDocker,
                GATK3Docker=GATK3Docker,
                GATK4Docker=GATK4Docker,
                PicardDocker=PicardDocker
        }
    }

    String FinalBamName="HapMap_${TotalReadDepth}x_.bam"

    call tasks.combineBams as CombineBams {
        input:
            Bams=Alignment.Bam,
            OutName=FinalBamName,
            Docker=SamtoolsDocker
    }

    call tasks.SortAndIndex as SortAndIndex {
        input:
            Bam=CombineBams.combinedFile,
            Docker=SamtoolsDocker
    }

    output {
        File SimulatedBam = SortAndIndex.OutBam
        File SimulatedIndex = SortAndIndex.Index
    }
}

task RunWGSim {
    input {
        File reference
        Float ReadDepth
        Int ReadPairDistance
        Int ReadPairDistanceSD
        Int ReadLength
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
    String basename = basename(reference, ".fa")
    
    command <<<
        genomeLenght=3099922541
        numReads=`awk -v LOG=$genomeLenght 'BEGIN{print (~{ReadDepth}*LOG)/(~{ReadLength}*2*2)}' | awk '{ printf "%.0f\n", $1 }'`
        /opt/conda/bin/wgsim -N$numReads -1~{ReadLength} -2~{ReadLength} -r0 -e0 -R0 -X0 -d~{ReadPairDistance} -s~{ReadPairDistanceSD} ~{reference} ~{basename}_read1.fa ~{basename}_read2.fa
    >>>

    output {
        File fowardFasta="~{basename}_read1.fa"
        File reverseFasta="~{basename}_read2.fa"
    }
}