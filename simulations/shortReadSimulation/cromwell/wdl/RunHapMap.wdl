version 1.0

import "Structs.wdl"
import "Alignment.wdl" as align
import "AlignmentTasks.wdl" as tasks

workflow HapMapSim {
    meta {
        allowNestedInputs: true
    }

    input {

        Int TotalReadDepth #Desired Read depth for the simulated bam file. Must be multiple of 500

        # Fasta Files of somatic cell lines
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

        # Bam file of HG005 Simulated data
        File HG005_Bam
        Int HG005_Depth

        # Paramters for WGSIM simuation process
        Int ReadPairDistance
        Int ReadPairDistanceSD
        Int ReadLength
        Int OPTICAL_DUPLICATE_PIXEL_DISTANCE
        Float ErrorRate = 0
        
        # Refefence Fasta Files
        File Reference
        File RefFAI
        File RefAMB
        File RefANN
        File RefBWT
        File RefPAC
        File RefSA
        File ReferenceDict

        # SNV and Indel List used for BQSR
        File KnownSNVs
        File SNVIndex
        File KnownIndels
        File IndelIndex

        # List ofDockers for the pipeline
        String WGSimDocker
        String FastpDocker
        String BwaDocker
        String GATK3Docker
        String GATK4Docker
        String PicardDocker
        String SamtoolsDocker
        
    }

    # Consolidate all info for each somatic cell line into an sample class
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
    # Calculated number of HG005 Simulated data needed to be in simuated data to reach desired read depth 
    # then combine the HG005 Bam file with itself that many timre
    Int num_fives = TotalReadDepth/HG005_Depth

    if (num_fives>1) {
        scatter (j in range(num_fives)) {
            File five=HG005_Bam
        }

        call tasks.combineBams as CombineHG005 {
            input:
                Bams=five,
                OutName="HG005_combined.bam",
                Docker=SamtoolsDocker
        }
    }

    File fives = select_first([CombineHG005.combinedFile, HG005_Bam])

    # In parallel generate somatic cell lines to 50x, which then can be combined to the desired read depth
    Int Runs=TotalReadDepth/50
    scatter (i in range(Runs)) {

        #For each somatic cell line in parallel use WGSIM to generate reads at desired read depth, and aling them to the inputed refence 
        Array[Sample] MySamples=[HG00438, HG002, HG02257, HG02486, HG02622]
        scatter (samples in MySamples) {
            #Simuated from both maternal and paternal haplotypes, will be combined prior to Alignment
            call RunWGSim as RunWGSimMother {
                input:
                    reference=samples.MotherFasta,
                    ReadDepth=samples.Depth,
                    ReadPairDistance=ReadPairDistance,
                    ReadPairDistanceSD=ReadPairDistanceSD,
                    ReadLength=ReadLength,
                    ErrorRate=ErrorRate,
                    Docker=WGSimDocker
            }

            call RunWGSim as RunWGSimFather {
                input:
                    reference=samples.FatherFasta,
                    ReadDepth=samples.Depth,
                    ReadPairDistance=ReadPairDistance,
                    ReadPairDistanceSD=ReadPairDistanceSD,
                    ReadLength=ReadLength,
                    ErrorRate=ErrorRate,
                    Docker=WGSimDocker
            }

            String FowardOut="${samples.SampleName}_foward.fa"
            String ReverseOut="${samples.SampleName}_reverse.fa"

            # Combine simuated reads from maternal and paternal haplotype
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
            # Align short read sequencing to the reference fasta
            call align.Alignment as Alignment {
                input:
                    sample=samples.SampleName,
                    FowardFastq=combineFoward.combinedFile,
                    ReverseFastq=combineReverse.combinedFile,
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
        #Combine alifned bam files across somatic cell lines
        call tasks.combineBams as CombineBams {
            input:
                Bams=Alignment.Bam,
                OutName="TestSim.bam",
                Docker=SamtoolsDocker
        }
    }
    #Combined 50x somatic cell line bam files 
    call tasks.combineBams as CombineSims {
        input:
            Bams=CombineBams.combinedFile,
            OutName="Hapmap_~{TotalReadDepth}x_simulation.bam",
            Docker=SamtoolsDocker
    }
    #combine somatic cell lines with HG005 bam file
    String FinalName = "${TotalReadDepth}x_Hapmap_simulation.bam"

    call tasks.combineBams as FinalHapMapBam {
        input:
            Bams=[CombineSims.combinedFile, fives],
            OutName=FinalName,
            Docker=SamtoolsDocker
    }

    # Sort and index and file bam file
    call tasks.SortAndIndex as SortAndIndex {
        input:
            Bam=FinalHapMapBam.combinedFile,
            Docker=SamtoolsDocker
    }

    output {
        File FinalCombinedBam = SortAndIndex.OutBam
        File FinalCombinedIndex = SortAndIndex.Index
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
        Float ErrorRate
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
    # Calculate the number of reads needed to be simulated given the lenght of the read and the lenght of the genome to achive the desired read depth, and then perform the simulation
    command <<<
        genomeLenght=3099922541
        numReads=`awk -v LOG=$genomeLenght 'BEGIN{print (~{ReadDepth}*LOG)/(~{ReadLength}*2*2)}' | awk '{ printf "%.0f\n", $1 }'`
    /opt/conda/bin/wgsim -N$numReads -1~{ReadLength} -2~{ReadLength} -r0 -e~{ErrorRate} -R0 -X0 -d~{ReadPairDistance} -s~{ReadPairDistanceSD} ~{reference} ~{basename}_read1.fa ~{basename}_read2.fa
    >>>

    output {
        File fowardFasta="~{basename}_read1.fa"
        File reverseFasta="~{basename}_read2.fa"
    }


}