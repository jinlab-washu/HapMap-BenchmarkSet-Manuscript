#################################################################
# Set up calling cigar code variants for each chromosome

## author: Andrew Ruttenberg
## contact: ruttenberg.andrew@wustl.edu
#################################################################

numberArguments=$#

if [ "$numberArguments" -ne 6 ]
then
    echo "Incorrect number of arguments"
    echo "${numberArguments} arguments given"
    echo "Correct usage: bash BetterMpileup.sh {Indel/SV} {Bam File} {Reference File} {if LRS or Not (True/False)} {run name} {truthset} {Out Directory}"
    exit 53
fi

type=$1
bam=$2
ref=$3
LRS=$4
name=$5
truth=$6
wkdir=$7

if [[ "$LRS" != "True" && "$LRS" != "False" ]]; then
    echo "LRS is neither True nor False"
    exit 893457
fi

if [[ "$type" != "Indel" && "$type" != "SV" ]]; then
    echo "please specify if you are running this on indels or SVs"
    exit 234678
fi

if [[ "$name" == *_* ]]; 
then
    echo "Run name cannot have _ in it"
    exit 53478
fi

outdir="${wkdir}/${name}"
mkdir $outdir
mkdir "$outdir/outputVCFs"
mkdir "$outdir/samfiles" 

indel=$([ "$type" == "Indel" ] && echo True || echo False)

mkdir $outdir/truthsets
for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
    bash /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/CigarCaller/CigarChecker.sh $truth $bam $ref $LRS $indel "${outdir}" "${outdir}/outputVCFs/${name}_${chrom}.vcf" $chrom
done

