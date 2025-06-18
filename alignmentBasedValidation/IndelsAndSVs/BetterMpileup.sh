numberArguments=$#

if [ "$numberArguments" -ne 5 ] && [ "$numberArguments" -ne 6 ]
then
    echo "Incorrect number of arguments"
    echo "${numberArguments} arguments given"
    echo "Correct usage: bash BetterMpileup.sh {Indel/SV} {Bam File} {Reference File} {if LRS or Not} {run name} {truthset (optional)}"
    exit 53
fi

type=$1
bam=$2
ref=$3
LRS=$4
name=$5
defaultTruth="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Truthsets/${type}/${type}_*.vcf"
truth=${6:-$defaultTruth}

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

wkdir="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/${type}/CigarOutput"
outdir="${wkdir}/${name}"
mkdir $outdir
mkdir "$outdir/outputVCFs"
samdir="$outdir/samfiles" 
mkdir $samdir

indel=$([ "$type" == "Indel" ] && echo True || echo False)

if [ "$numberArguments" -eq 5 ];
then
    echo "default truthset"
    for file in $truth; do
        chrom=`basename $file | awk -F'_' '{print $2}' | awk -F '.' '{print $1}'`
        mkdir $outdir/$chrom 
        bsub -G compute-jin810 -q general -g /ruttenberg.andrew/FiftyJobs -o "${outdir}/${chrom}/mpile.log" -e "${outdir}/${chrom}/mpile.err" -R 'rusage[mem=100GB]' -M 400B -a 'docker(anman1227/rutt_basic:vs4)' bash /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/CigarCaller/CigarChecker.sh $file $bam $ref $LRS $indel "${outdir}/${chrom}" "${outdir}/outputVCFs/${chrom}.vcf"
    done
elif [ "$numberArguments" -eq 6 ];
then
    echo "truthset given"
    mkdir $outdir/truthsets
    for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
        bsub -G compute-jin810 -q general -g /ruttenberg.andrew/FiftyJobs -o "${outdir}/${chrom}_mpile.log" -e "${outdir}/${chrom}_mpile.err" -R 'rusage[mem=100GB]' -M 400GB -a 'docker(anman1227/rutt_basic:vs4)' bash /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/CigarCaller/CigarChecker.sh $truth $bam $ref $LRS $indel "${outdir}" "${outdir}/outputVCFs/${name}_${chrom}.vcf" $chrom
    done
fi
