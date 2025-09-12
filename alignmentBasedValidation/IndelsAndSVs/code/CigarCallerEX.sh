#################################################################
# Make a sam file for each region of interest and call all cigar 
# code variants in each of the sam files

## author: Andrew Ruttenberg
## contact: ruttenberg.andrew@wustl.edu
#################################################################

numberArguments=$#

if [ "$numberArguments" -ne 7 ] && [ "$numberArguments" -ne 8 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {benchmarkset VCF} {Bam File} {Reference File} {if LRS or Not} {indels or not} {wkdir} {outname} {chrom}"
    exit 53
fi

benchmark=$1
bam=$2
ref=$3
LRS=$4
indels=$5
wkdir=$6
outfile=$7


if [[ "$LRS" != "True" && "$LRS" != "False" ]]; then
    echo "LRS is neither True nor False"
    exit 893457
fi

if [[ "$name" == *_* ]]; 
then
    echo "Run name cannot have _ in it"
    exit 53478
fi

echo "making sam files"

chrom=$8
bcftools view -r $chrom $benchmark -o "$wkdir/benchmarksets/$chrom.vcf" -O v
benchmark="$wkdir/benchmarksets/$chrom.vcf"
outsam="$wkdir/samfiles/$chrom"
mkdir $outsam


while IFS= read -r line
do
  if [[ $line == \#* ]]; then
    continue
  fi
  # Process each line (e.g., print the line)
  chrom=`echo "$line" | awk '{print $1}'`
  pos=`echo "$line" | awk '{print $2}'`
  length=`echo "$line" | awk '{print length($5)}'`
  lower=`echo "$pos 50" | awk '{sum = $1 - $2; print sum}'`
  upper=`echo "$pos 50 $length" | awk '{sum = $1 + $2 + $3; print sum}'`
  samtools view -@ 16 $bam "$chrom:$lower-$upper" > "${outsam}/${chrom}_${pos}.sam"
done < "$benchmark"

echo "Running Pileup"
/usr/bin/python3 /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/CigarCaller/CigarCaller.py -s $outsam -r $ref -l $LRS -o $outfile -i $indels
noVCF="${outfile%.vcf}"
sorted="${noVCF}_sorted.vcf"
bcftools sort $outfile -o $sorted
bgzip $sorted
tabix "${sorted}.gz"
rm -r $outsam
