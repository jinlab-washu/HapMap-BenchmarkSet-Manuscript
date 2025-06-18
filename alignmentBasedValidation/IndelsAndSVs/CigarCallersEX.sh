numberArguments=$#

if [ "$numberArguments" -ne 7 ] && [ "$numberArguments" -ne 8 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {truthset VCF} {Bam File} {Reference File} {if LRS or Not} {indels or not} {wkdir} {outname} {chrom (optional)}"
    exit 53
fi

truth=$1
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

if [ "$numberArguments" -eq 8 ]
then 
  echo "making truthset"
  chrom=$8
  echo $truth
  bcftools view -r $chrom $truth -o "$wkdir/truthsets/$chrom.vcf" -O v
  echo "truthset made"
  truth="$wkdir/truthsets/$chrom.vcf"
  outsam="$wkdir/samfiles/$chrom"
  mkdir $outsam
elif [ "$numberArguments" -eq 7 ]
then
  outsam="$wkdir/samfiles"
fi

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
done < "$truth"

echo "Running Pileup"
/usr/bin/python3 /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/CigarCaller/CigarMaker.py -s $outsam -r $ref -l $LRS -o $outfile -i $indels
noVCF="${outfile%.vcf}"
sorted="${noVCF}_sorted.vcf"
bcftools sort $outfile -o $sorted
bgzip $sorted
tabix "${sorted}.gz"
rm -r $outsam
