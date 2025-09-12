#################################################################
# Set up Running mpileup

## author: Andrew Ruttenberg
## contact: ruttenberg.andrew@wustl.edu
#################################################################

numberArguments=$#

if [ "$numberArguments" -ne 4 ] && [ "$numberArguments" -ne 5 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {bam File} {ref file} {run Name} {output Directory} {regions of interest (optional)}"
    exit 53
fi

#Make Dirs for the outputs
outpath=$4
outdir="${outpath}/${3}"
mkdir -p $outdir

#if vcf provided, make dir to store vcf by chromosome
if [ "$numberArguments" -eq 4 ]
then
    mkdir $outdir/benchmarksetByChrom
fi

file=readgroup.txt
echo "SED"
sed "s|BAM|${1}|g" $file
sed "s|BAM|${1}|g" $file > $outdir/readgroup.txt
# run mpileup by chromosome
for i in {1..24}; do
  bash /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/runMpileup/runMpileupEX.sh $1 $2 $3 $4 $i $5
done

