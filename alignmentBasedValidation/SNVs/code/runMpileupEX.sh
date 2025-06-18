#################################################################
# Running mpileup for an individual chromosome

## author: Andrew Ruttenberg
## contact: ruttenberg.andrew@wustl.edu
#################################################################

numberArguments=$#

if [ "$numberArguments" -ne 5 ] && [ "$numberArguments" -ne 6 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {bam File} {ref file} {run Name} {out Directroy} {chrom} {regions of interest (optional)}"
    exit 53
fi

outpath=$4
outdir="${outpath}/${3}"

chr=$5

if [ "$chr" -eq 23 ]; then
      chr="X"
    elif [ "$chr" -eq 24 ]; then
      chr="Y"
    elif [ "$chr" -eq 25 ]; then
      chr="M"
    else
      chr=$chr
    fi


outfile="${outdir}/${3}_chr${chr}.vcf"

# if vcf provided subset subset vcf the chrome and run just on those variants, otherwise run just by chromosome
if [ "$numberArguments" -eq 5 ]; then
  bcftools mpileup -a FORMAT/AD,FORMAT/DP --no-BAQ -d 20000 -q 0 -Q 0 -f $2 -G $outdir/readgroup.txt -Ov -o $outfile -R chr${chr} $1
elif [ "$numberArguments" -eq 6 ]; then
  bcftools view -r chr${chr} $6 -o $outdir/truthsetByChrom/chr${chr}.vcf -O v
  
  if [ -e "$1" ]; then
    echo "File exists."
  else
      echo "File does not exist."
  fi
  
  bcftools mpileup -a FORMAT/AD,FORMAT/DP --no-BAQ -d 20000 -q 0 -Q 0 -f $2 -G $outdir/readgroup.txt -Ov -o $outfile -R $outdir/truthsetByChrom/chr${chr}.vcf $1
  echo "test3"
  
fi
