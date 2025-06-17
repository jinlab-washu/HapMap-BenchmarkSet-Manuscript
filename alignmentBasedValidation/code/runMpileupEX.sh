numberArguments=$#

if [ "$numberArguments" -ne 3 ] && [ "$numberArguments" -ne 4 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {bam File} {ref file} {run Name} {regions of interest (optional)}"
    exit 53
fi

outpath="/scratch1/fs1/jin810/ruttenberg/mpileup/outputs"
outdir="${outpath}/${3}"

# convert job index to chromosome
#chr=$LSB_JOBINDEX

chr=$LSB_JOBINDEX

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
echo $1
# if vcf provided subset subset vcf the chrome and run just on those variants, otherwise run just by chromosome
if [ "$numberArguments" -eq 3 ]; then
  bcftools mpileup -a FORMAT/AD,FORMAT/DP --no-BAQ -d 20000 -q 0 -Q 0 -f $2 -G $outdir/readgroup.txt -Ov -o $outfile -R chr${chr} $1
elif [ "$numberArguments" -eq 4 ]; then
  echo "test1"
  bcftools view -r chr${chr} $4 -o $outdir/truthsetByChrom/chr${chr}.vcf -O v
  echo "test2"
  echo "bcftools mpileup -a FORMAT/AD,FORMAT/DP --no-BAQ -d 20000 -q 0 -Q 0 -f ${2} -G $outdir/readgroup.txt -Ov -o $outfile -R $outdir/truthsetByChrom/chr${chr}.vcf $1"
  
  if [ -e "$1" ]; then
    echo "File exists."
  else
      echo "File does not exist."
  fi
  
  bcftools mpileup -a FORMAT/AD,FORMAT/DP --no-BAQ -d 20000 -q 0 -Q 0 -f $2 -G $outdir/readgroup.txt -Ov -o $outfile -R $outdir/truthsetByChrom/chr${chr}.vcf $1
  echo "test3"
  
fi