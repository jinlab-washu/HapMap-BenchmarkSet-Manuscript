numberArguments=$#

if [ "$numberArguments" -ne 1 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {mPileup Out Dir}"
    exit 53
fi

indir=$1

if [ ${indir:0:1} != "/" ];
then    
    echo "must be absolute path to dir"
    exit 8
fi

if [ ${indir: -1} = "/" ];
then    
    echo "remove the / at the end of the path"
    exit 9803475
fi


LAST_SUBDIR=$(basename "$indir")
outfile=$indir/$LAST_SUBDIR.vcf
sorted=`echo $outfile | sed 's/.vcf/_sorted.vcf/g'`
clean=`echo $outfile | sed 's/.vcf/_sorted.clean.vcf/g'`

command="bcftools concat"
cd $indir
for file in $indir/*.vcf; do
    bcftools sort $file -o temp.vcf
    mv temp.vcf $file
    bgzip $file
    tabix $file.gz
    command="${command} $file.gz" 
done
command="${command} -Ov -o ${outfile}"

$command
bcftools sort $outfile -Ov -o $sorted
bgzip $sorted
tabix "${sorted}.gz"
bcftools norm -m - ${sorted}.gz -Ov | awk '$5!="<*>" {print $0}' > $clean
bgzip $clean
tabix "${clean}.gz"
