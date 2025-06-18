numberArguments=$#

if [ "$numberArguments" -ne 2 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {CigarCaller Output Dir} {output file name}"
    exit 53
fi

indir=$1
out=$2

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


outfile=$indir/$out.vcf
sorted=`echo $outfile | sed 's/.vcf/_sorted.vcf/g'`

command="bcftools concat"
cd $indir
for file in $indir/*.vcf.gz; do
    command="${command} $file" 
done
command="${command} -Ov -o ${outfile}"

$command
bcftools sort $outfile -Ov -o $sorted
bgzip $sorted
tabix "${sorted}.gz"
