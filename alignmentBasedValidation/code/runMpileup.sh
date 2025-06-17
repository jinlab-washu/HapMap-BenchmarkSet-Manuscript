numberArguments=$#

if [ "$numberArguments" -ne 3 ] && [ "$numberArguments" -ne 4 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {bam File} {ref file} {run Name} {regions of interest (optional)}"
    exit 53
fi

#Make Dirs for the outputs
outpath="/scratch1/fs1/jin810/ruttenberg/mpileup/outputs"
outdir="${outpath}/${3}"
logs="${outdir}/logsAndErrs"
mkdir -p $outdir
mkdir -p $logs

#if vcf provided, make dir to store vcf by chromosome
if [ "$numberArguments" -eq 4 ]
then
    mkdir $outdir/truthsetByChrom
fi

file=/storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/runMpileup/readgroup.txt
echo "SED"
sed "s|BAM|${1}|g" $file
sed "s|BAM|${1}|g" $file > $outdir/readgroup.txt
# run mpileup by chromosome
export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810/Active:/storage1/fs1/jin810/Active /storage1/fs1/smaht/Active:/storage1/fs1/smaht/Active /home/a.ruttenberg:/home/a.ruttenberg /scratch1/fs1/jin810/ruttenberg:/scratch1/fs1/jin810/ruttenberg /storage2/fs1/epigenome/Active:/storage2/fs1/epigenome/Active'
bsub -G compute-jin810 -q general -g /ruttenberg.andrew/HundredJobs -J 'mPileArray[1-24]' -e $logs/%J.%I.err -o $logs/%J.%I.out -M 500GB -R 'rusage[mem=25GB, tmp=25GB]' -a 'docker(anman1227/rutt_basic:vs5)' bash /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/runMpileup/runMpileupEX.sh $1 $2 $3 $4
