#Docker: biocontainers/python3-pysam:v0.10.0ds-2-deb_cv1

numberArguments=$#
if [ "$numberArguments" -ne 2 ]
then
    echo "Incorrect number of arguments"
    echo "${numberArguments} arguments given"
    echo "Correct usage: bash phasing.sh {False neg set} {Ourdir}"
    exit 53
fi

#Step1, make list of FN variants divided by cell line
echo "Step1"

outdir=$2
MissedVars=$1
mkdir $outdir/FalseNegitiveSet
for haplo in HG002_mat HG002_pat HG00438_mat HG00438_pat HG02257_mat HG02257_pat HG02486_mat HG02486_pat HG02622_mat HG02622_pat; do
    echo $haplo
    cell="${haplo%%_*}"
    parent="${haplo#*_}"
    if [ "$parent" == "mat" ]; then
        geno="1|."
    elif [ "$parent" == "pat" ]; then
        geno=".|1"
    else
        echo "error"
        exit 372
    fi
    genosearch="${cell}_GT=${geno}"
    zcat $MissedVars | grep "#" > "$outdir/FalseNegitiveSet/${haplo}_missedVars.vcf"
    zcat $MissedVars | grep "$genosearch" >> "$outdir/FalseNegitiveSet/${haplo}_missedVars.vcf"
    bgzip "$outdir/FalseNegitiveSet/${haplo}_missedVars.vcf"
    tabix "$outdir/FalseNegitiveSet/${haplo}_missedVars.vcf.gz"
done

#Steo2, use dipcall to rescue variants
echo "Step2"

MissedVarsDir=$outdir/FalseNegitiveSet
graph="/storage2/fs1/epigenome/Active/shared_smaht/SMaHT_MCGB_Graphs/new-v2.9.7/9188e8/9188e8.vcf.gz"
ref="/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
dipDir="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/Dipcall/GRCh38"
for haplo in HG002_mat HG002_pat HG00438_mat HG00438_pat HG02257_mat HG02257_pat HG02486_mat HG02486_pat HG02622_mat HG02622_pat; do
    echo $haplo
    mkdir "${outdir}/${haplo}"
    dipSomatic=`echo $haplo | tr '_' '-'`
    for germ in HG005_mat HG005_pat; do
        dipGerm=`echo $germ | tr '_' '-'`
        dipcall="${dipDir}/GRCh38-${dipGerm}ernal-${dipSomatic}ernal.dip.vcf.gz"
        /usr/bin/python3 /storage1/fs1/jin810/Active/testing/Ruttenberg/Code/FinalSMAHTCode/RegeionReconstruction/RegeionReconstruction.py \
        -m "${MissedVarsDir}/${haplo}_missedVars.vcf.gz" \
        -g $graph \
        -d $dipcall \
        -f $ref \
        --GraphSample $haplo --DipSample syndip \
        --distance 50 \
        -o "$outdir/${haplo}/${germ}_${haplo}"
    done
done

#step3: combine vcell liones
echo "Step3"

for haplo in HG002_mat HG002_pat HG00438_mat HG00438_pat HG02257_mat HG02257_pat HG02486_mat HG02486_pat HG02622_mat HG02622_pat; do
    for germ in HG005_mat HG005_pat; do
        bgzip "$outdir/${haplo}/${germ}_${haplo}/tp.vcf"
        gzip -dk "$outdir/${haplo}/${germ}_${haplo}/tp.vcf.gz"
        tabix ${outdir}/${haplo}/${germ}_${haplo}/tp.vcf.gz
    done
    bcftools merge -m none --force-samples ${outdir}/${haplo}/HG005_mat_${haplo}/tp.vcf.gz ${outdir}/${haplo}/HG005_pat_${haplo}/tp.vcf.gz -o ${outdir}/${haplo}/${haplo}_rescused.vcf -O v
    bgzip ${outdir}/${haplo}/${haplo}_rescused.vcf
    gzip -dk ${outdir}/${haplo}/${haplo}_rescused.vcf.gz
    tabix ${outdir}/${haplo}/${haplo}_rescused.vcf.gz
done
bcftools merge -m none --force-samples ${outdir}/HG002_mat/HG002_mat_rescused.vcf.gz ${outdir}/HG002_pat/HG002_pat_rescused.vcf.gz ${outdir}/HG00438_mat/HG00438_mat_rescused.vcf.gz ${outdir}/HG00438_pat/HG00438_pat_rescused.vcf.gz ${outdir}/HG02257_mat/HG02257_mat_rescused.vcf.gz ${outdir}/HG02257_pat/HG02257_pat_rescused.vcf.gz ${outdir}/HG02486_mat/HG02486_mat_rescused.vcf.gz ${outdir}/HG02486_pat/HG02486_pat_rescused.vcf.gz ${outdir}/HG02622_mat/HG02622_mat_rescused.vcf.gz ${outdir}/HG02622_pat/HG02622_pat_rescused.vcf.gz -o ${outdir}/rescused.vcf -O v
bgzip ${outdir}/rescused.vcf
gzip -dk ${outdir}/rescused.vcf.gz
tabix ${outdir}/rescused.vcf.gz