numberArguments=$#

if [ "$numberArguments" -ne 3 ]
then
    echo "Incorrect number of arguments"
    echo "Correct usage: bash runWGS.sh {Read Depth} {Job group} {jobName}"
    exit 53
fi

ReadDepth=$1
JobGroup=$2
jobName=$3
workdir=/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SRS_Sims/AligmentCromwell
database=${workdir}/output/${jobName}/temp/Database

echo "making dir"

if [ ! -d ${workdir}/output/${jobName} ]; then
    mkdir ${workdir}/output/${jobName}
    mkdir ${workdir}/output/${jobName}/temp
    mkdir $database
    mkdir $database/cromwell-db
    mkdir $database/cromwell-db/cromwell-db

    echo "Dirs made"

    cat /storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SRS_Sims/AligmentCromwell/templates/configTemplate.config | sed "s|DATABASE|$database/cromwell-db|g" | sed "s|JOBGROUP|$JobGroup|g" > ${workdir}/output/${jobName}/temp/config.config
    cat /storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SRS_Sims/AligmentCromwell/templates/Template.json | sed "s|READDEPTH|$ReadDepth|g" > ${workdir}/output/${jobName}/temp/inputs.json

    echo "imputs made"
fi

echo "running cromwell"

java -Xmx10000g -Dconfig.file=${workdir}/output/${jobName}/temp/config.config -jar /opt/cromwell-84.jar run ${workdir}/wdl/RunHapMap.wdl -i ${workdir}/output/${jobName}/temp/inputs.json -p ${workdir}/wdl/dep.zip --metadata-output ${workdir}/cromwell-metadata-output/GATKSVPipelineSingleSample_metadata.json > ${workdir}/output/${jobName}/${jobName}.log

ret_code=$?

cat  ${workdir}/output/${jobName}/${jobName}.log | grep 'HapMapSim.FinalCombined' | cut -d '"' -f4 | sort | uniq > ${workdir}/output/${jobName}/${jobName}_output_list.txt

while read line; do
  cp $line ${workdir}/output/${jobName}
done <${workdir}/output/${jobName}/${jobName}_output_list.txt

workflow=$(grep "SingleWorkflowRunnerActor: Workflow submitted" ${workdir}/output/${jobName}/${jobName}_output.log | sed 's/\x1B\[[0-9;]*m//g' |rev | cut -f1 -d' ' |rev)
rm -r /scratch1/fs1/jin810/andrew_cache_output/HapMapSim/${workflow}
echo $workflow

exit $ret_code