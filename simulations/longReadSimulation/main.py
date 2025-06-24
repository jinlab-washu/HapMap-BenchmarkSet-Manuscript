
#################################################################
# Main pipeline script for generating large-scale PacBio 
# simulation workflow with job scheduling and resource management.
#
# Author: Wenjin Zhang
# Contact: wenjin@wustl.edu
#################################################################

import os
import re
import sys
import json
import copy
import math
import numpy
import subprocess


import matplotlib
import matplotlib.pyplot as plt




"""
The uncompressed folder to after compression ratio is around 8.5 

Compressed folder size:
10x 166G
40x 655G
200x ~3.2T


Simulation time cost

CCS time cost (30 threads)
 0.5% 0.9275 days
   2% 3.7days
  10% 18.5days
83.5% 155days

"""






samples = ["HG002", "HG005", "HG00438", "HG02257", "HG02486", "HG02622"]
phases = ["mat", "pat"]


sample_ratios = {
    "HG005":   0.835,

    "HG02622": 0.10,

    "HG002":   0.02,
    "HG02257": 0.02,
    "HG02486": 0.02,

    "HG00438": 0.005,
}




contig_counts = {
    "HG002_mat":   445,
    "HG002_pat":   610,
    "HG005_mat":   503,
    "HG005_pat":   682,
    "HG00438_mat": 258,
    "HG00438_pat": 276,
    "HG02257_mat": 292,
    "HG02257_pat": 306,
    "HG02486_mat": 283,
    "HG02486_pat": 342,
    "HG02622_mat": 292,
    "HG02622_pat": 270,
}
contig_lengths = {}

bjob_groups = ["compute-chenguangwang", "compute-epigenome", "compute-epigenome-condo", "compute-hprc"]

total_coverage = 4000



# CCS runtime estimation:
# 30 threads, 1x, 1bp in hours
wlr = 1600 / 135796793 / 10 / 3600



def get_fasta_fp(sample, phase):
    assert sample in samples
    assert phase in phases
    return f"/HPRC_Assemblies/Original_Contig_Only_Assemblies/{sample}/{sample}.{phase}ernal.f1_assembly_v2_genbank.fa"










# Just to figure out how many contigs were there
for sample in samples:
    for phase in phases:
        fa = get_fasta_fp(sample, phase)

        cmd = f"python3 /simulation/fasta_stat/cs.py {fa} > ./stat/{sample}_{phase}.chrom.size &"
        # print(cmd)



for sample in samples:
    contig_lengths[sample] = {}
    for phase in phases:
        contig_count = contig_counts[f"{sample}_{phase}"]

        # print(f"{sample}_{phase}: {contig_count}")

        stat_fp = f"./stat/{sample}_{phase}.chrom.size"

        contig_length = {}
        with open(stat_fp) as f:
            for l in f:
                l = l.strip()
                if len(l) == 0:
                    continue

                ci, contig_name, cl = l.split("\t")
                ci = int(ci)
                cl = int(cl)
                contig_length[ci] = cl
                # print(f"{sample}_{phase} contig{ci}: {cl}")

                # if f"{sample}_{phase}" == "HG00438_mat" and ci == 2:
                # print(f"{sample}_{phase} contig{ci}: {cl}")

        # print(f"{sample}_{phase}: {contig_count} == {contig_count2}")
        assert contig_count == len(contig_length)
        contig_lengths[sample][phase] = contig_length










"""
STEP 1: RUN PBSIM3 in parallel 
"""



threads = 30
memory = 30


split_coverages = {}

# new pwd
working_dir = f"/simulation/"
bjobs_folder = f"./step1_bjobs"
os.makedirs(bjobs_folder, exist_ok=True)


total_ratio = 0
for sample in samples:

    sample_ratio = sample_ratios[sample]
    split_coverages[sample] = {}

    for phase in phases:
        sample_phase_ratio = sample_ratio * 0.5
        cov = total_coverage * sample_phase_ratio
        covi = int(math.floor(cov))
        assert covi - cov < 1e-6

        # print(f"{sample}_{phase}: {sample_phase_ratio} * {total_coverage} = {covi}")

        coverages = [covi]

        if covi == 1670:
            coverages = [200] * 8 + [70]
            assert sum(coverages) == covi
        split_coverages[sample][phase] = coverages

        for split_i, split_cov in enumerate(coverages):
            split_i += 1
            # print(f"\tSplit {split_i}: {split_cov}")

            fa = get_fasta_fp(sample, phase)
            # sbatch_fp = f"{sbatch_folder}/sim_{s}_{p}.sbatch"
            # sbatch_log = f"{sbatch_folder}/sim_{s}_{p}.log"

            # wd = f"{working_dir}/qs_model_run/{s}_{p}/"
            # os.makedirs(wd, exist_ok=True)

            wd = f"{working_dir}{sample}_{phase}/split{split_i}/"
            # print(fa)





            bjob_script = f"""#!/bin/bash
#BSUB -n {threads + 1}
#BSUB -q general
#BSUB -M {memory}GB
#BSUB -W 10080
#BSUB -N 
#BSUB -G compute-epigenome
#BSUB -oo {working_dir}/step1_bjobs/{sample}_{phase}_split{split_i}.out
#BSUB -R 'select[mem>{memory}000 && tmp>20] rusage[mem={memory}000, tmp=20] span[hosts=1]'
#BSUB -a 'docker(ubuntu:22.04)' 



"""


            cmd = ""
            cmd += f"mkdir -p {wd}\n"
            cmd += f"cd {wd}\n\n"
            cmd += f'python3 /script/pbsim3p.py -wd {wd} -fasta {fa} -id_ad_prefix {sample}{phase[0].upper()}S{split_i} -cleanup Y -t {threads} -pbsim_params "--strategy wgs --method qshmm --qshmm /pbsim3/data/QSHMM-RSII.model --depth {split_cov} --pass-num 7 --accuracy-mean 1 --difference-ratio 1000:0:0 --length-min 300 --length-max 80000 --length-mean 22000 --length-sd 5000" '
            cmd += "\n\n"


            with open(f"{bjobs_folder}/{sample}_{phase}_split{split_i}.bjob", "w") as f:
                f.write(bjob_script)
                f.write(cmd)















"""
STEP 2: Index bam files by pbindex, and generate CCS
"""



threads = 30
memory = 50




working_dir = f"/simulation/"
bjobs_folder = f"./step2_bjobs"
os.makedirs(bjobs_folder, exist_ok=True)



for sample in samples:

    for phase in phases:

        coverages = split_coverages[sample][phase]
        contig_length = contig_lengths[sample][phase]
        # print(f"{sample}_{phase}: {splits}")

        for split_i, cov in enumerate(coverages):
            split_i += 1


            wd = f"{working_dir}{sample}_{phase}/split{split_i}/pbsim/"
            # print(wd)

            partitions = []
            work_load = 0
            current_partition = []


            for ci in contig_length.keys():
                cl = contig_length[ci]
                work_load_tmp = work_load + cl * cov
                estimated_runtime = work_load_tmp * wlr

                if estimated_runtime <= 12:
                    current_partition.append(ci)
                    work_load = work_load_tmp
                else:
                    partitions.append(current_partition)
                    current_partition = [ci]
                    work_load = cl * cov
            partitions.append(current_partition)

            tmp = []
            for p in partitions:
                tmp+= p
            assert tmp == list(contig_length.keys())

            # print(f"{sample}_{phase} ({split_i}) {len(partitions)} partitions, ")
            # print("\t", partitions)

            os.makedirs(f"./step2_bjobs/{sample}_{phase}", exist_ok=True)
            for pi, partition in enumerate(partitions):

                jn = sample.split("0")[-1] + f"{phase[0].upper()}{split_i}-{pi}"
                bjob_script = f"""#!/bin/bash
#BSUB -n {threads + 1}
#BSUB -q general
#BSUB -M {memory}GB
#BSUB -W 10080
#BSUB -N 
#BSUB -G compute-epigenome
#BSUB -J {jn}
#BSUB -oo {working_dir}/step2_bjobs/{sample}_{phase}/split{split_i}_{pi}.out
#BSUB -R 'select[mem>{memory}000 && tmp>20] rusage[mem={memory}000, tmp=20] span[hosts=1]'
#BSUB -a 'docker(ctomlins/pacbio_tools:latest)' 




bash ./split{split_i}_{pi}.sh


"""
                bash_script = f"#!/bin/bash\n\ncd {wd}\n\n\n"

                contig_count = contig_counts[f"{sample}_{phase}"]
                for ci in partition:
                    # print(f"{sample}_{phase} split{split_i} contig{ci}")

                    bam1_fp = f"p{ci}_0001.bam"
                    bam2_fp = f"p{ci}.bam"
                    log = f"../logs/CCS_p{ci}.log"

                    cmd1 = f"/opt/conda/bin/pbindex -j {threads} {bam1_fp} > {log} 2>&1"
                    cmd2 = f"/opt/conda/bin/ccs -j {threads} {bam1_fp} {bam2_fp} >> {log} 2>&1"

                    bash_script += f"{cmd1}\n"
                    bash_script += f"{cmd2}\n\n"

                with open(f"{bjobs_folder}/{sample}_{phase}/split{split_i}_{pi}.bjob", "w") as f:
                    f.write(bjob_script)

                with open(f"{bjobs_folder}/{sample}_{phase}/split{split_i}_{pi}.sh", "w") as f:
                    f.write(bash_script)











"""
STEP 3: Move CCS bam files

It was supposed to be merging CCS bam files, but samtools cat is causing issues

"""


# Step 3
check_sh_folder = f"./step3_1/"
os.makedirs(check_sh_folder, exist_ok=True)

link_sh_folder = f"./step3_2/"
os.makedirs(link_sh_folder, exist_ok=True)

for sample in samples:

    wd = f"/separated/{sample}"

    cmd1 = f"#!/bin/bash\n\n\n"

    cmd2 = f"#!/bin/bash\n\n\nmkdir -p {wd}/\n\n"
    cmd2 += f"cd {wd}/\n\n"

    for phase in phases:

        coverages = split_coverages[sample][phase]
        contig_length = contig_lengths[sample][phase]
        # print(f"{sample}_{phase}: {splits}")

        for split_i, cov in enumerate(coverages):
            split_i += 1

            for ci in contig_length.keys():
                bam_fp = f"{working_dir}{sample}_{phase}/split{split_i}/pbsim/p{ci}.bam"
                ccs_report_fp = f"{working_dir}{sample}_{phase}/split{split_i}/pbsim/p{ci}.ccs_report.txt"


                cmd1 += f"python3 /script/ccs_report_check.py {ccs_report_fp}\n"
                cmd2 += f"cp {bam_fp} ./{phase[0].upper()}S{split_i}P{ci}.bam\n"
                # print(cmd1)
                # print(bam_fp)

        cmd1 += "\n"
        cmd2 += "\n"


    cmd1 += "\n\n"
    cmd2 += "\n\n"
    # merged bam files are causing issues for some reason
    # cmd += f"samtools cat -@ 30 ./*bam > ../merged.bam\n\n\n"



    with open(f"{check_sh_folder}/{sample}.sh", "w") as f:
        f.write(cmd1)

    with open(f"{link_sh_folder}/{sample}.sh", "w") as f:
        f.write(cmd2)









"""
STEP 4: CCS bam file QC

"""


# Step 4
bam_stat_folder = f"./step4/"
os.makedirs(bam_stat_folder, exist_ok=True)

os.makedirs(bam_stat_folder + "out", exist_ok=True)
os.makedirs(bam_stat_folder + "err", exist_ok=True)



for sample in samples:

    for phase in phases:

        coverages = split_coverages[sample][phase]
        contig_length = contig_lengths[sample][phase]
        # print(f"{sample}_{phase}: {splits}")

        for split_i, cov in enumerate(coverages):
            split_i += 1



            wd = f"{working_dir}bam_stat/"
            wds =  f"{wd}{sample}/"
            # print(f"mkdir -p {wds}")

            threads = 1
            memory = 2

            cmd1 = ""

            for ci in sorted(map(int, contig_length.keys())):
                ci = int(ci)
                bam_fp = f"{working_dir}{sample}_{phase}/split{split_i}/pbsim/p{ci}.bam"
                bam_stat_fp = f"{wds}{phase[0].upper()}S{split_i}P{ci}.txt"
                # print(bam_stat_fp)


                cmd1 += f"python3 /script/bam_stat.py {bam_fp} {bam_stat_fp}\n\n"


                if ci%10 == 0 or ci == len(contig_length):
                    xi = ci // 10

                    if ci == len(contig_length) and ci % 10 != 0:
                        xi += 1

                    bsub = f"""#!/bin/bash
#BSUB -n {threads + 1}
#BSUB -q general
#BSUB -M {memory}GB
#BSUB -W 10080
#BSUB -G compute-epigenome
#BSUB -o {working_dir}/step4/out/{sample}{phase[0].upper()}S{split_i}X{xi}.out
#BSUB -e {working_dir}/step4/err/{sample}{phase[0].upper()}S{split_i}X{xi}.err
#BSUB -R 'select[mem>{memory}000 && tmp>20] rusage[mem={memory}000, tmp=20] span[hosts=1]'
#BSUB -a 'docker(ubuntu:22.04)' 
                    \n\n\n""" + cmd1



                    with open(f"{bam_stat_folder}/{sample}_{phase}_{split_i}_{xi}.bjob", "w") as f:
                        f.write(bsub + "\n\n\n\n")

                    cmd1 = ""










































