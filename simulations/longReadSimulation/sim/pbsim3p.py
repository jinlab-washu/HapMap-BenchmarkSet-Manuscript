import os
import re
import sys
import gzip
import shutil
import argparse
import subprocess
import multiprocessing











# This program achieves parallelization by dividing the input fasta file into multiple parts, and then running PBSIM3 on each part in parallel.


# Steps
# 1. Divide the input fasta file into multiple parts
# 2. Run PBSIM3 on each part in parallel
# 3. Combine the results from each part into a single output file
# 4. Clean up the intermediate files and logs



# Folder structure
"""
- working_directory
    - fasta
        - p1.fasta
        - p2.fasta
        - ...
    - pbsim
        - p1
            - p1.sam or p1.fq
    - logs
        - split.log
        - p1.log
        - p2.log
        - ...
    - output.fastq or output.sam
    
"""








pbsim_bin_fp = f"/Users/wenjin/code/wustl_research/SMaHT/pbsim3p/pbsim3/bin/bin/pbsim"
# TODO
# pbsim_bin_fp = f"pbsim"






class Utility(object):

    def __int__(self):
        pass

    @staticmethod
    def is_gzip(fp):
        fpl = fp.lower()
        return fpl.endswith(".gz") or fpl.endswith(".gzip")


    @staticmethod
    def file_opener(fp, mode="r"):
        if not os.path.exists(fp) and mode == "r":
            raise FileNotFoundError(f"File not found: {fp}")

        if Utility.is_gzip(fp):
            return gzip.open(fp, mode + "t")
        else:
            return open(fp, mode)

    @staticmethod
    def arg2bool(arg):
        if arg.lower().startswith("y"):
            return True
        else:
            return False



class DirectoryManager(object):

    def __init__(self):
        pass

    def set_wd(self, wd):
        self.wd = wd

    def prepare_wd(self):
        if os.path.exists(self.wd):
            shutil.rmtree(self.wd)

        # Set up the working directory
        os.makedirs(self.wd, exist_ok=True)
        os.makedirs(self.get_fasta_folder(), exist_ok=True)
        os.makedirs(self.get_pbsim_directory(), exist_ok=True)
        os.makedirs(self.get_logs_folder(), exist_ok=True)


    def get_fasta_folder(self):
        return os.path.join(self.wd, "fasta")

    def get_fasta_fp(self, part):
        return os.path.join(self.get_fasta_folder(), f"p{part}.fasta")

    def get_pbsim_directory(self):
        return os.path.join(self.wd, "pbsim")

    def get_pbsim_workding_directory(self, part):
        return os.path.join(self.get_pbsim_directory(), f"p{part}")

    def get_logs_folder(self):
        return os.path.join(self.wd, "logs")













util = Utility()
directory_manager = DirectoryManager()






def divide_fasta(fasta_fp):
    num_seq = 0

    fasta_fh = util.file_opener(fasta_fp)
    seq_lengths = {}

    init = True
    split_fq_fh = None
    seg_length = 0
    for line in fasta_fh:
        if line.startswith(">"):
            num_seq += 1

            split_fq_fp = directory_manager.get_fasta_fp(num_seq)
            if init:
                split_fq_fh = open(split_fq_fp, "w")
                init = False
            else:
                # TODO
                # if num_seq > 1:
                #     break

                split_fq_fh.close()
                split_fq_fh = open(split_fq_fp, "w")
                seq_lengths[num_seq-1] = seg_length

            seg_length = 0

        else:
            seg_length += len(line.strip())

        split_fq_fh.write(line)
    split_fq_fh.close()

    seq_lengths[num_seq] = seg_length

    return seq_lengths




def system_call(cmd):
    # print(cmd)
    subprocess.run(cmd, shell=True)

    return 0






def main(args):
    wd = args.wd
    fasta = args.fasta
    threads = args.t
    pbsim_params = args.pbsim_params
    cleanup = util.arg2bool(args.cleanup)
    compress = util.arg2bool(args.compress)
    id_ad_prefix = ""
    if args.id_ad_prefix:
        id_ad_prefix = args.id_ad_prefix
    id_ad_prefix = id_ad_prefix.strip()



    threads = int(threads)


    set_pass_num = 1
    pns = re.compile(r"--pass-num (\d+)").findall(pbsim_params)
    if len(pns) == 0:
        set_pass_num = 1
    if len(pns) > 1:
        raise ValueError("Multiple pass-num parameters")
    if len(pns) == 1:
        set_pass_num = int(pns[0])

    # print(f"Pass num: {set_pass_num}")



    directory_manager.set_wd(wd)
    directory_manager.prepare_wd()



    seq_lengths = divide_fasta(fasta)
    num_seq = len(seq_lengths)

    # Sort dict based on value
    seq_lengths = {k: v for k, v in sorted(seq_lengths.items(), key=lambda item: item[1], reverse=True)}
    #for k, v in seq_lengths.items():
    #    print(f"Sequence {k}: {v} bp")


    # Generate PBSIM3 commands
    pbsim_cmds = []
    for i in seq_lengths.keys():

        fasta_fp = directory_manager.get_fasta_fp(i)
        pbsim_prefix = directory_manager.get_pbsim_workding_directory(i)

        pbsim_out = os.path.join(directory_manager.get_logs_folder(), f"PBSIM_p{i}.log")



        pbsim_cmd = f"{pbsim_bin_fp} {pbsim_params} --genome {fasta_fp} --prefix {pbsim_prefix} --id-prefix {id_ad_prefix}P{i}_ > {pbsim_out} 2>&1"
        print(pbsim_cmd)
        if compress:

            fmt = "fastq"
            if set_pass_num > 1:
                fmt = "sam"

            o1 = f"{pbsim_prefix}_0001.{fmt}"
            o2 = f"{pbsim_prefix}_0001.maf"
            o3 = f"{pbsim_prefix}_0001.ref"

            gzip_files = [o2, o3]

            if fmt == "sam":
                o4 = o1[:-4] + ".bam"
                cmd = f"samtools view -bS {o1} > {o4}"
                pbsim_cmd += f" && {cmd}"
                pbsim_cmd += f" && rm {o1}"
            else:
                gzip_files.append(o1)

            for f in gzip_files:
                pbsim_cmd += f" && gzip {f}"

        # print()
        # print(pbsim_cmd)
        pbsim_cmds.append(pbsim_cmd)


    # Run PBSIM3 in parallel
    pool = multiprocessing.Pool(threads)
    r = pool.imap(system_call, pbsim_cmds)
    for r0 in r:
        pass
    pool.close()


    if cleanup:
        shutil.rmtree(directory_manager.get_fasta_folder())








if __name__ == "__main__":

    arg = argparse.ArgumentParser(description="PBSIM3P: Parallelized PBSIM3")
    arg.add_argument("-wd", help="The working directory", required=True)
    arg.add_argument("-id_ad_prefix", help="Additional ID prefix in the output file. To avoid ID collision.", required=False)
    arg.add_argument("-fasta", help="The input fasta file", required=True)
    arg.add_argument("-t", help="Number of threads", default=1)
    arg.add_argument("-pbsim_params", help="PBSIM3 parameter", required=True)

    arg.add_argument("-cleanup", help="Delete intermediate files", default="y")
    arg.add_argument("-compress", help="Compress the output file", default="y")


    arg = arg.parse_args()



    main(arg)




























































































