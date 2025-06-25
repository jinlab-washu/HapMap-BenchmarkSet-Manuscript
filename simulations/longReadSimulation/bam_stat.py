import os
import sys
import subprocess
import argparse




def sam_bam_opener(sam_fp):
    
    if not os.path.exists(sam_fp):
        print(f"Input file {sam_fp} does not exist")
        sys.exit(1)

    bam_flag = False
    if sam_fp.lower().endswith(".bam"):
        bam_flag = True

    if not bam_flag:
        fh = open(sam_fp)

    else:
        cmd = f"samtools view -h {sam_fp}"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        fh = p.stdout

    return fh


def main_report(sam_fp, report_fp):

    # sam_fp = "/scratch/wzhang/projects/SMaHT/long_read_sim/run_test2/sd_0001.sam"
    # report_fp = "/scratch/wzhang/projects/SMaHT/long_read_sim/run_test2/report.txt"

    report_fh = open(report_fp, "w")
    sbam_fh = sam_bam_opener(sam_fp)



    li = 0
    for l in sbam_fh:

        if isinstance(l, bytes):
            l = l.decode()

        if l.startswith("@"):
            continue

        li += 1

        l = l.strip().split()
        regular_fields = l[:11]
        tag_list = l[11:]

        tags = {}
        for t in tag_list:
            tag, ttype, val = t.split(":", 2)
            if ttype == "i":
                val = int(val)
            elif ttype == "f":
                val = float(val)
            tags[tag] = val




        qname = regular_fields[0]

        seq = regular_fields[9]
        regular_fields[9] = ""

        qual = regular_fields[10]
        regular_fields[10] = ""

        seq_freq = {}
        for i in range(len(seq)):
            letter = seq[i]
            if letter not in seq_freq:
                seq_freq[letter] = 0
            seq_freq[letter] += 1

        phred_freq = {}
        for i in range(len(qual)):
            phred = qual[i]
            if phred not in phred_freq:
                phred_freq[phred] = 0
            phred_freq[phred] += 1

        pass_num = tags.get("np", 0)

        # print(regular_fields)
        # print(tags)

        # print(seq_freq)
        # print(phred_freq)
        l = f"{qname}\t{len(seq)}\t{seq_freq}\t{phred_freq}\t{pass_num}\n"
        # print(l)
        report_fh.write(l)

        if li > 1000:
            # continue
            # break
            pass

    report_fh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='QC pipeline for long read simulation')
    parser.add_argument('input', type=str, help='Input bam')
    parser.add_argument('output', type=str, help='Output bam')

    args = parser.parse_args()

    input_fp = args.input
    output_fp = args.output


    if input_fp is not None and output_fp is not None:
        main_report(input_fp, output_fp)
    else:
        # Local test
        main_report("./p238.bam", "./p238_report.txt")











































































































































































