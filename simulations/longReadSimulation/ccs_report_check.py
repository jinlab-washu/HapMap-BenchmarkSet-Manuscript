

#################################################################
# Check CCS (Circular Consensus Sequencing) report files to 
# validate ZMW pass filter rates are above 90% threshold.
#
# Author: Wenjin Zhang
# Contact: wenjin@wustl.edu
#################################################################

import os
import re
import sys
import gzip




if __name__ == "__main__":

    arg = sys.argv[1:]
    if len(arg) == 0:
        print("Usage: python ccs_report_check.py <CCS report file path>")
        sys.exit(1)

    ccs_fp = arg[0]
    if not os.path.exists(ccs_fp):
        print(f"Error: file not found ({ccs_fp})")
        sys.exit(1)


    ccs_report = ""
    with open(ccs_fp) as f:
        ccs_report = f.read()


    p1 = re.compile(r"ZMWs pass filters.*: (\d*) \((\d*.\d*)%\)")
    r1 = p1.findall(ccs_report)


    try:
        assert len(r1) == 1
        zmw_pass_filter = float(r1[0][1])
    except:
        print(f"ZMW pass not found: ({ccs_fp})")
        sys.exit(1)

    if zmw_pass_filter < 90:
        print(f"ZMW pass filter {zmw_pass_filter}% : ({ccs_fp})")










