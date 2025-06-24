#!/bin/bash

cd /simulation/HG002_mat/split1/pbsim/


/opt/conda/bin/pbindex -j 30 p1_0001.bam > ../logs/CCS_p1.log 2>&1
/opt/conda/bin/ccs -j 30 p1_0001.bam p1.bam >> ../logs/CCS_p1.log 2>&1

/opt/conda/bin/pbindex -j 30 p2_0001.bam > ../logs/CCS_p2.log 2>&1
/opt/conda/bin/ccs -j 30 p2_0001.bam p2.bam >> ../logs/CCS_p2.log 2>&1

/opt/conda/bin/pbindex -j 30 p3_0001.bam > ../logs/CCS_p3.log 2>&1
/opt/conda/bin/ccs -j 30 p3_0001.bam p3.bam >> ../logs/CCS_p3.log 2>&1

/opt/conda/bin/pbindex -j 30 p4_0001.bam > ../logs/CCS_p4.log 2>&1
/opt/conda/bin/ccs -j 30 p4_0001.bam p4.bam >> ../logs/CCS_p4.log 2>&1

/opt/conda/bin/pbindex -j 30 p5_0001.bam > ../logs/CCS_p5.log 2>&1
/opt/conda/bin/ccs -j 30 p5_0001.bam p5.bam >> ../logs/CCS_p5.log 2>&1

/opt/conda/bin/pbindex -j 30 p6_0001.bam > ../logs/CCS_p6.log 2>&1
/opt/conda/bin/ccs -j 30 p6_0001.bam p6.bam >> ../logs/CCS_p6.log 2>&1

/opt/conda/bin/pbindex -j 30 p7_0001.bam > ../logs/CCS_p7.log 2>&1
/opt/conda/bin/ccs -j 30 p7_0001.bam p7.bam >> ../logs/CCS_p7.log 2>&1

/opt/conda/bin/pbindex -j 30 p8_0001.bam > ../logs/CCS_p8.log 2>&1
/opt/conda/bin/ccs -j 30 p8_0001.bam p8.bam >> ../logs/CCS_p8.log 2>&1

/opt/conda/bin/pbindex -j 30 p9_0001.bam > ../logs/CCS_p9.log 2>&1
/opt/conda/bin/ccs -j 30 p9_0001.bam p9.bam >> ../logs/CCS_p9.log 2>&1

/opt/conda/bin/pbindex -j 30 p10_0001.bam > ../logs/CCS_p10.log 2>&1
/opt/conda/bin/ccs -j 30 p10_0001.bam p10.bam >> ../logs/CCS_p10.log 2>&1

/opt/conda/bin/pbindex -j 30 p11_0001.bam > ../logs/CCS_p11.log 2>&1
/opt/conda/bin/ccs -j 30 p11_0001.bam p11.bam >> ../logs/CCS_p11.log 2>&1

/opt/conda/bin/pbindex -j 30 p12_0001.bam > ../logs/CCS_p12.log 2>&1
/opt/conda/bin/ccs -j 30 p12_0001.bam p12.bam >> ../logs/CCS_p12.log 2>&1

/opt/conda/bin/pbindex -j 30 p13_0001.bam > ../logs/CCS_p13.log 2>&1
/opt/conda/bin/ccs -j 30 p13_0001.bam p13.bam >> ../logs/CCS_p13.log 2>&1

/opt/conda/bin/pbindex -j 30 p14_0001.bam > ../logs/CCS_p14.log 2>&1
/opt/conda/bin/ccs -j 30 p14_0001.bam p14.bam >> ../logs/CCS_p14.log 2>&1

/opt/conda/bin/pbindex -j 30 p15_0001.bam > ../logs/CCS_p15.log 2>&1
/opt/conda/bin/ccs -j 30 p15_0001.bam p15.bam >> ../logs/CCS_p15.log 2>&1

/opt/conda/bin/pbindex -j 30 p16_0001.bam > ../logs/CCS_p16.log 2>&1
/opt/conda/bin/ccs -j 30 p16_0001.bam p16.bam >> ../logs/CCS_p16.log 2>&1

/opt/conda/bin/pbindex -j 30 p17_0001.bam > ../logs/CCS_p17.log 2>&1
/opt/conda/bin/ccs -j 30 p17_0001.bam p17.bam >> ../logs/CCS_p17.log 2>&1

/opt/conda/bin/pbindex -j 30 p18_0001.bam > ../logs/CCS_p18.log 2>&1
/opt/conda/bin/ccs -j 30 p18_0001.bam p18.bam >> ../logs/CCS_p18.log 2>&1

