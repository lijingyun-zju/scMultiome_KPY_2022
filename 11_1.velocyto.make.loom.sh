#!/bin/bash

#SBATCH -c 6                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-12:00                          # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=60000                          # Memory total in MB (for all cores)
#SBATCH -o scenic.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e scenic.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jingyun.li@childrens.harvard.edu    # Email to which notifications will be sent

#velocyto run10x -@ 4 /n/data2/bch/hemonc/ckim/JINGYUN/MultiOmic.20220429/0.data.raw/220425_A00794_0599_BHVY2YDRXY_SUB12290/count/KPY-RNA /n/data2/bch/hemonc/ckim/JINGYUN/software/refdata-gex-mm10-2020-A/genes/genes.gtf
velocyto run -@ 6 -b /n/data2/bch/hemonc/ckim/JINGYUN/MultiOmic.20220429/0.data.raw/220425_A00794_0599_BHVY2YDRXY_SUB12290/count/KPY-RNA/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o output_KPY /n/data2/bch/hemonc/ckim/JINGYUN/MultiOmic.20220429/0.data.raw/220425_A00794_0599_BHVY2YDRXY_SUB12290/count/KPY-RNA/outs/gex_possorted_bam.bam /n/data2/bch/hemonc/ckim/JINGYUN/software/refdata-gex-mm10-2020-A/genes/genes.gtf 


