#Code script for cellranger pipeline which is executed on Linux System. 
#Example is taken on the basis of adrenal sample 1 under FD condition  (FD_1_adrenal). 
#/../bin/cellranger: cellranger installation location.
#/../mm10: reference file location.
#/../raw/adrenal: raw fastq files location.
---------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------
#!/bin/sh
#$ -cwd 
#$ -M lifei6@sibcb.ac.cn
#$ -S /bin/bash
/../bin/cellranger count --id=FD_1_adrenal \
                   --transcriptome=/../mm10 \
                   --fastqs=/../raw/adrenal \
                   --sample=FD-1-adrenal \
                   --expect-cells=10000 \
                   --localcores=8 \
                   --localmem=64
