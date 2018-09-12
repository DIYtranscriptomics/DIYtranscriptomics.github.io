# to run this 'shell script' you will need to open terminal and navigate to the directory where this script resides on your computer,
# change permissions on your computer so that you can run a shell script by typing: 'chmod u+x readMapping.sh' (without the quotes) at the prompt 
# then type './readMapping.sh' (without the quotes) at the prompt.  
# this will begin the process of running each line of code in the shell script as if you had manually entered it

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o uninf_rep1 -t 4 -b 30 --single -l 250 -s 30 Control1_mergedLanes_mergedRuns.fastq.gz &> uninf_rep1.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o uninf_rep2 -t 4 -b 30 --single -l 250 -s 30 Control2_mergedLanes_mergedRuns.fastq.gz &> uninf_rep2.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o uninf_rep3 -t 4 -b 30 --single -l 250 -s 30 Control3_mergedLanes_mergedRuns.fastq.gz &> uninf_rep3.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o crypto.wt_rep1 -t 4 -b 30 --single -l 250 -s 30 WT1_mergedLanes_mergedRuns.fastq.gz &> crypto.wt_rep1.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o crypto.wt_rep2 -t 4 -b 30 --single -l 250 -s 30 WT2_mergedLanes_mergedRuns.fastq.gz &> crypto.wt_rep2.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o crypto.wt_rep3 -t 4 -b 30 --single -l 250 -s 30 WT3_mergedLanes_mergedRuns.fastq.gz &> crypto.wt_rep3.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o crypto.mut_rep1 -t 4 -b 30 --single -l 250 -s 30 Trans1_mergedLanes_mergedRuns.fastq.gz &> crypto.mut_rep1.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o crypto.mut_rep2 -t 4 -b 30 --single -l 250 -s 30 Trans2_mergedLanes_mergedRuns.fastq.gz &> crypto.mut_rep2.log

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o crypto.mut_rep3 -t 4 -b 30 --single -l 250 -s 30 Trans3_mergedLanes_mergedRuns.fastq.gz &> crypto.mut_rep3.log

echo "Finished"

