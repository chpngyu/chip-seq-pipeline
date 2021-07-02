#!/bin/bash
# Calcuate represneted motifs for a TF that has done in mutiple expriments in one or more biosamples.
# To run this script please create a new empty folder and run the script in the folder and
# This script includes 3 steps

### Declarations block
# please modify the values bellow to fit your environment
# provide the full path for the python scripts
ScriptRoot=/home/cpyu/temp/20210701

# one mandatory table for ChIP-Seq peaks for each expriments 
PeakTable=peaks_for_experiments
GENOME='/home/cpyu/Animals/human/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
BlackList='/home/cpyu/Animals/human/ChIP-seq/blacklist/ENCFF023CZC_sorted.bed'
### End of Declarations


#### 
# [Step 1] Obtain significant peaks for a TF in one biosample which includes multiple expriments, and
# a TF in a biosample with only a single expriment will not be merged and passed to the next step 
DIR="Biosample"
if [ -d "$DIR" ]; then
  echo "${DIR} already exists, data in the folder will be rewritten"
else
  mkdir $DIR
fi
cd $DIR
# biosample_table is a ouput file where a TF with mutiple expriments were merged
python $ScriptRoot/group_by_biosample.py -t ../$PeakTable -o batch_run -b ../biosample_table -p $ScriptRoot --genome $GENOME --blacklist $BlackList
bash batch_run

# [Step 2] Obtain significant peaks for a TF with mutiple biosamples
cd ..
DIR="Group_by_TF"
if [ -d "$DIR" ]; then
  echo "${DIR} already exists, data in the folder will be rewritten"
else
  mkdir $DIR
fi
cd $DIR
python $ScriptRoot/group_by_tf.py -t ../biosample_table -o batch_run -p $ScriptRoot --genome $GENOME --blacklist $BlackList
bash batch_run

# [Step 3] gather all inferred PWMs form each biosamples
python $ScriptRoot/gather_pwm.py ranking.meme ranking_summary meme_ranking_out
echo "Done. please check ranking.meme ranking_summary"
