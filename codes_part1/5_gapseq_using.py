# download gapseq from https://github.com/jotech/gapseq
# git clone https://github.com/jotech/gapseq
# Prediction of network candidate reactions, building of a draft model and gap filling:
#
# ./gapseq doall toy/myb71.fna
# Do the same but with a defined medium for gap filling:
#
# ./gapseq doall toy/ecoli.fna.gz dat/media/MM_glu.csv


import os
import shutil
fasta='myb71.fna'
shutil.copy(fasta,'/Users/xluhon/gapseq/toy')
os.chdir('/Users/xluhon/gapseq')
os.system(f'./gapseq doall toy/{fasta}')

# task 1 - build GEMs using gapseq for one species without GEMs



