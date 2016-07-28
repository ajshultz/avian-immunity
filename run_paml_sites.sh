#!/bin/bash

##PAML SCRIPT / PIPELINE##

#get info, write to log
#move to correct dir

module load python/2.7.6-fasrc01

WORKDIR=$(pwd -P)
HOG=$1
SUBDIRNUM=`expr ${HOG##HOG} % 100`
printf -v SUBDIR "%04d" $SUBDIRNUM
PAMLDIR=$SUBDIR/$HOG
mkdir -p $PAMLDIR
cd $PAMLDIR

#set up paml runs

BASECTL=$WORKDIR/codeml.sites.ctl

CTL=${BASECTL##*/}
cp $BASECTL $HOG.$CTL
MODEL=${CTL##codeml.}
sed -i "s/SEQINPUT/$HOG.phy/" $HOG.$CTL
sed -i "s/TREEINPUT/$HOG.tree1.nwk/" $HOG.$CTL; #name of whatever output tree, make sure to count number of species
sed -i "s/OUTPUT/${MODEL%%.ctl}.out/" $HOG.$CTL

CURTIME=$(date)
echo "$CURTIME: Starting PAML runs."

#Get alignment and tree from other directory structure
cp /n/edwards_lab/Users/tsackton/ratites/analysis/paml/$PAMLDIR/$HOG.codeml.ancrec.ctl.out/$HOG.phy $HOG.phy
cp /n/edwards_lab/Users/tsackton/ratites/analysis/paml/$PAMLDIR/$HOG.codeml.ancrec.ctl.out/ancrec.out $HOG.ancrec.out
python $WORKDIR/paml_tree_extract.py $HOG.ancrec.out $HOG



echo "$(pwd -P)"
codeml $HOG.codeml.sites.ctl
CURTIME=$(date)
echo "$CURTIME: Finished Tree1 codeml.sites.ctl"
touch DONE