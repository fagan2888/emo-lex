#! /bin/bash
#
# experiment.ksh: runs experiment for paper.
# based on run.ksh from  Eckart Zitzler, ETH Zurich, Feb 18, 2005
# written by William La Cava, UPenn, Jan 2018

# to be edited by the user
#selectors="ibea nsga2 spea2 lex"

variators="dtlz1_m5 dtlz1_m25 dtlz1_m50 dtlz1_m75 dtlz1_m100 dtlz2_m5 dtlz2_m25 dtlz2_m50 dtlz2_m75 dtlz2_m100 dtlz3_m5 dtlz3_m25 dtlz3_m50 dtlz3_m75 dtlz3_m100 dtlz4_m5 dtlz4_m25 dtlz4_m50 dtlz4_m75 dtlz4_m100"

if [ $# -eq 0 ] ; then
    selectors="lex nsga2 hype epsmoea"
elif [ $# -eq 1 ] ; then 
    selectors=$1
elif [ $# -eq 2 ] ; then
    selectors=$1
    variators=$2
fi

trials=30
os="linux";

export os
# load newer version of gcc 
module load gcc/5.2.0
generation=1000
path=${PWD}
for sel in $selectors ; do
    for var in $variators ; do
        until [ $trials -eq 0 ] ; do
            m=${var#*_m}
#            echo $var $sel $t
          
            out_file="$path/runs/bsub/$sel_$var_t$t_%J.out"
            err_file="$path/runs/bsub/$sel_$var_t$t%J.err"
#            echo "out_file:$out_file"
            echo "bsub -o $out_file -e $err_file -n 4 -q moore7_normal -R "span[hosts=1]" -J "$sel"_"$var" "./compute.sh run $sel $var $generation $trials""
           # bsub -o $out_file -e $err_file -n 4 -q moore7_normal -R "span[hosts=1]" -J "$sel"_"$var" "./compute.sh run $sel $var $generation $trials"
           ((trials--))
       done ;
   done ;
done ;


