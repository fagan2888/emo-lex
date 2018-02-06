#! /bin/bash
#
# experiment.ksh: runs experiment for paper.
# based on run.ksh from  Eckart Zitzler, ETH Zurich, Feb 18, 2005
# written by William La Cava, UPenn, Jan 2018

# to be edited by the user
#selectors="ibea nsga2 spea2 lex"

variators="dtlz1_m3 dtlz2_m3 dtlz3_m3 dtlz4_m3 dtlz1_m5 dtlz1_m25 dtlz1_m50 dtlz1_m75 dtlz1_m100 dtlz2_m5 dtlz2_m25 dtlz2_m50 dtlz2_m75 dtlz2_m100 dtlz3_m5 dtlz3_m25 dtlz3_m50 dtlz3_m75 dtlz3_m100 dtlz4_m5 dtlz4_m25 dtlz4_m50 dtlz4_m75 dtlz4_m100"

if [ $# == 0 ] ; then
    selectors="lex nsga2 hype"
elif [ $# == 1 ] ; then 
    selectors=$1
elif [ $# == 2 ] ; then
    selectors=$1
    variators=$2
fi


os="linux";

export os
# load newer version of gcc 
module load gcc/5.2.0

path=${PWD}
for sel in $selectors ; do
    for var in $variators ; do
        m=${var#*_m}
        prob=${var%_*}
        
        # set generation
        if [ $prob == "dtlz2" ] ; then 
            if [ $m == 3 ] ; then
                generation=500
            elif [ $m == 5 ] ; then
                generation=1000
            elif [ $m == 25 ] ; then
                generation=1600
            elif [ $m == 50 ] ; then
                generation=2400
            elif [ $m == 75 ] ; then
                generation=3000
            elif [ $m == 100 ] ; then
                generation=4000
            fi 
        else             
            if [ $m == 3 ] ; then
                generation=500
            elif [ $m == 5 ] ; then
                generation=1000
            elif [ $m == 25 ] ; then
                generation=2000
            elif [ $m == 50 ] ; then
                generation=3000
            elif [ $m == 75 ] ; then
                generation=4000
            elif [ $m == 100 ] ; then
                generation=5000
            fi
        fi 

        echo "prob: $prob, m: $m, generation: $generation"
#            echo $var $sel $t
        out_file="$path/runs/bsub/$sel-$var-%J.out"
        err_file="$path/runs/bsub/$sel-$var-%J.err"
#            echo "out_file:$out_file"
        echo "bsub -o $out_file -e $err_file -n 1 -q moore7_normal -R "span[hosts=1]" -J "$sel"_"$var" "./compute.sh run $sel $var $generation""
        bsub -o $out_file -e $err_file -n 1 -q moore7_normal -R "span[hosts=1]" -J "$sel"_"$var" "./compute.sh run $sel $var $generation"
   done ;
done ;


