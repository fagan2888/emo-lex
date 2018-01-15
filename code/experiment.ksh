#! /usr/bin/ksh
#
# experiment.ksh: runs experiment for paper.
# based on run.ksh from  Eckart Zitzler, ETH Zurich, Feb 18, 2005
# written by William La Cava, UPenn, Jan 2018

# to be edited by the user
#selectors="ibea nsga2 spea2 lex"
selectors="lex"
#variators="dtlz2 knapsack zdt6"
variators="dtlz1_m25 dtlz1_m50 dtlz1_m75 dtlz1_m100 dtlz2_m25 dtlz2_m50 dtlz2_m75 dtlz2_m100 dtlz3_m25 dtlz3_m50 dtlz3_m75 dtlz3_m100 dtlz4_m25 dtlz4_m50 dtlz4_m75 dtlz4_m100"

os="linux";

export os


if [[ $1 != "compare" ]] ; then
    for sel in $selectors ; do
	for var in $variators ; do
        m=${var#*_m}
        echo $m
        # set population size
        if [ $m == 25 ]; then
            generation=2000
        elif [ $m == 50 ] ; then
            generation=3000
        elif [ $m == 75 ] ; then
            generation=4000
        elif [ $m == 100 ] ; then
            generation=5000
        fi
	    . ./compute.ksh run $sel $var $generation
	done ;
    done ;
fi

if [[ $1 != "train" ]]; then
        
    for var in $variators ; do
        . ./compute.ksh bounds dummy $var $generation
    done ;
        
    for sel in $selectors ; do
        for var in $variators ; do
        . ./compute.ksh indicators $sel $var $generation
        done ;
    done ;
        
    for var in $variators ; do
        . ./compute.ksh tests dummy $var $generation
    done ;

    for sel in $selectors ; do
        for var in $variators ; do
        . ./compute.ksh eafs $sel $var $generation
        done ;
    done ;

    for sela in $selectors ; do
        for selb in $selectors ; do
            for var in $variators ; do
            . ./compute.ksh eaftests $sela $selb $var $generation
        done ;
        done ;
    done ;

    for sela in $selectors ; do
        for selb in $selectors ; do
            for var in $variators ; do
            . ./compute.ksh ranktests $sela $selb $var $generation
        done ;
        done ;
    done ;
fi
