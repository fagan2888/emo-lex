#! /bin/bash
#
# compute.sh: UNIX/LINUX complement to compute.bat
#
# Eckart Zitzler, ETH Zurich, Feb 18, 2005

#=================================================
# DO NOT MAKE ANY CHANGES BELOW THIS LINE
#=================================================
os=linux
# 1: task, 2: selector, 3: variator, 4: generation, 5: trial #
case $1 in
    run)
    
	cd "$2_$os" ;
	cp "../$3_$os/$2_$3_cfg" . ;
	./"$2" "$2"_param.txt "$2"_"$3"_ 0.1 &
	cd ../"$3"_"$os" ;
	./"$3" "$3"_param.txt "$2"_"$3"_ 0.1 &
	cd ../monitor_"$os"
    
    echo "./monitor ../"$3"_"$os"/"$3"_param.txt ../"$3"_"$os"/"$2"_"$3"_
    ../"$2"_"$os"/"$2"_param.txt ../"$2"_"$os"/"$2"_"$3"_ monitor_param_"$4".txt ../runs/"$3"_"$2" 0.1 ;"
    
    ./monitor ../"$3"_"$os"/"$3"_param.txt ../"$3"_"$os"/"$2"_"$3"_ ../"$2"_"$os"/"$2"_param.txt ../"$2"_"$os"/"$2"_"$3"_ monitor_param_"$4".txt ../runs/"$3"_"$2" 0.1 ;
	cd ..
	;;
    bounds) 
	cd runs ;
	cat "$3"_*."$4" > ../tests/"$3"."$4" ;
	cd ../tools_"$os" ;
	./bound ../tests/"$3"."$4" ../tests/"$3"_bound."$4" ;
	./normalize ../tests/"$3"_bound."$4" ../tests/"$3"."$4" ../tests/"$3"_norm."$4" ;
	./filter  ../tests/"$3"_norm."$4" ../tests/"$3"_ref."$4" ;
	cd ..
	;;
    indicators)
	cd tools_"$os" ;
	./normalize ../tests/"$3"_bound."$4" ../runs/"$3"_"$2"."$4" ../tests/"$3"_"$2"_norm."$4" ;
	cd ../indicators_"$os" ;
	./hyp_ind ../tests/"$3"_"$2"_norm."$4" ../tests/"$3"_ref."$4" ../tests/"$3"_"$2"_hyp."$4" ;
	echo >> ../tests/"$3"_"$2"_hyp."$4" ;
	./eps_ind ../tests/"$3"_"$2"_norm."$4" ../tests/"$3"_ref."$4" ../tests/"$3"_"$2"_eps."$4" ;
	echo >> ../tests/"$3"_"$2"_eps."$4" ;
	./r_ind ../tests/"$3"_"$2"_norm."$4" ../tests/"$3"_ref."$4" ../tests/"$3"_"$2"_r."$4" ;
	echo >> ../tests/"$3"_"$2"_r."$4" ;
	cd ..
	;;
    tests)
	cd tests ;
	cat "$3"_*_eps."$4" > "$3"_eps."$4" ;
	cat "$3"_*_hyp."$4" > "$3"_hyp."$4" ;
	cat "$3"_*_r."$4" > "$3"_r."$4" ;
	cd ../statistics_"$os";
	./kruskal-wallis ../tests/"$3"_eps."$4" kruskalparam.txt ../tests/"$3"_eps_kruskal."$4" ;
	./kruskal-wallis ../tests/"$3"_hyp."$4" kruskalparam.txt ../tests/"$3"_hyp_kruskal."$4" ;
	./kruskal-wallis ../tests/"$3"_r."$4" kruskalparam.txt ../tests/"$3"_r_kruskal."$4" ;
	./mann-whit ../tests/"$3"_eps."$4" emptyparam.txt ../tests/"$3"_eps_mann."$4" ;
	./mann-whit ../tests/"$3"_hyp."$4" emptyparam.txt ../tests/"$3"_hyp_mann."$4" ;
	./mann-whit ../tests/"$3"_r."$4" emptyparam.txt ../tests/"$3"_r_mann."$4" ;
        ./wilcoxon-sign ../tests/"$3"_eps."$4" emptyparam.txt ../tests/"$3"_eps_wilcoxon."$4" ;
        ./wilcoxon-sign ../tests/"$3"_hyp."$4" emptyparam.txt ../tests/"$3"_hyp_wilcoxon."$4" ;
        ./wilcoxon-sign ../tests/"$3"_r."$4" emptyparam.txt ../tests/"$3"_r_wilcoxon."$4" ;
        ./fisher-matched ../tests/"$3"_eps."$4" fisherparam.txt ../tests/"$3"_eps_fisherm."$4" ;
        ./fisher-matched ../tests/"$3"_hyp."$4" fisherparam.txt ../tests/"$3"_hyp_fisherm."$4" ;
        ./fisher-matched ../tests/"$3"_r."$4" fisherparam.txt ../tests/"$3"_r_fisherm."$4" 
        ./fisher-indep ../tests/"$3"_eps."$4" fisherparam.txt ../tests/"$3"_eps_fisheri."$4" ;
        ./fisher-indep ../tests/"$3"_hyp."$4" fisherparam.txt ../tests/"$3"_hyp_fisheri."$4" ;
        ./fisher-indep ../tests/"$3"_r."$4" fisherparam.txt ../tests/"$3"_r_fisheri."$4"
	cd ..
	;;
    eafs)
        cd attainment_"$os" ;
	./eaf -o ../tests/"$3"_"$2"_eaf."$4" ../tests/"$3"_"$2"_norm."$4" ;
	cd ..
	;;
    eaftests)
    if [[ $2 < "$3" ]] ; then 
        cd attainment_"$os" ;
	eaf -i ../tests/tmp ../tests/"$4"_"$2"_norm ../tests/"$4"_"$3"_norm > ../tests/dump ;
	eaf-test ../tests/tmp > ../tests/"$4"_"$2"_"$3"_eaftest ;
	rm ../tests/dump ;
	rm ../tests/tmp ;
	cd .. ;
    fi
	;;
    ranktests)
    if [[ $2 < $3 ]] ; then 
    	cd indicators_"$os" ;
	dominance-rank ../tests/"$4"_"$2"_norm ../tests/"$4"_"$3"_norm ../tests/dump > ../tests/tmp ;
	cd ../statistics_"$os" ;
	mann-whit ../tests/tmp emptyparam.txt ../tests/"$4"_"$2"_"$3"_ranktest ;
	rm ../tests/dump ;
	rm ../tests/tmp ;
	cd .. ;
    fi
    ;;
esac
