#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d         : usual run"
    echo " [0] --d_all    : usual run all"
    echo " [0] -g         : generate in file list"
    echo " [0] --g_all    : generate in file list all"
    echo " [0] --copy_all : copy all"
    echo " [0] -c         : recompile"
    echo " [0] -h         : print help"
}

fileList="file.list"
fileListAll="fileAll.list"
histOut="hist.root"
histOutAll="histAll.root"
merg_pdf_sh="merg_pdf.sh"

#source D_cord.sh -c; source D_cord.sh -g; source D_cord.sh -d; source merg_pdf.sh

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./D_cord 0 $fileList $histOut
    elif [ "$1" = "--d_all" ]; then
	./D_cord 0 $fileListAll $histOutAll
    elif [ "$1" = "-g" ]; then
	rm -rf $fileList
	rm -rf $merg_pdf_sh
	for i in $(seq 1 100); do
	    #fname="./test_outdir_0/hist_phi_csvName"$i".csv"
	    fname="./outdir_0/hist_phi_csvName"$i".csv"
	    if [ -f $fname ]; then
		echo $fname | tee -a $fileList
	    fi
	done
	out_pdf_list=$(cat $fileList | awk '{print $1".pdf"}' | xargs)
	echo "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=c1_merge.pdf $out_pdf_list" > $merg_pdf_sh
    elif [ "$1" = "--g_all" ]; then
	rm -rf $fileListAll
	#
	for j in $(seq 0 9); do
	    for i in $(seq 1 10000); do
		fname="./outdir_"$j"/hist_phi_csvName"$i".csv"
		if [ -f $fname ]; then
		    echo $fname | tee -a $fileListAll
		fi
	    done
	done
    elif [ "$1" = "--copy_all" ]; then
	for j in $(seq 0 9); do
	    fname="../../work/CTA/cta-lstchain/cta-lstchain_mu_images_ana/outdir_"$j
	    cmd="cp -r $fname ."
	    $cmd
	done
    elif [ "$1" = "-c" ]; then
	make clean; make;
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"


