#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d  : usual run"
    echo " [0] -g  : generate in file list"
    echo " [0] -c  : recompile"
    echo " [0] -h  : print help"
}

fileList="file.list"
histOut="hist.root"
merg_pdf_sh="merg_pdf.sh"

#source D_cord.sh -c; source D_cord.sh -g; source D_cord.sh -d; source merg_pdf.sh

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./D_cord 0 $fileList $histOut
    elif [ "$1" = "-g" ]; then
	rm -rf $fileList
	rm -rf $merg_pdf_sh
	for i in $(seq 1 100); do
	    fname="./test_outdir_0/hist_phi_csvName"$i".csv"
	    if [ -f $fname ]; then
		echo $fname | tee -a $fileList
	    fi
	done
	out_pdf_list=$(cat $fileList | awk '{print $1".pdf"}' | xargs)
	echo "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=c1_merge.pdf $out_pdf_list" > $merg_pdf_sh
    elif [ "$1" = "-c" ]; then
	make clean; make;
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"


