#!/bin/bash
#PBS -l nodes=1:ppn=4

GENOME="mm9"
PROCESSOR=1

#### usage ####
usage() {
	echo Program: "motifStatAna (compute motif statistics across multiple samples)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: motifStatAna -i <dir> -k <file> -o <dir> [OPTIONS]"
	echo "Options:"
    echo " -i <dir>    [directory containing result from prior run of motifAna script]"
    echo "             [if multiple, seperate them by a comma]"
    echo " -k <file>   [motif file for which enrichment will be analyzed]"
    echo "             [if multiple, separate them by a comma]"
    echo " -o <dir>    [output directory to keep pdf image file]"
    echo "[OPTIONS]"
    echo " -g <string> [genome for which to perform the analysis (default: mm9)]"
    echo " -p <int>    [number of processors to use (default: 1)]"
    echo " -j <string> [name of motifs that are of interest]"
    echo "             [if multiple separate them by a comma]"
    echo "             [If not provided analysis will be done using all motifs]" 
    echo " -h          [help]"
    echo
    exit 0
}

#### parse options ####
while getopts i:k:o:g:p:j:h ARG; do
	case "$ARG" in
		i) INDIR=$OPTARG;;
        o) OUTDIR=$OPTARG;;
        k) MOTIF_FILE=$OPTARG;;
        g) GENOME=$OPTARG;;
        p) PROCESSOR=$OPTARG;;
        j) MOTIF_NAME=$OPTARG;;
        h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$INDIR" -o -z "$MOTIF_FILE" -o -z "$OUTDIR" -o "$HELP" ]; then
	usage
fi

###################
#helperfunction
function wait_for_jobs_to_finish {
    for job in `jobs -p`
    do
        echo $job
        wait $job
    done
    echo $1
}
###############

<<"COMMENT1"
COMMENT1

echo -n "Populating files based on input genome, $GENOME (`date`).. "
if [ "$GENOME" == "mm9" ]; then
    GENOME_FILE="/home/pundhir/project/genome_annotations/mouse.mm9.genome"
    REPEAT_FILE="/home/pundhir/project/genome_annotations/mouse.mm9.simpleRepeat.gz"
    GENOME_MOTIF="mm9r"
elif [ "$GENOME" == "hg19" ]; then
    GENOME_FILE="/home/pundhir/project/genome_annotations/human.hg19.genome"
    REPEAT_FILE="/home/pundhir/project/genome_annotations/human.hg19.simpleRepeat.gz"
    GENOME_MOTIF="hg19r"
elif [ "$GENOME" == "danRer7" ]; then
    GENOME_FILE="/home/pundhir/project/genome_annotations/zebrafish.danRer7.genome"
    REPEAT_FILE="/home/pundhir/project/genome_annotations/zebrafish.danRer7.simpleRepeat.gz"
    GENOME_MOTIF="danRer7r"
else
    echo "Presently the program only support analysis for mm9, hg19 or danRer7"
    echo
    usage
fi
echo done

## parse input directories
oIFS=$IFS
IFS=","
INDIRS=($INDIR)
IFS=$oIFS

<<"COMMENT"
COMMENT

STATDIR=""
for (( i=0; i<${#INDIRS[@]}; i++ )); do
    echo -n "Create output directory.. "
    if [ ! -d "${INDIRS[$i]}/stat" ]; then
        mkdir -p ${INDIRS[$i]}/stat
    fi
    echo "done"

    echo -n "Compute motif statistics.. "
    if [ ! -z "$MOTIF_NAME" ]; then
        motifAna -i ${INDIRS[$i]}/REGIONS_INTEREST.bed -o ${INDIRS[$i]}/stat -m 2 -l $MOTIF_FILE -j $MOTIF_NAME
    else
        motifAna -i ${INDIRS[$i]}/REGIONS_INTEREST.bed -o ${INDIRS[$i]}/stat -m 2 -l $MOTIF_FILE
    fi
    STATDIR="$STATDIR,${INDIRS[$i]}/stat"
    echo "done"
done
STATDIR=$(echo $STATDIR | perl -ane '$_=~s/^\,//g; print $_;')

echo -n "Plot motif statistics.. "
motifStatAna.R -i $STATDIR -o $OUTDIR
echo "done"
