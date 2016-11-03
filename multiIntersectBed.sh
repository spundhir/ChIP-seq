#!/bin/bash
#PBS -l nodes=1:ppn=4

#### usage ####
usage() {
    echo
	echo Program: "multiIntersectBed.sh (intersect genomic coordinates from multiple BED files)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: multiIntersectBed.sh -i <files> [OPTIONS]"
    echo "Note: Different from multiIntersectBed (bedtools) since, this only selects the single most observed coordinate among the consecutive overlapping coordinates"
	echo "Options:"
    echo " -i <files>  [input BED files seperated by a comma]"
    echo "[OPTIONS]"
    echo " -j <string> [names to describe each input files seperated by a comma]"
    echo " -f <string> [filter input BED files for this input string parameter]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:j:f:h ARG; do
	case "$ARG" in
		i) BEDFILE=$OPTARG;;
        j) NAME=$OPTARG;;
        f) FILTER=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$BEDFILE" -o "$HELP" ]; then
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

## parse input bam files in an array
IFS=","
BEDFILES=($BEDFILE)
BEDFILES_COUNT=${#BEDFILES[@]}
IFS=" "

## initialize name parameter, if provided
if [ ! -z "$NAME" ]; then
    IFS=","
    NAMES=($NAME)
    NAMES_COUNT=${#NAMES[@]}
    IFS=" "
else
    NAMES_COUNT=0
fi

if [ "$BEDFILES_COUNT" -lt 2 -o ! -z "$NAME" -a "$BEDFILES_COUNT" -ne "$NAMES_COUNT" ]; then
    echo -n "minimum two input bed files are required as input. Also provide name for each input bed file";
    usage
fi

COMMAND_BED=""
COMMAND_NAME=""
for (( i=0; i<$BEDFILES_COUNT; i++ )); do
    TMP_NAME[i]=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)
    #TMP_NAME[i]=$RANDOM
    if [ ! -z "$FILTER" ]; then
        zless ${BEDFILES[$i]} | perl -ane 'for($i=0; $i<scalar(@F); $i++) { if($F[$i]=~/^'$FILTER'$/) { print $_; last; } }' |  bedtools sort -i - > ${TMP_NAME[$i]}.bed
    else
        bedtools sort -i ${BEDFILES[$i]} > ${TMP_NAME[$i]}.bed
    fi
    COMMAND_BED="$COMMAND_BED ${TMP_NAME[$i]}.bed"
    if [ ! -z "$NAME" ]; then
        COMMAND_NAME="$COMMAND_NAME ${NAMES[$i]}"
    fi
done
wait

#echo "bedtools multiinter -i $COMMAND_BED -names $COMMAND_NAME"; exit

if [ ! -z "$NAME" ]; then
    bedtools multiinter -i $COMMAND_BED -names $COMMAND_NAME | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_coor!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
    #bedtools multiinter -i $COMMAND_BED -names $COMMAND_NAME | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_counter!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
else
    bedtools multiinter -i $COMMAND_BED | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_coor!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
    #bedtools multiinter -i $COMMAND_BED | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_counter!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
fi | while read line; do 
    if [ ! -z "$FILTER" ]; then 
        echo $line$'\t'$FILTER; 
    else
        echo "$line";
    fi
done

## remove temporary files
for (( i=0; i<$BEDFILES_COUNT; i++ )); do
    rm ${TMP_NAME[$i]}.bed
done
