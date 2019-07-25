#!/bin/bash
#PBS -l nodes=1:ppn=4

CONFIG_FILTER="multiIn"

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
    echo "             **OR**"
    echo "             [input configuration file containing bed file information]"
    echo "             [<id> <uniqueId> <bed file> (id should be multiIn]"
    echo "[OPTIONS]"
    echo " -j <string> [unique name to describe each input bed file separated by a comma]"
    echo "             **OR**"
    echo "             [can be specified in the config file]"
    echo " -f <string> [filter in input BED files for this input string parameter]"
    echo " -F <string> [filter out input BED files for this input string parameter]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:j:f:F:h ARG; do
	case "$ARG" in
		i) BEDFILE=$OPTARG;;
        j) NAME=$OPTARG;;
        f) FILTER_IN=$OPTARG;;
        F) FILTER_OUT=$OPTARG;;
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

## check if input is BED files or configuration file containing BED file information
INPUT=$(echo $BEDFILE | perl -ane '$_=~s/\,.*//g; print $_;')
if [ "$(sortBed -i $INPUT 2>/dev/null | wc -l)" -le 0 ]; then
    ## read configuration file
    NAME=$(cat $BEDFILE | perl -ane '
        if($_=~/^'$CONFIG_FILTER'/) {
            if(!$seen{$F[1]}) { 
                $file_name.="$F[1],";
                $seen{$F[1]}=1;
            }
        } END {
            $file_name=~s/\,$//g;
            print "$file_name\n";
        }'
    )

    INPUT=$(cat $BEDFILE | perl -ane '
        if($_=~/^'$CONFIG_FILTER'/) {
            $file.="$F[2],";
        } END {
            $file=~s/\,$//g;
            print "$file\n";
        }'
    )
    BEDFILE=$INPUT
fi
#echo -e "$BEDFILE\t$NAME"; exit

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
    if [ ! -z "$FILTER_IN" ]; then
        zless ${BEDFILES[$i]} | perl -ane 'for($i=0; $i<scalar(@F); $i++) { if($F[$i]=~/^'$FILTER_IN'$/) { print $_; last; } }' |  bedtools sort -i - > ${TMP_NAME[$i]}.bed
    elif [ ! -z "$FILTER_OUT" ]; then
        zless ${BEDFILES[$i]} | perl -ane '$found=0; for($i=0; $i<scalar(@F); $i++) { if($F[$i]=~/^'$FILTER_OUT'$/i) { $found=1; last; } } if($found==0) { print $_; }' |  bedtools sort -i - > ${TMP_NAME[$i]}.bed
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
    bedtools multiinter -i $COMMAND_BED -names $COMMAND_NAME | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; } } elsif($last_coor!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
    ## last_counter=$F[3] removed from else condition (Feb 22, 2017)
    #bedtools multiinter -i $COMMAND_BED -names $COMMAND_NAME | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_coor!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
    #bedtools multiinter -i $COMMAND_BED -names $COMMAND_NAME | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_counter!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
else
    bedtools multiinter -i $COMMAND_BED | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; } } elsif($last_coor!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
    ## last_counter=$F[3] removed from else condition (Feb 22, 2017)
    #bedtools multiinter -i $COMMAND_BED | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_coor!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
    #bedtools multiinter -i $COMMAND_BED | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_counter!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
fi | while read line; do 
    if [ ! -z "$FILTER_IN" ]; then 
        echo $line$'\t'$FILTER_IN; 
    else
        echo "$line";
    fi
done

## remove temporary files
for (( i=0; i<$BEDFILES_COUNT; i++ )); do
    rm ${TMP_NAME[$i]}.bed
done
