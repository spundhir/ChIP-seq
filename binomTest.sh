#!/bin/bash

#### usage ####
usage() {
	echo Program: "binorm.sh (perform binomial test)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: binorm.sh -i <file> -j <file> -k <string>"
	echo "Options:"
    echo " -i <file>   [input file containing reference features, eg. genome segementation in BED format]"
    echo " -j <file>   [input file containing features of interest, eg. NFRs in BED format (can be stdin)]"
	echo "[optional]"
    echo " -k <string> [list of specific reference features to search for eg. E, TSS etc]"
    echo "             [if multiple, seperate them by a comma]"
    echo " -l <string> [list of specific features to filter regions of interest]"
	echo " -h          [help]"
	echo
	exit 0
}

MAPPING_FREQUENCY=1

#### parse options ####
while getopts i:j:k:l:h ARG; do
	case "$ARG" in
        i) FEATURES_REF=$OPTARG;;
        j) FEATURES_INT=$OPTARG;;
        k) FILTER_REF=$OPTARG;;
        l) FILTER_INT=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$FEATURES_REF" -o -z "$FEATURES_INT" -o "$HELP" ]; then
	usage
fi

## create temporary BED file if input is from stdin
if [ "$FEATURES_INT" == "stdin" ]; then
    TMP=$(date | md5sum | cut -f 1 -d " ")
    #TMP=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)
    while read LINE; do
    echo ${LINE}
    done | perl -ane '$line=""; foreach(@F) { $line.="$_\t"; } $line=~s/\t$//g; print "$line\n";' > $TMP
    FEATURES_INT=$TMP
else
    TMP=$(date | md5sum | cut -f 1 -d " ")
    #TMP=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)
    scp $FEATURES_INT $TMP
    FEATURES_INT=$TMP
fi

if [ ! -z "$FILTER_INT" ]; then
    perl -ane 'if($_=~/\s+'$FILTER_INT'\s+/) { print $_; }' $FEATURES_INT > $FEATURES_INT.tmp;
    mv $FEATURES_INT.tmp $FEATURES_INT
fi

#echo "$FEATURES_REF\t$FEATURES_INT\t$FILTER_REF";

if [ -z "$FILTER_REF" ]; then
    N=`zless $FEATURES_INT | wc -l | cut -f 1 -d " "`;
    Pr=`zless $FEATURES_REF | perl -ane '$cov+=($F[2]-$F[1])+1; END { printf("%0.3f", $cov/3000000000); }'`;
    mean=`echo $Pr | perl -ane '$mean='$N'*$_; printf("%0.3f", $mean);'`;
    stdev=`echo $Pr | perl -ane '$stdev=sqrt('$N'*$_*(1-$_)); printf("%0.3f", $stdev);'`;

    overlap=`intersectBed -a $FEATURES_INT -b $FEATURES_REF -u | wc -l`;
    pvalue=`Rscript /home/pundhir/software/myScripts/PredictNFR_v0.01/pnorm.R $overlap $mean $stdev | cut -f 2 -d " "`;
    per=`perl -e '$per=('$overlap'*100)/'$N'; printf("%0.2f", $per);'`;
    #echo -e "$entity\t$N\t$Pr\t$mean\t$stdev\t$overlap\t$per"; exit;
    file=`echo $FEATURES_REF | sed 's/^.*\///g'`;
    echo -e "$file\tNA\t$N\t$overlap\t$mean\t$stdev\t$pvalue\t$per";
else
    IFS=","
    FEATURES=($FILTER_REF)
    IFS=""
    FIELD=$(intersectBed -a $FEATURES_INT -b $FEATURES_REF -wo | head -n 1 | perl -ane '$field=scalar(@F); printf("%drn,%d", $field, $field);')
    NCOL=$(intersectBed -a $FEATURES_INT -b $FEATURES_REF -wo | head -n 1 | perl -ane 'print scalar(@F);')
    #echo -e "$FIELD\t$NCOL"

    for (( i=0; i<${#FEATURES[@]}; i++ )); do
        entity=${FEATURES[$i]};
        ENTITY_COL=$(intersectBed -a $FEATURES_INT -b $FEATURES_REF -wo | grep -w $entity | head -n 1 | perl -ane 'BEGIN { $col=1; } foreach(@F) { if($_=~/^'$entity'$/) { print $col; last; } $col++; }')
        N=`zless $FEATURES_INT | wc -l | cut -f 1 -d " "`;
        Pr=`zgrep -w $entity $FEATURES_REF | perl -ane '$cov+=($F[2]-$F[1])+1; END { printf("%0.3f", $cov/3000000000); }'`;
        mean=`echo $Pr | perl -ane '$mean='$N'*$_; printf("%0.3f", $mean);'`;
        stdev=`echo $Pr | perl -ane '$stdev=sqrt('$N'*$_*(1-$_)); printf("%0.3f", $stdev);'`;

        if [ ! -z "$ENTITY_COL" ]; then
            #echo "intersectBed -a $FEATURES_INT -b $FEATURES_REF -wao | sort -k 1,1 -k 2n,2 -k 3n,3 -k $FIELD | perl -ane '\$key=\"\$F[0]_\$F[1]_\$F[2]\"; if(!defined(\$seen{\$key})) { print \$_; \$seen{\$key}=1;}' | cut -f $ENTITY_COL | sort | uniq -c | sed -E 's/^\s+//g' | grep -w $entity | cut -f 1 -d \" \""; exit;
            overlap=`~/software/bedtools2-2.19.1/bin/intersectBed -a $FEATURES_INT -b $FEATURES_REF -wao | sort -k 1,1 -k 2n,2 -k 3n,3 -k $FIELD | perl -ane '$key="$F[0]_$F[1]_$F[2]"; if(!defined($seen{$key})) { print $_; $seen{$key}=1;}' | cut -f $ENTITY_COL | sort | uniq -c | sed -E 's/^\s+//g' | perl -ane 'if($F[1]=~/^'$entity'$/) { print $_; }' | cut -f 1 -d " "`;
            if [ -z "$overlap" ]; then
                overlap=0
            fi
       else 
            overlap=0
        fi
        #echo -e "Entity: $entity; N: $N; Pr: $Pr; Mean: $mean; Stdev: $stdev; Overlap: $overlap"; exit;

        pvalue=`Rscript /home/pundhir/software/myScripts/PredictNFR_v0.01/pnorm.R $overlap $mean $stdev | cut -f 2 -d " "`;
        per=`perl -e '$per=('$overlap'*100)/'$N'; printf("%0.2f", $per);'`;
        #echo -e "$entity\t$N\t$Pr\t$mean\t$stdev\t$overlap\t$pvalue";
        file=`echo $FEATURES_REF | sed 's/^.*\///g'`;
        echo -e "$file\t${FEATURES[$i]}\t$N\t$overlap\t$mean\t$stdev\t$pvalue\t$per";
    done

    ## performing enrichment analysis for features that do not overlap with reference
    N=`zless $FEATURES_INT | wc -l | cut -f 1 -d " "`;
    Pr=`zless $FEATURES_REF | perl -ane '$cov+=($F[2]-$F[1])+1; END { printf("%0.3f", (3000000000-$cov)/3000000000); }'`;
    mean=`echo $Pr | perl -ane '$mean='$N'*$_; printf("%0.3f", $mean);'`;
    stdev=`echo $Pr | perl -ane '$stdev=sqrt('$N'*$_*(1-$_)); printf("%0.3f", $stdev);'`;
    overlap=`intersectBed -a $FEATURES_INT -b $FEATURES_REF -v | wc -l`;
    pvalue=`Rscript /home/pundhir/software/myScripts/PredictNFR_v0.01/pnorm.R $overlap $mean $stdev | cut -f 2 -d " "`;
    per=`perl -e '$per=('$overlap'*100)/'$N'; printf("%0.2f", $per);'`;
    #echo -e "$entity\t$N\t$Pr\t$mean\t$stdev\t$overlap\t$pvalue";
    file=`echo $FEATURES_REF | sed 's/^.*\///g'`;
    echo -e "$file\tother\t$N\t$overlap\t$mean\t$stdev\t$pvalue\t$per";
fi

if [ ! -z "$TMP" ]; then
    rm $TMP
fi
