#!/usr/bin/perl -w

use strict;
use warnings;
use Tie::IxHash;
use Getopt::Long;
use Statistics::Basic qw(:all);
use perlModule;

###############################################################################
## parse input options
use vars qw($summitFile $bamFile $outDir $sizeFactor $option $winUp $winDown $minClusterHeight $minBlockHeight $distance $scale $blockHeight $noMergeOverlapBlocks $minNFRLength $fileSuffix $help);
$winUp=200;
$winDown=1300;
$minClusterHeight=20;
#$minBlockHeight=5;
$minBlockHeight=20;
$distance=70;
$scale=0.6;
#$blockHeight="rel";
$blockHeight="abs";
$minNFRLength=20;
$fileSuffix="";

GetOptions ("s=s"  => \$summitFile,
            "b=s"  => \$bamFile,
            "o=s"  => \$outDir,
            "z=s"  => \$sizeFactor,
            "p=s"  => \$option,
            "u=i"  => \$winUp,
            "d=i"  => \$winDown,
            "c=s"  => \$minClusterHeight,
            "k=s"  => \$minBlockHeight,
            "x=s"  => \$distance,
            "l=s"  => \$scale,
            "g=s"  => \$blockHeight,
            "m"    => \$noMergeOverlapBlocks,
            "n=s"  => \$minNFRLength,
            "f=s"  => \$fileSuffix,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$summitFile || !$bamFile || !$outDir || !$sizeFactor);

###############################################################################
sub usage {
	print STDERR "\nProgram: findNFR.pl (determine Nucleosome Free Regions (NFR) using ChIP-seq data for histone marks)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: findNFR.pl -s <file> -b <file> -o <dir> -z <float> [OPTIONS]\n";
	print STDERR " -s <file>         [peak summit file(s) from macs2/IDR]\n";
    print STDERR "                   [if multiple, seperate them by a comma]\n";
	print STDERR " -b <file>         [histone ChIP-seq file in BAM format]\n";
	print STDERR " -o <dir>          [directory where output files will be kept]\n";
	print STDERR " -z <float>        [size factor to normalize read expression]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -p <string>       [computation option (default: all):]\n";
    print STDERR "                   a: define blocks and block groups in histone enriched (summit) region\n";
    print STDERR "                   b: define nuclesome free regions\n";
    print STDERR "                   c: determine significant nucleosome free regions\n";
	print STDERR " -u <int>          [nucleotides upstream to summit (default: 200)]\n";
	print STDERR " -d <int>          [nucleotides downstream to summit (default: 1300)]\n";
	print STDERR " -c <int>          [mininum number of read in the block group (default: 20)]\n";
	print STDERR " -k <int>          [mininum number of read in the block (default: 20)]\n";
	print STDERR " -x <int>          [maximum distance between the blocks (default: 70)]\n";
	print STDERR " -l <float>        [scale to define blocks (default: 0.6)]\n";
	print STDERR " -g <string>       [relative block height (abs or rel) (default: abs)]\n";
	print STDERR " -m                [do not merge overlapping blocks]\n";
    print STDERR " -n <int>          [minimum length of nucleosome free region (default: 20)]\n";
    print STDERR " -f <string>       [a string added at the end of output files. useful when running in parallel]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

my $ID=$bamFile;
$ID=~s/^.*\///g;
$ID=~s/\.gz$//g;
my $start=(); my $end=(); my $coor=(); my @data=();

## create output directory, if does not exist
if ( ! -d $outDir) {
    system("mkdir $outDir");
}

## Step-1: define blocks and block groups in histone enriched (summit) region
if(!defined($option) || $option=~/[aA]+/) {
    $summitFile=~s/\,/ /g;
    chomp($summitFile);
    @data=`zless $summitFile | sortBed -i stdin | mergeBed -i stdin`;
    #@data=openFile($summitFile);

    # remove block group file, if already exists
    if(-e "$outDir/$ID.bg$fileSuffix") {
        system("rm $outDir/$ID.bg$fileSuffix");
    }

    foreach my $l(@data) {
        my @F=split(/\s+/,$l);
        $start=$F[1]-$winUp;
        $end=$start+($winUp + $winDown);
        my $coorDown="$F[0]:$start-$end";
        # retrieve reads corresponding to summit region (downstream)
        system("samtools view -b $bamFile $coorDown | bedtools bamtobed -i - | perl -ane '\$F[5]=~s/\-\$/+/g; \$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[4]\\t\$F[5]\\n\";' > $outDir/$ID.tmp$fileSuffix");

        # define block group and blocks corresponding to summit region (downstream)
        #system("blockbuster.x -minClusterHeight $minClusterHeight -minBlockHeight $minBlockHeight -distance $distance -scale $scale -blockHeight $blockHeight $outDir/$ID.tmp$fileSuffix | grep -v \"^>\" | perl -ane 'print \"\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[5]\\t\$F[6]\\t\$F[4]\\n\";' | sortBed -i stdin | perl -ane 'print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t$coorDown\\t\$F[3]\\t+\\n\";' >> $outDir/$ID.bg$fileSuffix");
        system("blockbuster.x -minClusterHeight $minClusterHeight -minBlockHeight $minBlockHeight -distance $distance -scale $scale -blockHeight $blockHeight $outDir/$ID.tmp$fileSuffix | grep -v \"^>\" | perl -ane 'print \"\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[5]\\t\$F[6]\\t\$F[4]\\n\";' | sortBed -i stdin | bedtools merge -c 4 -o collapse -i - | perl -ane 'print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t$coorDown\\t\$F[3]\\t+\\n\";' >> $outDir/$ID.bg$fileSuffix");

        $end=$F[1]+$winUp;
        $start=$end-($winUp + $winDown);
        my $coorUp="$F[0]:$start-$end";
        # retrieve reads corresponding to summit region (upstream)
        system("samtools view -b $bamFile $coorUp | bedtools bamtobed -i - | perl -ane '\$F[5]=~s/\-\$/+/g; \$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[4]\\t\$F[5]\\n\";' > $outDir/$ID.tmp$fileSuffix");

        # define block group and blocks corresponding to summit region (upstream)
        #system("blockbuster.x -minClusterHeight $minClusterHeight -minBlockHeight $minBlockHeight -distance $distance -scale $scale -blockHeight $blockHeight $outDir/$ID.tmp$fileSuffix | grep -v \"^>\" | perl -ane 'print \"\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[5]\\t\$F[6]\\t\$F[4]\\n\";' | sortBed -i stdin | sort -k 2rn,2 -k 3rn,3 | perl -ane 'print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t$coorUp\\t\$F[3]\\t-\\n\";' >> $outDir/$ID.bg$fileSuffix");
        system("blockbuster.x -minClusterHeight $minClusterHeight -minBlockHeight $minBlockHeight -distance $distance -scale $scale -blockHeight $blockHeight $outDir/$ID.tmp$fileSuffix | grep -v \"^>\" | perl -ane 'print \"\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[5]\\t\$F[6]\\t\$F[4]\\n\";' | sortBed -i stdin | bedtools merge -c 4 -o collapse -i - | sort -k 2rn,2 -k 3rn,3 | perl -ane 'print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t$coorUp\\t\$F[3]\\t-\\n\";' >> $outDir/$ID.bg$fileSuffix");

        if($coorUp=~/^$coorDown$/) {
            print STDERR "ERROR: the first and second block do not have same strand despite having same ID\n";
            print STDERR "$l";
            exit(-1);
        }
    }

    system("rm $outDir/$ID.tmp$fileSuffix");
}

## Step-2: define nuclesome free regions
if(!defined($option) || $option=~/[bB]+/) {
    if(-e "$outDir/$ID.bg$fileSuffix") {
        @data=openFile("$outDir/$ID.bg$fileSuffix");
        #@data=openFile("$outDir/test");
    }
    else {
        print STDERR "Cannot find $ID.bg$fileSuffix file. Please run the program with option a first\n";
        usage();
    }

    ## create output file for writing
    open(OUTFILE, ">$outDir/$ID.nfr$fileSuffix") || die $!;

    my %bInfo=(); my @t=(); my %NFR=();
    for(my $i=0; $i<scalar(@data)-2; $i++) {
        my @F1=split(/\s+/, $data[$i]);
        my @F2=split(/\s+/, $data[$i+1]);

        if($F1[3]=~/^$F2[3]$/) {
            $bInfo{'first'}{'chr'}=$F1[0];
            $bInfo{'first'}{'start'}=$F1[1];
            $bInfo{'first'}{'end'}=$F1[2];
            $bInfo{'first'}{'strand'}=$F1[5];
            @t=split(/\,/,$F1[4]);
            $bInfo{'first'}{'expr'}=0;
            foreach(@t) {
                $bInfo{'first'}{'expr'}+=$_;
            }

            $bInfo{'second'}{'chr'}=$F2[0];
            $bInfo{'second'}{'start'}=$F2[1];
            $bInfo{'second'}{'end'}=$F2[2];
            $bInfo{'second'}{'strand'}=$F2[5];
            @t=split(/\,/,$F2[4]);
            $bInfo{'second'}{'expr'}=0;
            foreach(@t) {
                $bInfo{'second'}{'expr'}+=$_;
            }

            if($bInfo{'first'}{'strand'}=~/\+/ && $bInfo{'second'}{'strand'}=~/\+/) {
                $NFR{'chr'}=$bInfo{'first'}{'chr'};
                $NFR{'start'}=$bInfo{'first'}{'end'}+1;
                $NFR{'end'}=$bInfo{'second'}{'start'}-1;
                $NFR{'strand'}="+";
                $NFR{'length'}=($NFR{'end'}-$NFR{'start'})+1;
                $coor="$NFR{'chr'}:$NFR{'start'}-$NFR{'end'}";
                $NFR{'startBlock'}="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
                $NFR{'startBlockExpr'}=$bInfo{'first'}{'expr'};
                $NFR{'endBlock'}="$bInfo{'second'}{'chr'}:$bInfo{'second'}{'start'}-$bInfo{'second'}{'end'}";
                $NFR{'endBlockExpr'}=$bInfo{'second'}{'expr'};
                $NFR{'id'}=$F1[3];

                $NFR{'expr'}=`samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane 'BEGIN{\$expr=sprintf(\"%0.2f\", 1/$sizeFactor);} \$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'`;

                #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
                #print "DEBUG: $bInfo{'first'}{'expr'}\t$bInfo{'adjacent'}{'expr'}\t$NFR{$COUNTER}{'expr'}\t$NFR{$COUNTER}{'length'}\n";
                $NFR{'stddev'}=stddev($bInfo{'first'}{'expr'}, $bInfo{'second'}{'expr'});
                $NFR{'stddev'}=~s/\,//g;
                $NFR{'score'}=sprintf("%0.2f", (($bInfo{'first'}{'expr'}+$bInfo{'second'}{'expr'}))/($NFR{'expr'}));
            }
            elsif($bInfo{'first'}{'strand'}=~/\-/ && $bInfo{'second'}{'strand'}=~/\-/) {
                $NFR{'chr'}=$bInfo{'first'}{'chr'};
                $NFR{'start'}=$bInfo{'second'}{'end'}+1;
                $NFR{'end'}=$bInfo{'first'}{'start'}-1;
                $NFR{'strand'}="-";
                $NFR{'length'}=($NFR{'end'}-$NFR{'start'})+1;
                $coor="$NFR{'chr'}:$NFR{'start'}-$NFR{'end'}";
                $NFR{'startBlock'}="$bInfo{'second'}{'chr'}:$bInfo{'second'}{'start'}-$bInfo{'second'}{'end'}";
                $NFR{'startBlockExpr'}=$bInfo{'second'}{'expr'};
                $NFR{'endBlock'}="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
                $NFR{'endBlockExpr'}=$bInfo{'first'}{'expr'};
                $NFR{'id'}=$F1[3];

                $NFR{'expr'}=`samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane 'BEGIN{\$expr=sprintf(\"%0.2f\", 1/$sizeFactor);} \$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'`;

                #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
                #print "DEBUG: $bInfo{'first'}{'expr'}\t$bInfo{'adjacent'}{'expr'}\t$NFR{$COUNTER}{'expr'}\t$NFR{$COUNTER}{'length'}\n";
                $NFR{'stddev'}=stddev($bInfo{'first'}{'expr'}, $bInfo{'second'}{'expr'});
                $NFR{'stddev'}=~s/\,//g;
                $NFR{'score'}=sprintf("%0.2f", (($bInfo{'first'}{'expr'}+$bInfo{'second'}{'expr'}))/($NFR{'expr'}));
            }
            else {
                print STDERR "WARNING: the first and second block do not have same strand despite having same ID\n";
                print STDERR "--> first block: $data[$i]\n";
                print STDERR "--> second block: $data[$i+1]\n";
            }

            if($NFR{'length'} >= $minNFRLength) {
                print OUTFILE "$NFR{'chr'}\t$NFR{'start'}\t$NFR{'end'}\t$NFR{'id'}\t$NFR{'score'}\t$NFR{'strand'}\t$NFR{'startBlock'}\t$NFR{'startBlockExpr'}\t$NFR{'endBlock'}\t$NFR{'endBlockExpr'}\t$NFR{'stddev'}\t$NFR{'length'}\n";
            }
        }
    }
    close OUTFILE;
}

## Step-3: define non-overlapping and significant nuclesome free regions
if(!defined($option) || $option=~/[cC]+/) {
    if(-e "$outDir/$ID.nfr") {
        @data=`zless $outDir/$ID.nfr$fileSuffix | sort -k 1,1 -k 2n,2 -k 3n,3`;
    }
    else {
        print STDERR "Cannot find $ID.nfr file. Please run the program with option b first\n";
        usage()
    }

    ## create output file for writing non-overlapping NFR
    open(OUTFILE, ">$outDir/$ID.nfr.sig") || die $!;

    ## create output file for writing UCSC tracks of NFR
    open(UCSCFILE, ">$outDir/$ID.nfr.sig.ucsc") || die $!;
    print UCSCFILE "track name=\"Predicted NFR ($ID.nfr.sig)\" description=\"Predicted NFR ($ID.nfr.sig)\" itemRgb=\"On\"\n";
    
    my @F1=(); my @F2=();
    for(my $i=0; $i<(scalar(@data)-1); $i++) {
        if(scalar(@F1)==0) {
            chomp($data[$i]);
            @F1=split(/\s+/, $data[$i]);
        }
        @F2=split(/\s+/, $data[$i+1]);

        ## check overlap between current and next NFR
        my $overlap=checkOverlap($F1[1], $F1[2], $F2[1], $F2[2], 4, $F1[0], $F2[0]);

        ## choose the highest scoring one, if overlapped
        if($overlap) {
            if($F1[4]<$F2[4]) {
                @F1=@F2; 
            }
        }
        else {
            ## print the non-overlapping and highest scoring NFR region
            my $NFR=();
            foreach(@F1) { $NFR.="$_\t"; } $NFR=~s/\t$//g;
            print OUTFILE "$NFR\n";

            my @startBlock=split(/[\:\-]+/,$F1[6]);
            my @endBlock=split(/[\:\-]+/,$F1[8]);
            print UCSCFILE "$startBlock[0]\t$startBlock[1]\t$startBlock[2]\t$F1[6]\t$F1[7]\t.\t$startBlock[1]\t$startBlock[2]\t0,255,0\n";
            print UCSCFILE "$F1[0]\t$F1[1]\t$F1[2]\t$F1[3]\t$F1[4]\t$F1[5]\t$F1[1]\t$F1[2]\t255,0,0\n";
            print UCSCFILE "$endBlock[0]\t$endBlock[1]\t$endBlock[2]\t$F1[8]\t$F1[9]\t.\t$endBlock[1]\t$endBlock[2]\t0,255,0\n";
            @F1=();
        }
    }

    close UCSCFILE;
    close OUTFILE;

    print "$outDir/$ID.nfr.sig\n";
}

## Step-2: define nuclesome free regions (old)
if(defined($option) && $option=~/[dD]+/) {
    if(-e "$outDir/$ID.bg$fileSuffix") {
        @data=openFile("$outDir/$ID.bg$fileSuffix");
        #@data=openFile("$outDir/test");
    }
    else {
        print STDERR "Cannot find $ID.bg$fileSuffix file. Please run the program with option a first\n";
        usage()
    }

    ## create output file for writing
    open(OUTFILE, ">$outDir/$ID.nfr$fileSuffix") || die $!;

    my %bInfo=(); my %NFR=(); my $COUNTER=(); my %seen=();
    foreach my $l(@data) {
        #last if($COUNTER==100);
        my @F=split(/\s+/,$l);
        if(!defined($seen{$F[3]})) {
            ## print predicted nucleosome free regions (NFR)
            if(keys(%NFR)) {
                foreach(reverse sort { $NFR{$a}{'score'} <=> $NFR{$b}{'score'} } keys(%NFR)) {
                        print OUTFILE "$NFR{$_}{'chr'}\t$NFR{$_}{'start'}\t$NFR{$_}{'end'}\t$NFR{$_}{'id'}\t$NFR{$_}{'score'}\t$NFR{$_}{'strand'}\t$NFR{$_}{'startBlock'}\t$NFR{$_}{'startBlockExpr'}\t$NFR{$_}{'endBlock'}\t$NFR{$_}{'endBlockExpr'}\t$NFR{$_}{'stddev'}\t$NFR{$_}{'length'}\n";
                }
            }

            # determine information for summit peak
            %bInfo=(); %seen=(); %NFR=(); $seen{$F[3]}=1; $COUNTER=1;
            $bInfo{'summit'}{'chr'}=$F[0];
            $bInfo{'summit'}{'start'}=$F[1];
            $bInfo{'summit'}{'end'}=$F[2];
            $bInfo{'summit'}{'strand'}=$F[5];
            my @t=split(/\,/,$F[4]);
            $bInfo{'summit'}{'expr'}=0;
            foreach(@t) {
                $bInfo{'summit'}{'expr'}+=$_;
            }
            #$bInfo{'summit'}{'expr'}=sprintf("%0.2f", $bInfo{'summit'}{'expr'}/scalar(@t));
        }
        else {
            # determine information for adjacent peak
            $bInfo{'adjacent'}{'chr'}=$F[0];
            $bInfo{'adjacent'}{'start'}=$F[1];
            $bInfo{'adjacent'}{'end'}=$F[2];
            $bInfo{'adjacent'}{'strand'}=$F[5];
            my @t=split(/\,/,$F[4]);
            $bInfo{'adjacent'}{'expr'}=0;
            foreach(@t) {
                $bInfo{'adjacent'}{'expr'}+=$_;
            }
            #$bInfo{'adjacent'}{'expr'}=sprintf("%0.2f", $bInfo{'adjacent'}{'expr'}/scalar(@t));

            # determine information for putative nucleosome free region
            #printf("%0.2f\n", (($bInfo{'adjacent'}{'start'}-1)-($bInfo{'summit'}{'end'}+1)+1));
            if($bInfo{'adjacent'}{'strand'}=~/^\+$/ && (($bInfo{'adjacent'}{'start'}-1)-($bInfo{'summit'}{'end'}+1)+1) > $minNFRLength) {
                $NFR{$COUNTER}{'chr'}=$bInfo{'summit'}{'chr'};
                $NFR{$COUNTER}{'start'}=$bInfo{'summit'}{'end'}+1;
                $NFR{$COUNTER}{'end'}=$bInfo{'adjacent'}{'start'}-1;
                $NFR{$COUNTER}{'strand'}="+";
                $NFR{$COUNTER}{'length'}=($NFR{$COUNTER}{'end'}-$NFR{$COUNTER}{'start'})+1;
                $coor="$NFR{$COUNTER}{'chr'}:$NFR{$COUNTER}{'start'}-$NFR{$COUNTER}{'end'}";
                $NFR{$COUNTER}{'startBlock'}="$bInfo{'summit'}{'chr'}:$bInfo{'summit'}{'start'}-$bInfo{'summit'}{'end'}";
                $NFR{$COUNTER}{'startBlockExpr'}=$bInfo{'summit'}{'expr'};
                $NFR{$COUNTER}{'endBlock'}="$bInfo{'adjacent'}{'chr'}:$bInfo{'adjacent'}{'start'}-$bInfo{'adjacent'}{'end'}";
                $NFR{$COUNTER}{'endBlockExpr'}=$bInfo{'adjacent'}{'expr'};
                $NFR{$COUNTER}{'id'}=$F[3];

                $NFR{$COUNTER}{'expr'}=`samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane 'BEGIN{\$expr=sprintf(\"%0.2f\", 1/$sizeFactor);} \$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'`;

                #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
                #print "DEBUG: $bInfo{'summit'}{'expr'}\t$bInfo{'adjacent'}{'expr'}\t$NFR{$COUNTER}{'expr'}\t$NFR{$COUNTER}{'length'}\n";
                $NFR{$COUNTER}{'stddev'}=stddev($bInfo{'summit'}{'expr'}, $bInfo{'adjacent'}{'expr'});
                $NFR{'stddev'}=~s/\,//g;
                #$NFR{$COUNTER}{'stddev'}=2;
                #$NFR{$COUNTER}{'score'}=sprintf("%0.2f", (($bInfo{'summit'}{'expr'}+$bInfo{'adjacent'}{'expr'}))/($NFR{$COUNTER}{'expr'}/$NFR{$COUNTER}{'length'}));
                $NFR{$COUNTER}{'score'}=sprintf("%0.2f", (($bInfo{'summit'}{'expr'}+$bInfo{'adjacent'}{'expr'}))/($NFR{$COUNTER}{'expr'}));

                $COUNTER++;
            }
            elsif($bInfo{'adjacent'}{'strand'}=~/^\-$/ && (($bInfo{'summit'}{'start'}-1)-($bInfo{'adjacent'}{'end'}+1)+1) > $minNFRLength) {
                $NFR{$COUNTER}{'chr'}=$bInfo{'summit'}{'chr'};
                $NFR{$COUNTER}{'start'}=$bInfo{'adjacent'}{'end'}+1;
                $NFR{$COUNTER}{'end'}=$bInfo{'summit'}{'start'}-1;
                $NFR{$COUNTER}{'strand'}="-";
                $NFR{$COUNTER}{'length'}=($NFR{$COUNTER}{'end'}-$NFR{$COUNTER}{'start'})+1;
                $coor="$NFR{$COUNTER}{'chr'}:$NFR{$COUNTER}{'start'}-$NFR{$COUNTER}{'end'}";
                $NFR{$COUNTER}{'startBlock'}="$bInfo{'adjacent'}{'chr'}:$bInfo{'adjacent'}{'start'}-$bInfo{'adjacent'}{'end'}";
                $NFR{$COUNTER}{'startBlockExpr'}=$bInfo{'adjacent'}{'expr'};
                $NFR{$COUNTER}{'endBlock'}="$bInfo{'summit'}{'chr'}:$bInfo{'summit'}{'start'}-$bInfo{'summit'}{'end'}";
                $NFR{$COUNTER}{'endBlockExpr'}=$bInfo{'summit'}{'expr'};
                $NFR{$COUNTER}{'id'}=$F[3];

                $NFR{$COUNTER}{'expr'}=`samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane 'BEGIN{\$expr=sprintf(\"%0.2f\", 1/$sizeFactor);} \$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'`;

                #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
                #print "DEBUG: $bInfo{'summit'}{'expr'}\t$bInfo{'adjacent'}{'expr'}\t$NFR{$COUNTER}{'expr'}\t$NFR{$COUNTER}{'length'}\n";
                $NFR{$COUNTER}{'stddev'}=stddev($bInfo{'summit'}{'expr'}, $bInfo{'adjacent'}{'expr'});
                $NFR{'stddev'}=~s/\,//g;
                #$NFR{$COUNTER}{'stddev'}=2;
                #$NFR{$COUNTER}{'score'}=sprintf("%0.2f", (($bInfo{'summit'}{'expr'}+$bInfo{'adjacent'}{'expr'}))/($NFR{$COUNTER}{'expr'}/$NFR{$COUNTER}{'length'}));
                $NFR{$COUNTER}{'score'}=sprintf("%0.2f", (($bInfo{'summit'}{'expr'}+$bInfo{'adjacent'}{'expr'}))/($NFR{$COUNTER}{'expr'}));

                $COUNTER++;
            }

            # reassign current peak as summit peak, if its expression is higher than original summit peak
            if($bInfo{'adjacent'}{'expr'} > $bInfo{'summit'}{'expr'}) {
                %{$bInfo{'summit'}}=%{$bInfo{'adjacent'}};
            }
        }
    }

    ## print predicted nucleosome free regions (NFR)
    if(keys(%NFR)) {
        foreach(reverse sort { $NFR{$a}{'score'} <=> $NFR{$b}{'score'} } keys(%NFR)) {
                print OUTFILE "$NFR{$_}{'chr'}\t$NFR{$_}{'start'}\t$NFR{$_}{'end'}\t$NFR{$_}{'id'}\t$NFR{$_}{'score'}\t$NFR{$_}{'strand'}\t$NFR{$_}{'startBlock'}\t$NFR{$_}{'startBlockExpr'}\t$NFR{$_}{'endBlock'}\t$NFR{$_}{'endBlockExpr'}\t$NFR{$_}{'stddev'}\t$NFR{$_}{'length'}\n";
        }
    }
    close OUTFILE;
}

exit(0);
