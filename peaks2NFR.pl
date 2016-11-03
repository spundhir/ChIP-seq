#!/usr/bin/perl -w

use strict;
use warnings;
use Tie::IxHash;
use Getopt::Long;
use Statistics::Basic qw(:all);
use perlModule;

###############################################################################
## parse input options
use vars qw($peakFile $bamFile $sizeFactor $pvalue $minNFRLength $maxNFRLength $extend $genome $fileSuffix $help);
$pvalue=0.05;
$minNFRLength=20;
$maxNFRLength=1000;
$extend=0;
$genome="mm9";
$fileSuffix="";

GetOptions ("i=s"  => \$peakFile,
            "b=s"  => \$bamFile,
            "z=s"  => \$sizeFactor,
            "p=s"  => \$pvalue,
            "n=s"  => \$minNFRLength,
            "v=s"  => \$maxNFRLength,
            "e=s"  => \$extend,
            "y=s"  => \$genome,
            "f=s"  => \$fileSuffix,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$peakFile || !$bamFile || !$sizeFactor);

###############################################################################
sub usage {
	print STDERR "\nProgram: peaks2NFR.pl (determine Nucleosome Free Regions (NFR) using peaks called from bed2peaks script)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: peaks2NFR.pl -i <file> -b <file> -o <dir> -z <float> [OPTIONS]\n";
	print STDERR " -i <file>         [input file containing peaks region in BED format (bed2peaks output format)]\n";
    print STDERR "                   [can be stdin]\n";
	print STDERR " -b <file>         [histone ChIP-seq file in BAM format]\n";
	print STDERR " -o <file>         [output file containing predicted NFRs]\n";
	print STDERR " -z <float>        [size factor to normalize read expression]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -p <float>        [pvalue for read enrichment (default: 0.05)]\n";
    print STDERR " -n <int>          [minimum length of nucleosome free region (default: 20)]\n";
    print STDERR " -v <int>          [maximum length of nucleosome free region (default: 1000)]\n";
    print STDERR " -e <int>          [extend 3' end of reads by input number of bases (default: 0)]\n";
    print STDERR " -y <string>       [genome (default: mm9)]\n";
    print STDERR " -f <string>       [a string added at the end of output files. useful when running in parallel]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

my $ID=$bamFile;
$ID=~s/^.*\///g;
$ID=~s/\.gz$//g;
my $start=(); my $end=(); my $coor=(); my @data=();

## Step-1: define nuclesome free regions
if(-e "$peakFile") {
    @data=openFile($peakFile);
}
elsif($peakFile=~/^stdin$/) {
    @data=openFile();
    $peakFile="stdin"
}
else {
    #print STDERR "Cannot find $peakFile\n";
    #usage();
}

## create output file for writing
my $tmpFile=sprintf("%s.tmp$fileSuffix", $peakFile);
open(OUTFILE, ">$tmpFile") || die $!;

my %bInfo=(); my @t=(); my %NFR=();
for(my $i=0; $i<=scalar(@data)-2; $i++) {
    my @F1=split(/\s+/, $data[$i]);
    my @F2=split(/\s+/, $data[$i+1]);

    if($F1[3]=~/^$F2[3]$/ || ($F2[1]-$F1[2])<$maxNFRLength) {
        #print "DEBUG: $data[$i]\nDEBUG: $data[$i+1]\n";
        ## collect first block group information
        $bInfo{'first'}{'chr'}=$F1[0];
        $bInfo{'first'}{'start'}=$F1[1];
        $bInfo{'first'}{'end'}=$F1[2];
        $bInfo{'first'}{'strand'}=$F1[5];
        $bInfo{'first'}{'length'}=($bInfo{'first'}{'end'}-$bInfo{'first'}{'start'})+1;
        #@t=split(/\,/,$F1[4]);
        #$bInfo{'first'}{'expr'}=0;
        #foreach(@t) {
        #    $bInfo{'first'}{'expr'}+=$_;
        #}
        $coor="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
        $bInfo{'first'}{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
        chomp($bInfo{'first'}{'expr'});

        ## collect second block group information
        $bInfo{'second'}{'chr'}=$F2[0];
        $bInfo{'second'}{'start'}=$F2[1];
        $bInfo{'second'}{'end'}=$F2[2];
        $bInfo{'second'}{'strand'}=$F2[5];
        $bInfo{'second'}{'length'}=($bInfo{'second'}{'end'}-$bInfo{'second'}{'start'})+1;
        #@t=split(/\,/,$F2[4]);
        #$bInfo{'second'}{'expr'}=0;
        #foreach(@t) {
        #    $bInfo{'second'}{'expr'}+=$_;
        #}
        $coor="$bInfo{'second'}{'chr'}:$bInfo{'second'}{'start'}-$bInfo{'second'}{'end'}";
        $bInfo{'second'}{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
        chomp($bInfo{'second'}{'expr'});

        ## define NFR based on first and second block group information
        if($bInfo{'first'}{'strand'}=~/\+/ && $bInfo{'second'}{'strand'}=~/\+/) {
            $NFR{'chr'}=$bInfo{'first'}{'chr'};
            $NFR{'start'}=$bInfo{'first'}{'end'}+1;
            $NFR{'end'}=$bInfo{'second'}{'start'}-1;
            $NFR{'strand'}="+";
            $NFR{'length'}=($NFR{'end'}-$NFR{'start'})+1;
            $coor="$NFR{'chr'}:$NFR{'start'}-$NFR{'end'}";
            $NFR{'startBlock'}="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
            $NFR{'startBlockExpr'}=$bInfo{'first'}{'expr'};
            $NFR{'startBlockLen'}=$bInfo{'first'}{'length'};
            $NFR{'endBlock'}="$bInfo{'second'}{'chr'}:$bInfo{'second'}{'start'}-$bInfo{'second'}{'end'}";
            $NFR{'endBlockExpr'}=$bInfo{'second'}{'expr'};
            $NFR{'endBlockLen'}=$bInfo{'second'}{'length'};
            if($F1[3]=~/^$F2[3]$/) {
                $NFR{'id'}=$F1[3];
            }
            else {
                my @t1=split(/[\:\-]+/,$F1[3]);
                my @t2=split(/[\:\-]+/,$F2[3]);
                $NFR{'id'}="$t1[0]:$t1[1]-$t2[2]";
            }

            #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
            #print "DEBUG: $bInfo{'first'}{'expr'}\t$bInfo{'second'}{'expr'}\t$NFR{'expr'}\t$NFR{'length'}\n";

            $NFR{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
            chomp($NFR{'expr'});
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
            $NFR{'startBlockLen'}=$bInfo{'second'}{'length'};
            $NFR{'endBlock'}="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
            $NFR{'endBlockExpr'}=$bInfo{'first'}{'expr'};
            $NFR{'endBlockLen'}=$bInfo{'first'}{'length'};
            $NFR{'id'}=$F1[3];

            #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
            #print "DEBUG: $bInfo{'first'}{'expr'}\t$bInfo{'second'}{'expr'}\t$NFR{'expr'}\t$NFR{'length'}\n";

            $NFR{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
            chomp($NFR{'expr'});
        }
        else {
            print STDERR "WARNING: the first and second block do not have same strand despite having same ID\n";
            print STDERR "--> first block: $data[$i]\n";
            print STDERR "--> second block: $data[$i+1]\n";
        }

        if($NFR{'length'} >= $minNFRLength) {
            ## old scoring and output scheme
            #$NFR{'stddev'}=stddev($bInfo{'first'}{'expr'}, $bInfo{'second'}{'expr'});
            #$NFR{'stddev'}=~s/\,//g;
            #print("DEBUG: (($bInfo{'first'}{'expr'}/$bInfo{'first'}{'length'})+($bInfo{'second'}{'expr'}/$bInfo{'second'}{'length'}))/($NFR{'expr'}/$NFR{'length'}))\n");
            #$NFR{'score'}=((($NFR{'startBlockExpr'}/$NFR{'startBlockLen'})+($NFR{'endBlockExpr'}/$NFR{'endBlockLen'}))/($NFR{'expr'}/$NFR{'length'}));
            #$NFR{'score'}=sprintf("%0.4f", ($NFR{'score'}));
            #$NFR{'score'}=sprintf("%0.4f", log($NFR{'score'}));

            #print OUTFILE "$NFR{'chr'}\t$NFR{'start'}\t$NFR{'end'}\t$NFR{'id'}\t$NFR{'score'}\t$NFR{'strand'}\t$NFR{'startBlock'}\t$NFR{'startBlockExpr'}\t$NFR{'endBlock'}\t$NFR{'endBlockExpr'}\t$NFR{'expr'}\t$NFR{'length'}\n";

            ## new scoring and output scheme
            $NFR{'score'}=((($NFR{'startBlockExpr'}+$NFR{'endBlockExpr'})/($NFR{'startBlockLen'}+$NFR{'endBlockLen'}))-($NFR{'expr'}/$NFR{'length'}));
            $NFR{'score'}=sprintf("%0.4f", ($NFR{'score'}));

            print OUTFILE "$NFR{'chr'}\t$NFR{'start'}\t$NFR{'end'}\t$NFR{'id'}\t$NFR{'score'}\t$NFR{'strand'}\t$NFR{'startBlock'}\t$NFR{'endBlock'}\t$NFR{'startBlockExpr'}\t$NFR{'startBlockLen'}\t$NFR{'endBlockExpr'}\t$NFR{'endBlockLen'}\t$NFR{'expr'}\t$NFR{'length'}\n";
        }
    }
}
close OUTFILE;

## Step-3: define non-overlapping and significant nuclesome free regions
if(-e "$tmpFile") {
    @data=`zless $tmpFile | sort -k 1,1 -k 2n,2 -k 3n,3`;
}
else {
    print STDERR "Cannot find $tmpFile file\n";
    usage()
}

## create output file for writing non-overlapping NFR
if(scalar(@data)>=2) {
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
            print "$NFR\n";

            @F1=();
        }
    }
}
else {
    system("cat $tmpFile");
}

## remove temporary file
system("rm $tmpFile");

exit(0);
