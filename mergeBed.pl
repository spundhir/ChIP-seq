#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;

###############################################################################
## parse input options
use vars qw($bedFile $percentOverlap $help);

GetOptions ("i=s"  => \$bedFile,
            "v=s"  => \$percentOverlap,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

###############################################################################
sub usage {
	print STDERR "\nProgram: mergeBed.pl (pick coordinate in BED format by selecting the highest scoring among all the overlapping coordinates)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: mergeBed.pl -i <file | STDIN> [OPTIONS]\n";
	print STDERR " -i <file>         [BED file]\n";
    print STDERR "                   [if STDIN, sort them by chrom start and end]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -v <int>          [define overlapping only if overlap is >= input percentage]\n";
	print STDERR " -h                [help]\n";
    print STDERR "[NOTE]\n";
    print STDERR " 1. input BED file must be sorted by chrom, then start and end\n";
    print STDERR " 2. the output coordinates differ when compared to mergeBed from BEDTools due to the reason that later progress by extending the end coordinate of overlapping coordinates leading to much more number of coordinates merged\n\n";
	exit(-1);
}

###############################################################################

my @data=();
if(defined($bedFile)) {
    $bedFile=~s/\,/ /g;
    chomp($bedFile);
    @data=`zless $bedFile | sortBed -i stdin`;
}
else {
    my $INFILE=();
    $INFILE=*STDIN;
    @data=<$INFILE>;
}

if($percentOverlap) {
    my @F1=(); my @F2=(); my $i=(); my $overlap=();
    for($i=0; $i<(scalar(@data)-1); $i++) {
        if(scalar(@F1)==0) {
            chomp($data[$i]);
            @F1=split(/\s+/, $data[$i]);
        }
        @F2=split(/\s+/, $data[$i+1]);

        ## check overlap between current and next coordinate
        $overlap=checkOverlap($F1[1], $F1[2], $F2[1], $F2[2], 5, $F1[0], $F2[0]);
        #print "$F1[1]\t$F1[2]\t$F2[1]\t$F2[2]\t$overlap\t$percentOverlap\n";

        ## choose the highest scoring one, if overlapped
        if($overlap>=$percentOverlap) {
            if($F1[4]<$F2[4]) {
                @F1=@F2; 
            }
        }
        else {
            ## print the non-overlapping and highest scoring coordinate
            my $COORDINATE=();
            foreach(@F1) { $COORDINATE.="$_\t"; } $COORDINATE=~s/\t$//g;
            print "$COORDINATE\n";
            @F1=();
        }
    }


    ## print for the last record
    if($overlap<$percentOverlap) {
        @F1=$data[$i];
    }

    ## print the non-overlapping and highest scoring coordinate
    my $COORDINATE=();
    foreach(@F1) { $COORDINATE.="$_\t"; } $COORDINATE=~s/\t$//g;
    print "$COORDINATE";
}
else {
    my @F1=(); my @F2=(); my $i=(); my $overlap=();
    for($i=0; $i<(scalar(@data)-1); $i++) {
        if(scalar(@F1)==0) {
            chomp($data[$i]);
            @F1=split(/\s+/, $data[$i]);
        }
        @F2=split(/\s+/, $data[$i+1]);

        ## check overlap between current and next coordinate
        $overlap=checkOverlap($F1[1], $F1[2], $F2[1], $F2[2], 4, $F1[0], $F2[0]);

        ## choose the highest scoring one, if overlapped
        if($overlap) {
            if($F1[4]<$F2[4]) {
                @F1=@F2; 
            }
        }
        else {
            ## print the non-overlapping and highest scoring coordinate
            my $COORDINATE=();
            foreach(@F1) { $COORDINATE.="$_\t"; } $COORDINATE=~s/\t$//g;
            print "$COORDINATE\n";
            @F1=();
        }
    }


    ## print for the last record
    if($overlap==0) {
        @F1=$data[$i];
    }

    ## print the non-overlapping and highest scoring coordinate
    my $COORDINATE=();
    foreach(@F1) { $COORDINATE.="$_\t"; } $COORDINATE=~s/\t$//g;
    print "$COORDINATE";
}
exit(0);
