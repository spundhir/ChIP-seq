#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use perlModule;

###############################################################################
## parse input options
use vars qw($bedFile $percentOverlap $closest $scoreCol $help);
$scoreCol=5;

GetOptions ("i=s"  => \$bedFile,
            "v=s"  => \$percentOverlap,
            "c"    => \$closest,
            "l=s"  => \$scoreCol,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

###############################################################################
sub usage {
	print STDERR "\nProgram: mergeBed.pl (pick coordinate in BED format by selecting the highest scoring among all the overlapping coordinates)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: mergeBed.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>         [BED file (can be stdin, use -)]\n";
    print STDERR "                   [if STDIN, sort them by chrom start and end]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -v <int>          [define overlapping only if overlap is >= input percentage]\n";
    print STDERR " -c                [instead of highest scoring, choose the closest (useful for enhancer to promoter association)]\n";
    print STDERR " -l <int>          [column in input file containing score information (default: 5)]\n";
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

if(defined($closest)) {
    my $i=(); my @F=();
    my $col_start=(); my $col_end=(); my $col_strand=();

    ## determine file column containing start, end and strand information for gene
    for($i=0; $i<scalar(@data); $i++) {
        last if($data[$i]!~/\s+\.\s+/);
    }
    @F=split(/\s+/,$data[$i]);
    for($i=6; $i<=scalar(@F); $i++) {
        if($F[$i]=~/^chr/) {
            $col_start=$i+1;
            $col_end=$i+2;
        }
        elsif($F[$i]=~/^\+$/ || $F[$i]=~/^\-$/) {
            $col_strand=$i;
            last;
        }
    }
    if(!defined($col_strand)) {
        print STDERR "Cannot find strand information for gene in input BED file\n";
        exit(-1);
    }

    print "$col_start\t$col_end\t$col_strand\n"; exit;
    ## determine closest gene coordinate for each enhancer
    my $key=(); my $key_previous=();
    my $dist=(); my %enhancer=();
    for($i=0; $i<(scalar(@data)); $i++) {
        chomp($data[$i]);
        ## next, if no gene is found proximal to the enhancer
        #if($data[$i]=~/\s+\.\s+/) { 
        #    print "$data[$i]\n";
        #    next;
        #}
        #print "\t>>$data[$i]\n";
        ## compute, if gene(s) are found proximal to the enhancer
        @F=split(/\s+/, $data[$i]);
        $key="$F[0]_$F[1]_$F[2]";
        if($F[$col_strand]=~/\+/) { $dist=$F[$col_end]-$F[2]; }
        else { $dist=$F[1]-$F[$col_start]; }
        if(!defined($key_previous)) {
            #print "\t>>1\t$key\n";
            $enhancer{$key}{'line'}=$data[$i];
            $enhancer{$key}{'dist'}=$dist;
            $key_previous=$key;
        }
        elsif(!defined($enhancer{$key})) {
            #print "\t>>2\t$key\n";
            print "$enhancer{$key_previous}{'line'}\n";
            $enhancer{$key}{'line'}=$data[$i];
            $enhancer{$key}{'dist'}=$dist;
            $key_previous=$key;
        }
        else {
            if($dist < $enhancer{$key_previous}{'dist'}) {
                #print "\t>>3\t$key\n";
                #print "\t>>$enhancer{$key}{'dist'}\t$enhancer{$key_previous}{'dist'}\n";
                $enhancer{$key}{'line'}=$data[$i];
                $enhancer{$key}{'dist'}=$dist;
                $key_previous=$key;
            }
            #else {
                #print "\t>>4\t$key\n";
                #print "\t>>$enhancer{$key}{'dist'}\t$enhancer{$key_previous}{'dist'}\n";
            #}
        }
    }

    ## print last coordinate
    print "$enhancer{$key_previous}{'line'}\n";
}
elsif($percentOverlap) {
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
            if($F1[$scoreCol-1]<$F2[$scoreCol-1]) {
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
            #print $F1[$scoreCol-1]."\t".$F2[$scoreCol-1]."\n";
            if($F1[$scoreCol-1]<$F2[$scoreCol-1]) {
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
