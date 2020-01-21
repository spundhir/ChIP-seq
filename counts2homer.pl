#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

###############################################################################
## parse input options
use vars qw($inFile $db $motifName $help);
$db="NA";

GetOptions ("i=s"  => \$inFile,
            "j=s"  => \$db,
            "k=s"  => \$motifName,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

###############################################################################
sub usage {
	print STDERR "\nProgram: counts2homer.pl (convert from COUNTS to HOMER format)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: counts2homer.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>         [input file in NUCLEOTIDE COUNT format (can be stdin)]\n";
    print STDERR "[OPTIONS]\n";
    print STDERR " -j <string>       [database name (default: NA)]\n";
    print STDERR " -k <string>       [motif name (default: as defined by input matrix)]\n";
	print STDERR " -h                [help]\n";
	exit(-1);
}

###############################################################################

my @counts=();

if(!defined($inFile) || $inFile=~/^stdin$/) {
    my $INFILE=();
    $INFILE=*STDIN;
    @counts=<$INFILE>;
}
else {
    chomp($inFile);
    @counts=`zless $inFile`;
}


# conver from counts to probabilities in MEME format
my @data=(); my $prob=();
my @F=(); my $sum=();

push(@data, "MEME version 4\n");
push(@data, "ALPHABET= ACGT\n");
push(@data, "strands: + -\n");
push(@data, "Background letter frequencies (from uniform background):\n");
push(@data, "A 0.25000 C 0.25000 G 0.25000 T 0.25000");

foreach my $l(@counts) {
    my @F=split(/\s+/,$l);
    chomp($l);
    if($l=~/^>/) {
        push(@data, "");
        if(defined($motifName)) {
            push(@data, sprintf("MOTIF\t%s", $motifName));
        }
        else {
            $F[0]=~s/\>//g;
            push(@data, sprintf("MOTIF\t%s", $F[0]));
        }
        push(@data, "letter-probability matrix: alength= NA w= NA nsites= NA");
    }
    else {
        $sum=0;
        foreach(@F) {
            $sum+=$_;
        }
        $prob="";
        foreach(@F) {
            $prob.=sprintf("%0.6f\t", $_/$sum);
        }
        push(@data, $prob);
    }
}

#foreach(@data) { print "$_\n"; } exit;

## convert from probabilities in MEME format to homer matrix
my $id=(); my $des=();
my $start=(); my $score=();
my @freq=(); my @matrix=();
my $max=(); @F=();

foreach my $l(@data) {
    @F=split(/\s+/,$l);
    if($l=~/MOTIF/) {
        if(defined($motifName)) {
            $id=$motifName; $des="$motifName/$db";
        }
        elsif(defined($F[2])) {
            $id=$F[2]; $des="$F[2]/$db";
        }
        else {
            $id=$F[1]; $des="$F[1]/$db";
        }
    }
    elsif($l=~/letter-pro/) {
        $start=1; $score=0;
        @freq=(); @matrix=();
    }
    elsif($start && ($l=~/^$/ || $l=~/[a-zA-Z]+/)) {
        foreach(@freq) {
            $score+=log($_/0.25);
        }
        $score=$score-4;
        print ">$id\t$des\t$score\n";
        foreach(@matrix) {
            print "$_\n";
        }
        $start=0;
    }
    elsif($start) {
        $l=~s/^\s+//g;
        push(@matrix, $l);
        $max=0;
        my @T=split(/\s+/,$l);
        foreach(@T) {
            if($_>$max) {
                $max=$_;
            }
        }
        push(@freq, $max);
    }
}

## print output for the last motif read
if($start) {
    foreach(@freq) {
        $score+=log($_/0.25);
    }
    $score=$score-4;
    print ">$id\t$des\t$score\n";
    foreach(@matrix) {
        print "$_\n";
    }
    $start=0;
}
exit;
