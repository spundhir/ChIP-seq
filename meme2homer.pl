#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

###############################################################################
## parse input options
use vars qw($inFile $db $motifName $extractMotif $help);
$db="NA";

GetOptions ("i=s"  => \$inFile,
            "j=s"  => \$db,
            "k=s"  => \$motifName,
            "e=s"  => \$extractMotif,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

###############################################################################
sub usage {
	print STDERR "\nProgram: meme2homer.pl (convert from MEME to HOMER format)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: meme2homer.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>         [input file in MEME format (can be stdin)]\n";
    print STDERR "[OPTIONS]\n";
    print STDERR " -j <string>       [database name (default: NA)]\n";
    print STDERR " -k <string>       [motif name (default: as defined by MEME)]\n";
    print STDERR " -e <string>       [only extract the matrix for given motif name]\n";
	print STDERR " -h                [help]\n";
	exit(-1);
}

###############################################################################

my @data=();

if(!defined($inFile) || $inFile=~/^stdin$/) {
    my $INFILE=();
    $INFILE=*STDIN;
    @data=<$INFILE>;
}
else {
    chomp($inFile);
    @data=`zless $inFile`;
}

my $id=(); my $des=();
my $start=(); my $score=();
my @freq=(); my @matrix=();
my $max=(); my @F=();

if(defined($extractMotif)) {
    print "MEME version 4\n\n";
    print "ALPHABET= ACGT\n\n";
    print "strands: + -\n\n";
    print "Background letter frequencies (from uniform background):\n";
    print "A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n";
    $start=0;
    foreach my $l(@data) {
        @F=split(/\s+/,$l);
        if($l=~/MOTIF/ && $F[2]=~/$extractMotif/) {
            print "$l";
            $start=1;
        } elsif($start==1 && $l=~/^URL/) {
            print "$l";
            $start=0;
            last;
        } elsif($start) {
            print "$l";
        }
    }
} else {
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
                print $_;
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
            print $_;
        }
        $start=0;
    }
}
exit;
