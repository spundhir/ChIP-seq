#!/usr/bin/perl -w

use strict;
use warnings;
use Tie::IxHash;
use Getopt::Long;
use Statistics::Basic qw(:all);
use perlModule;

###############################################################################
## parse input options
use vars qw($bedFile $bedFileId $genome $threshold $fileSuffix $help);

my %bamFile=();
$bedFileId="NA";
$genome="mm9";
$threshold=0.25;
$fileSuffix="";

GetOptions ("i=s"  => \$bedFile,
            "a=s"  => \$bamFile{'me1'}{'rep1'},
            "b=s"  => \$bamFile{'me1'}{'rep2'},
            "c=s"  => \$bamFile{'me3'}{'rep1'},
            "d=s"  => \$bamFile{'me3'}{'rep2'},
            "g=s"  => \$genome,
            "t=s"  => \$threshold,
            "u=s"  => \$bedFileId,
            "f=s"  => \$fileSuffix,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$bedFile || !$bamFile{'me1'}{'rep1'} || !$bamFile{'me1'}{'rep2'} || !$bamFile{'me3'}{'rep1'} || !$bamFile{'me3'}{'rep2'});

###############################################################################
sub usage {
	print STDERR "\nProgram: bed2direction.pl (perform direction analysis around BED coordinates using H3K4me1 and H3K4me3 modifications)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: commonNFR.pl -i <file> -j <file> -k <file> -l <file> -m <file> -o <dir> [OPTIONS]\n";
	print STDERR " -i <file>         [input file in BED format]\n";
    print STDERR "                   [can be stdin using a '-']\n";
    print STDERR " -a <file>         [input BAM file corresponding to H3K4me1 (replicate 1)]\n";
    print STDERR " -b <file>         [input BAM file corresponding to H3K4me1 (replicate 2)]\n";
    print STDERR " -c <file>         [input BAM file corresponding to H3K4me3 (replicate 1)]\n";
    print STDERR " -d <file>         [input BAM file corresponding to H3K4me3 (replicate 2)]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -g <string>       [genome (default: mm9)]\n";
    print STDERR " -t <float>        [directionality threshold (default: 0.8)]\n";
    print STDERR " -u <string>       [unique identifier for input NFR file (default: NA)]\n";
    print STDERR " -f <string>       [a string added at the end of output files. useful when running in parallel]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

## populate genome file based on input genome
my $GENOME_FILE;
if($genome=~/^mm9$/) {
    $GENOME_FILE="/home/pundhir/project/genome_annotations/mouse.mm9.genome"
}
elsif($genome=~/^hg19$/) {
    $GENOME_FILE="/home/pundhir/project/genome_annotations/human.hg19.genome"
}
else {
    print STDERR "Presently the program only support analysis for mm9 or hg19\n";
    usage();
}

open(INFILE, "<$bedFile") || die $!;

my %coor=(); my %len=(); my %expr=(); my %D=();
foreach my $l(<INFILE>) {
    chomp($l);
    my @F=split(/\s+/, $l);
    my $mid=sprintf("%0.0f", ($F[1]+$F[2])/2);
    $coor{'nfr'}=sprintf("%s:%d-%d", $F[0], $mid-10, $mid+10);
    $coor{'upStream'}=sprintf("%s:%d-%d", $F[0], ($mid-10)-240, $mid-11);
    $coor{'downStream'}=sprintf("%s:%d-%d", $F[0], $mid+11, ($mid+10)+240);
    #$coor{'nfr'}="$F[0]:$F[1]-$F[2]";
    #$coor{'upStream'}=$F[6];
    #$coor{'downStream'}=$F[7];
    #$len{'nfr'}=$F[13];
    #$len{'upStream'}=$F[9];
    #$len{'downStream'}=$F[11];

    ## H3K4me1 expression
    $expr{'me1'}{'nfr'}=`coor2expr -i $coor{'nfr'} -j $bamFile{'me1'}{'rep1'},$bamFile{'me1'}{'rep2'} -m -d -g $genome -p 0`;
    $expr{'me1'}{'upStream'}=`coor2expr -i $coor{'upStream'} -j $bamFile{'me1'}{'rep1'},$bamFile{'me1'}{'rep2'} -m -d -g $genome -p 0`;
    $expr{'me1'}{'downStream'}=`coor2expr -i $coor{'downStream'} -j $bamFile{'me1'}{'rep1'},$bamFile{'me1'}{'rep2'} -m -d -g $genome -p 0`;
    $expr{'me1'}{'total'}=$expr{'me1'}{'nfr'}+$expr{'me1'}{'upStream'}+$expr{'me1'}{'downStream'};
    chomp($expr{'me1'}{'nfr'});
    chomp($expr{'me1'}{'upStream'});
    chomp($expr{'me1'}{'downStream'});
    chomp($expr{'me1'}{'total'});
    $expr{'me1'}{'nfr'}=$expr{'me1'}{'nfr'}+0.01;
    $expr{'me1'}{'upStream'}=$expr{'me1'}{'upStream'}+0.01;
    $expr{'me1'}{'downStream'}=$expr{'me1'}{'downStream'}+0.01;
    $expr{'me1'}{'all'}=$expr{'me1'}{'nfr'}+$expr{'me1'}{'upStream'}+$expr{'me1'}{'downStream'}
    #$expr{'me1'}{'nfr'}=sprintf("%0.2f", $expr{'me1'}{'nfr'}/$len{'nfr'});
    #$expr{'me1'}{'upStream'}=sprintf("%0.2f", $expr{'me1'}{'upStream'}/$len{'upStream'});
    #$expr{'me1'}{'downStream'}=sprintf("%0.2f", $expr{'me1'}{'downStream'}/$len{'downStream'});

    ## H3K4me3 expression
    $expr{'me3'}{'nfr'}=`coor2expr -i $coor{'nfr'} -j $bamFile{'me3'}{'rep1'},$bamFile{'me3'}{'rep2'} -m -d -g $genome -p 0`;
    $expr{'me3'}{'upStream'}=`coor2expr -i $coor{'upStream'} -j $bamFile{'me3'}{'rep1'},$bamFile{'me3'}{'rep2'} -m -d -g $genome -p 0`;
    $expr{'me3'}{'downStream'}=`coor2expr -i $coor{'downStream'} -j $bamFile{'me3'}{'rep1'},$bamFile{'me3'}{'rep2'} -m -d -g $genome -p 0`;
    $expr{'me3'}{'total'}=$expr{'me3'}{'nfr'}+$expr{'me3'}{'upStream'}+$expr{'me3'}{'downStream'};
    chomp($expr{'me3'}{'nfr'});
    chomp($expr{'me3'}{'upStream'});
    chomp($expr{'me3'}{'downStream'});
    chomp($expr{'me3'}{'total'});
    $expr{'me3'}{'nfr'}=$expr{'me3'}{'nfr'}+0.01;
    $expr{'me3'}{'upStream'}=$expr{'me3'}{'upStream'}+0.01;
    $expr{'me3'}{'downStream'}=$expr{'me3'}{'downStream'}+0.01;
    $expr{'me3'}{'all'}=$expr{'me3'}{'nfr'}+$expr{'me3'}{'upStream'}+$expr{'me3'}{'downStream'}
    #$expr{'me3'}{'nfr'}=sprintf("%0.2f", $expr{'me3'}{'nfr'}/$len{'nfr'});
    #$expr{'me3'}{'upStream'}=sprintf("%0.2f", $expr{'me3'}{'upStream'}/$len{'upStream'});
    #$expr{'me3'}{'downStream'}=sprintf("%0.2f", $expr{'me3'}{'downStream'}/$len{'downStream'});

    ## H3K4me1 directionality score
    $D{'me1'}=($expr{'me1'}{'upStream'}-$expr{'me1'}{'downStream'})/($expr{'me1'}{'upStream'}+$expr{'me1'}{'downStream'});
    $D{'me1'}=sprintf("%0.2f", $D{'me1'});

    ## H3K4me3 directionality score
    $D{'me3'}=($expr{'me3'}{'upStream'}-$expr{'me3'}{'downStream'})/($expr{'me3'}{'upStream'}+$expr{'me3'}{'downStream'});
    $D{'me3'}=sprintf("%0.2f", $D{'me3'});
    if($D{'me1'} > $threshold && $D{'me3'} < -1*$threshold) {
        print "$l\t$bedFileId\t$expr{'me1'}{'upStream'}\t$expr{'me1'}{'nfr'}\t$expr{'me1'}{'downStream'}\t$D{'me1'}\t$expr{'me3'}{'upStream'}\t$expr{'me3'}{'nfr'}\t$expr{'me3'}{'downStream'}\t$D{'me3'}\tUS\n";
    }
    elsif($D{'me1'} < -1*$threshold && $D{'me3'} > $threshold) {
        print "$l\t$bedFileId\t$expr{'me1'}{'upStream'}\t$expr{'me1'}{'nfr'}\t$expr{'me1'}{'downStream'}\t$D{'me1'}\t$expr{'me3'}{'upStream'}\t$expr{'me3'}{'nfr'}\t$expr{'me3'}{'downStream'}\t$D{'me3'}\tSU\n";
    }
    elsif($expr{'me1'}{'total'} > $expr{'me3'}{'total'}) {
        print "$l\t$bedFileId\t$expr{'me1'}{'upStream'}\t$expr{'me1'}{'nfr'}\t$expr{'me1'}{'downStream'}\t$D{'me1'}\t$expr{'me3'}{'upStream'}\t$expr{'me3'}{'nfr'}\t$expr{'me3'}{'downStream'}\t$D{'me3'}\tUU\n";
    }
    else {
        print "$l\t$bedFileId\t$expr{'me1'}{'upStream'}\t$expr{'me1'}{'nfr'}\t$expr{'me1'}{'downStream'}\t$D{'me1'}\t$expr{'me3'}{'upStream'}\t$expr{'me3'}{'nfr'}\t$expr{'me3'}{'downStream'}\t$D{'me3'}\tSS\n";
    }
    #print "\t\t$coor{'upStream'}\t$coor{'nfr'}\t$coor{'downStream'}\n";
}
exit;
