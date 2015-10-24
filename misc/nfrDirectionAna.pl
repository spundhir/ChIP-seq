#!/usr/bin/perl -w

use strict;
use warnings;
use Tie::IxHash;
use Getopt::Long;
use Statistics::Basic qw(:all);
use perlModule;

###############################################################################
## parse input options
use vars qw($dirMe1 $dirMe3 $outDir $minNFRLength $maxNFRLength $nfrThreshold $genome $fileSuffix $help);

my %bamFile=();
my %sizeFactor=();
my %extendRead=();
$genome="mm9";
$fileSuffix="";

GetOptions ("i=s"  => \$dirMe1,
            "j=s"  => \$dirMe3,
            "o=s"  => \$outDir,
            "g=s"  => \$genome,
            "f=s"  => \$fileSuffix,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$dirMe1 || !$dirMe3 || !$outDir);

###############################################################################
sub usage {
	print STDERR "\nProgram: nfrDirectionAna.pl (perform direction analysis around NFRs using H3K4me1 and H3K4me3 based predictions)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: commonNFR.pl -i <dir> -j <dir> -o <dir> [OPTIONS]\n";
	print STDERR " -i <file>         [input directory containing NFR analysis results based on H3K4me1]\n";
	print STDERR " -j <file>         [input directory containing NFR analysis results based on H3K4me3]\n";
	print STDERR " -o <dir>          [directory where output files will be kept]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -g <string>       [genome (default: mm9)]\n";
    print STDERR " -f <string>       [a string added at the end of output files. useful when running in parallel]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

## create output directory, if does not exist
if ( ! -d $outDir) {
    system("mkdir $outDir");
}

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

readParameterFile($dirMe1, "me1");
readParameterFile($dirMe3, "me3");
=cu
foreach(keys(%bamFile)) {
    print "$_\t$bamFile{$_}{'rep1'}\n";
    print "$_\t$bamFile{$_}{'rep2'}\n";
    print "$_\t$extendRead{$_}{'rep1'}\n";
    print "$_\t$extendRead{$_}{'rep2'}\n";
    print "$_\t$sizeFactor{$_}{'rep1'}\n";
    print "$_\t$sizeFactor{$_}{'rep2'}\n";
    print "$minNFRLength\n";
    print "$maxNFRLength\n";
}
=cut

nfrDirectionAna("rep1");
#nfrDirectionAna("rep2");

sub nfrDirectionAna {
    my($rep)=@_;

    system("cat $dirMe1/nfr/$rep/*bg* > $outDir/h3k4me1.bg");
    system("cat $dirMe3/nfr/$rep/*bg* > $outDir/h3k4me3.bg");

    my @data=`closestBed -a $outDir/h3k4me1.bg -b $outDir/h3k4me3.bg -d -t first | perl -ane 'if(\$F[scalar(\@F)-1]<=3000) { print \$_; }'`;

    foreach my $l(@data) {
        chomp($l);
        my @F=split(/\s+/,$l);
        print "$l\n";

        my %coor=();
        ($coor{'left'}, $coor{'mid'}, $coor{'right'}, $coor{'order'})=organizeOverlapCoor($F[1], $F[2], $F[7], $F[8], $F[0], $F[6]);
        print("\t\t$coor{'left'}\t$coor{'mid'}\t$coor{'right'}\t$coor{'order'}\n");

        ## check if coordinates completely overlap
        next if($coor{'order'}==0);

        ## H3K4me1 block group is first
        if($coor{'order'}==1) {
            
        }
        ## H3K4me3 block group is first
        elsif($coor{'order'}==2) {
        
        }
    }
=cu
        $NFR{'overlapCoor'}=returnOverlapCoor($F[1], $F[2], $F[15], $F[16], $F[0], $F[14], "intersect", 0);
        $startBlock{'overlapCoor'}=returnOverlapCoor($startBlock{'rep1Coor'}[1], $startBlock{'rep1Coor'}[2], $startBlock{'rep2Coor'}[1], $startBlock{'rep2Coor'}[2], $startBlock{'rep1Coor'}[0], $startBlock{'rep2Coor'}[0], "union", 1);
        $endBlock{'overlapCoor'}=returnOverlapCoor($endBlock{'rep1Coor'}[1], $endBlock{'rep1Coor'}[2], $endBlock{'rep2Coor'}[1], $endBlock{'rep2Coor'}[2], $endBlock{'rep1Coor'}[0], $endBlock{'rep2Coor'}[0], "union", 2);

        ## check if overlap has been done correctly
        @{$NFR{'overlapCoorSplit'}}=split(/[\:\-]+/,$NFR{'overlapCoor'});
        @{$startBlock{'overlapCoorSplit'}}=split(/[\:\-]+/,$startBlock{'overlapCoor'});
        @{$endBlock{'overlapCoorSplit'}}=split(/[\:\-]+/,$endBlock{'overlapCoor'});

        if($NFR{'overlapCoorSplit'}[1]-$startBlock{'overlapCoorSplit'}[2]>1 || $endBlock{'overlapCoorSplit'}[1]-$NFR{'overlapCoorSplit'}[2]>1) {
            print STDERR "start or end block are not adjacent to NFR\n";
            print STDERR "--> $NFR{'overlapCoorSplit'}[1]\t$startBlock{'overlapCoorSplit'}[2]\t$endBlock{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\n";
            #print "$l\n";
            print STDERR "$F[0]:$F[1]-$F[2]\t$F[6]\t$F[7]\n";
            print STDERR "$F[14]:$F[15]-$F[16]\t$F[20]\t$F[21]\n";
            print STDERR "$NFR{'overlapCoor'}\t$startBlock{'overlapCoor'}\t$endBlock{'overlapCoor'}\n\n";
            exit(-1);
        }

        $NFR{'expr'}=`coor2expr -i $NFR{'overlapCoor'} -j $bamFileRep1,$bamFileRep2 -k $sizeFactorRep1,$sizeFactorRep2 -d -e $extendRep1,$extendRep2 -g $genome`;
        chomp($NFR{'expr'});
        $startBlock{'expr'}=`coor2expr -i $startBlock{'overlapCoor'} -j $bamFileRep1,$bamFileRep2 -k $sizeFactorRep1,$sizeFactorRep2 -d -e $extendRep1,$extendRep2 -g $genome`;
        chomp($startBlock{'expr'});
        $endBlock{'expr'}=`coor2expr -i $endBlock{'overlapCoor'} -j $bamFileRep1,$bamFileRep2 -k $sizeFactorRep1,$sizeFactorRep2 -d -e $extendRep1,$extendRep2 -g $genome`;
        chomp($endBlock{'expr'});

        #$NFR{'stddev'}=stddev($startBlock{'expr'}, $endBlock{'expr'});
        #$NFR{'stddev'}=~s/\,//g;

        $NFR{'length'}=($NFR{'overlapCoorSplit'}[2]-$NFR{'overlapCoorSplit'}[1])+1;
        $startBlock{'length'}=($startBlock{'overlapCoorSplit'}[2]-$startBlock{'overlapCoorSplit'}[1])+1;
        $endBlock{'length'}=($endBlock{'overlapCoorSplit'}[2]-$endBlock{'overlapCoorSplit'}[1])+1;

        ## old scoring scheme
        #$NFR{'score'}=((($startBlock{'expr'}/$startBlock{'length'})+($endBlock{'expr'}/$endBlock{'length'}))/($NFR{'expr'}/$NFR{'length'}));
        #$NFR{'score'}=sprintf("%0.4f", $NFR{'score'});
        #$NFR{'score'}=sprintf("%0.4f", log($NFR{'score'}));

        ## new scoring scheme
        $NFR{'score'}=((($startBlock{'expr'}+$endBlock{'expr'})/($startBlock{'length'}+$endBlock{'length'}))-($NFR{'expr'}/$NFR{'length'}));
        $NFR{'score'}=sprintf("%0.4f", $NFR{'score'});

        if($NFR{'length'} >= $minNFRLength) {
            print OUTFILE "$NFR{'overlapCoorSplit'}[0]\t$NFR{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\t$F[3]\t$NFR{'score'}\t$F[5]\t$startBlock{'overlapCoor'}\t$endBlock{'overlapCoor'}\t$startBlock{'expr'}\t$startBlock{'length'}\t$endBlock{'expr'}\t$endBlock{'length'}\t$NFR{'expr'}\t$NFR{'length'}\n";
            if(defined($nfrThreshold) && $NFR{'score'}>$nfrThreshold) {
                print SIGFILE "$NFR{'overlapCoorSplit'}[0]\t$NFR{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\t$F[3]\t$NFR{'score'}\t$F[5]\t$startBlock{'overlapCoor'}\t$endBlock{'overlapCoor'}\t$startBlock{'expr'}\t$startBlock{'length'}\t$endBlock{'expr'}\t$endBlock{'length'}\t$NFR{'expr'}\t$NFR{'length'}\n";
            }
        }
    }

    close OUTFILE;
    close SIGFILE if(defined($nfrThreshold));
=cut
}

sub readParameterFile {
    my($dir, $meth)=@_;

    ## read chosen parameters for H3K4me1 or me3 run
    if ( -e "$dir/nfr/PARAMETERS" ) {
        open(INFILE, "$dir/nfr/PARAMETERS") || die $!;
        foreach(<INFILE>) {
            chomp($_);
            if($_=~/^\#input BAM file \(Rep1\)\:/) {
                $_=~s/^.*\:\s+//g;
                $bamFile{$meth}{'rep1'}=$_;
            }
            elsif($_=~/^\#input BAM file \(Rep2\)\:/) {
                $_=~s/^.*\:\s+//g;
                $bamFile{$meth}{'rep2'}=$_;
            }
            elsif($_=~/^\#extend 3' end of reads \(Rep1\)\:/) {
                $_=~s/^.*\:\s+//g;
                $extendRead{$meth}{'rep1'}=$_;
            }
            elsif($_=~/^\#extend 3' end of reads \(Rep2\)\:/) {
                $_=~s/^.*\:\s+//g;
                $extendRead{$meth}{'rep2'}=$_;
            }
            elsif($_=~/^\#minimum length of NFR\:/) {
                $_=~s/^.*\:\s+//g;
                $minNFRLength=$_;
            }
            elsif($_=~/^\#maximum length of NFR\:/) {
                $_=~s/^.*\:\s+//g;
                $maxNFRLength=$_;
            }
        }
    } else {
        print STDERR "file $dir/nfr/PARAMETERS does not exist\n";
        exit(-1);
    }

    if ( -e "$dir/nfr/sizeFactor" ) {
        open(INFILE, "$dir/nfr/sizeFactor") || die $!;
        foreach(<INFILE>) {
            chomp($_);
            if($_=~/Rep1\s+/) {
                $_=~s/^.*\s+//g;
                $sizeFactor{$meth}{'rep1'}=$_;
            }
            elsif($_=~/Rep2\s+/) {
                $_=~s/^.*\s+//g;
                $sizeFactor{$meth}{'rep2'}=$_;
            }
        }
    }
    else {
        print STDERR "file $dir/nfr/sizeFactor does not exist\n";
        exit(-1);
    }
}

exit(0);
