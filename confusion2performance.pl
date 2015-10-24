#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Tie::IxHash;

##################################################################################################
# parse input options
use vars qw($help $tp $fp $tn $fn);

GetOptions ("i=i"  => \$tp,
            "j=i"  => \$tn,
            "k=i"  => \$fp,
            "l=i"  => \$fn,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !defined($tp) || !defined($tn) || !defined($fp) || !defined($fn));

#################################################################################################
sub usage {
	print STDERR "\nProgram: confusion2performance.pl (given a confusion matrix, evaluate the performance measures)\n";
	print STDERR "Author: RTH, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: confusion2performance.pl -i <int> -j <int> -k <int> -l <int>\n";
	print STDERR " -i <int>    [number of true positives]\n";
	print STDERR " -j <int>    [number of true negatives]\n";
	print STDERR " -k <int>    [number of false positives]\n";
	print STDERR " -l <int>    [number of false negatives]\n";
	print STDERR " -h          [this useful help message]\n";
	print STDERR "\n";
	exit(-1);
}
##################################################################################################

my %performance=();
tie %performance, 'Tie::IxHash';

$performance{'sensitivity'}=($tp)/($tp+$fn);
$performance{'specificity'}=($tn)/($tn+$fp);
if(sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn))>0) {
	$performance{'mcc'}=(($tp*$tn)-($fp*$fn))/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn));
}
else {
	$performance{'mcc'}="inf";
}
$performance{'tpr'}=($tp)/($tp+$fn);
$performance{'fpr'}=1-$performance{'specificity'};

foreach(keys(%performance)) {
	printf("%s\t%0.2f\n", $_, $performance{$_});
}

exit;
