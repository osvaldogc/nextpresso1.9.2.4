#!/usr/bin/perl
# FileName : preprocess.pm
# Author : Osvaldo Grana
# Description: performs trimming / downsampling of reads in fastq files
# v0.1		9abr2014
# v0.2		oct2017, removed ulimit filtering

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";

package preprocess;

sub trimming($$$$$){
	my($seqtkPath,$nNucleotidesLeftEnd,$nNucleotidesRightEnd,$inputFile,$outputFile)=@_;

	my $command=$seqtkPath."seqtk trimfq -b ".$nNucleotidesLeftEnd." -e ".$nNucleotidesRightEnd." ".$inputFile." > ".$outputFile;
	print "\n\t[executing] ".$command."\n";	
	system($command);

}

sub downsampling($$$$$){
	my($seqtkPath,$seed,$nReads,$inputFile,$outputFile)=@_;

	#my $command="ulimit -v 2000000; ".$seqtkPath."seqtk sample -s".$seed." ".$inputFile." ".$nReads." > ".$outputFile;
	my $command=$seqtkPath."seqtk sample -s".$seed." ".$inputFile." ".$nReads." > ".$outputFile;
	print "\n\t[executing] ".$command."\n";	
	system($command);

}

sub fromBamToFastq($$$){
	my($bedtoolsPath,$inputFile,$outputFile)=@_;
	
	my $command=$bedtoolsPath."bamToFastq -i ".$inputFile." -fq ".$outputFile;
	print "\n\t[executing] ".$command."\n";	
	system($command);
}


1
