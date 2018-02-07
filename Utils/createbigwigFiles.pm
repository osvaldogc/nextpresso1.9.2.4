#!/usr/bin/perl

# Author : Osvaldo Grana
# Description: creates bigWigFiles from bam files
# v0.1		Aug2014

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";


package createbigwigFiles;

sub runSamtools(){
	my($samtoolsPath,$bamFile,$sortedBamFile)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	my $command="export PATH=\$PATH:".$samtoolsPath."; ";
	$command.="samtools sort ".$bamFile." ".$sortedBamFile;
	
	print "\n\t[executing] ".$command."\n";	
	system($command);	
}

sub runGenomeCoverageBed(){
	my($chromosomeSizesFile,$bedtoolsPath,$sortedBamFile,$bedGraphFile,$bedGraphFileHeader)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;

	my $command="export PATH=\$PATH:".$bedtoolsPath."; ";
	$command.="genomeCoverageBed -split -trackline -trackopts ".$bedGraphFileHeader." -bga -ibam ".$sortedBamFile." -g ".$chromosomeSizesFile." > ".$bedGraphFile;

	print "\n\t[executing] ".$command."\n";	
	system($command);
	
}

sub runBedGraphToBigWig(){
	my($chromosomeSizesFile,$bedGraphToBigWigPath,$bedGraphFile,$bigWigFile)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;

	my $command="export PATH=\$PATH:".$bedGraphToBigWigPath."; ";
	$command.="bedGraphToBigWig ".$bedGraphFile." ".$chromosomeSizesFile." ".$bigWigFile;

	print "\n\t[executing] ".$command."\n";	
	system($command);
	
}




1
