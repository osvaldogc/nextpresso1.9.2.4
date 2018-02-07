#!/usr/bin/perl

# Author : Osvaldo Grana
# Description: runs tophatFusion to predict gene fusions 
# v0.1		7Oct2014

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";


package fusions;

sub runTophatFusion(){
	my ($numFusionReads,$numFusionPairs,$numFusionBoth,$fusionReadMismatches,$fusionMultireads,$nonHuman,
	$nTophatFusionThreads,$alignmentsDir,$tophatPath,$bowtiePath,$samtoolsPath,$pathToBlastn,$pathToBlastAll,$referenceSequence)=@_;

	use Env qw(PATH);
	
	my $command="export PATH=\$PATH:".$tophatPath.":".$bowtiePath.":".$samtoolsPath.":".$pathToBlastn.":".$pathToBlastAll."; ";
	$command.="cd ".$alignmentsDir."; ";
	$command.="tophat-fusion-post ";
	$command.="-p ".$nTophatFusionThreads;
	$command.=" --num-fusion-reads ".$numFusionReads;
	$command.=" --num-fusion-pairs ".$numFusionPairs;
	$command.=" --num-fusion-both ".$numFusionBoth;
	$command.=" --fusion-read-mismatches ".$fusionReadMismatches;
	$command.=" --fusion-multireads ".$fusionMultireads;
	if($nonHuman eq "true"){$command.=" --non-human"}
	$referenceSequence=~ s/\.fa$//;
	$command.=" ".$referenceSequence;
	
	print "\n\t[executing] ".$command."\n";	
	system($command);

}


1
