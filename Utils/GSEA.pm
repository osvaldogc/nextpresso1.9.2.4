#!/usr/bin/perl

# Author: Osvaldo Grana
# Description: runs GSEA on selected gene sets
# v0.1		30Sep2014

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";


package GSEA;

sub runGSEA(){
	my($collapse,$mode,$norm,$nperm,$scoring_scheme,$include_only_symbols,$make_sets,
	$plot_top_x,$rnd_seed,$set_max,$set_min,$zip_report,$GSEAoutDir,$inputFile,$comparisonFile,$gseaPath,$gseaChip,$gseamaxMemory,$comparison,$tmpDir)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;

	my $command="export TMPDIR=".$tmpDir." && java -cp .:".$gseaPath." -Xmx".$gseamaxMemory." xtools.gsea.GseaPreranked ";
	$command.="-gmx ".$inputFile." ";
	$command.="-collapse ".$collapse." ";
	$command.="-mode ".$mode." ";
	$command.="-norm ".$norm." ";
	$command.="-nperm ".$nperm." ";
	
	$comparisonFile=~ s/-/_/;
	$comparisonFile=~ s/-/_/;
	$command.="-rnk ".$comparisonFile." ";
	$command.="-scoring_scheme ".$scoring_scheme." ";
	
	my @tokens=split('/',$inputFile);
	my $geneSetName=$tokens[@tokens-1];
	$geneSetName=~ s/-/_/;
	$geneSetName=~ s/ /_/;
	
	$command.="-rpt_label ".$comparison."_".$geneSetName." ";
	$command.="-chip ".$gseaChip." ";
	$command.="-include_only_symbols ".$include_only_symbols." ";
	$command.="-make_sets ".$make_sets." ";
	$command.="-plot_top_x ".$plot_top_x." ";
	$command.="-rnd_seed ".$rnd_seed." ";
	$command.="-set_max ".$set_max." ";
	$command.="-set_min ".$set_min." ";
	#$command.="-zip_report ".$zip_report." ";
	$command.="-out ".$GSEAoutDir." ";
	$command.="-gui true ";
	
	print "\n\t[executing] ".$command."\n";	
	system($command);	
}

1
