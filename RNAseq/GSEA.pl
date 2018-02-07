#!/usr/bin/perl -w

# Author: Osvaldo Grana
# Description: runs GSEA on selected gene sets
# v0.1		sep2014
# v0.2		dic2016 adds a method to run GSEA on DESeq2 preranked files

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use GSEA;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;


sub main();
sub runGSEA($$$$$$$$$$$$$$$$$$$$$);
sub checkThreads($$$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;

	my($collapse,$mode,$norm,$nperm,$scoring_scheme,$include_only_symbols,$make_sets,
	$plot_top_x,$rnd_seed,$set_max,$set_min,$zip_report,$GSEAoutDir,$genesets,$nThreads,$comparison,$diffExp_outDir,$gseaPath,$gseaChip,$gseamaxMemory,$tmpDir);

	undef($collapse);
	undef($mode);
	undef($norm);
	undef($nperm);
	undef($scoring_scheme);
	undef($include_only_symbols);
	undef($make_sets);
	undef($plot_top_x);
	undef($rnd_seed);
	undef($set_max);
	undef($set_min);
	undef($zip_report);
	undef($GSEAoutDir);
	undef($genesets);
	undef($comparison);
	undef($diffExp_outDir);
	undef($gseaPath);
	undef($gseaChip);
	undef($gseamaxMemory);
	undef($tmpDir);
	$nThreads=2;

	GetOptions(
		"mode=s"=>\$mode,
		"collapse=s"=>\$collapse,
		"norm=s"=>\$norm, #string
		"nperm=s"=>\$nperm, #string
		"numberOfThreads=i"=>\$nThreads,	#integer
		"scoring_scheme=s"=>\$scoring_scheme,	#string
		"include_only_symbols=s"=>\$include_only_symbols,	#string
		"make_sets=s"=>\$make_sets,	#string
		"plot_top_x=s"=>\$plot_top_x,	#string
		"rnd_seed=s"=>\$rnd_seed,	#string
		"set_max=s"=>\$set_max,	#string
		"set_min=s"=>\$set_min,
		"zip_report=s"=>\$zip_report,
		"GSEAoutDir=s"=>\$GSEAoutDir,
		"diffExp_outDir=s"=>\$diffExp_outDir,
		"genesets=s"=>\$genesets,
		"comparison=s"=>\$comparison,
		"gseaPath=s"=>\$gseaPath,
		"gseaChip=s"=>\$gseaChip,
		"gseamaxMemory=s"=>\$gseamaxMemory,
		"tmpDir=s"=>\$tmpDir	
		
	);

#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 

	runGSEA($collapse,$mode,$norm,$nperm,$scoring_scheme,$include_only_symbols,$make_sets,
	$plot_top_x,$rnd_seed,$set_max,$set_min,$zip_report,$GSEAoutDir,$genesets,$nThreads,$comparison,$diffExp_outDir,$gseaPath,$gseaChip,$gseamaxMemory,$tmpDir);
	
	
		
}

sub runGSEA($$$$$$$$$$$$$$$$$$$$$){	

	my($collapse,$mode,$norm,$nperm,$scoring_scheme,$include_only_symbols,$make_sets,
	$plot_top_x,$rnd_seed,$set_max,$set_min,$zip_report,$GSEAoutDir,$genesets,$nThreads,$comparison,$diffExp_outDir,$gseaPath,$gseaChip,$gseamaxMemory,$tmpDir)=@_;

	my @files=split(',',$genesets);
	#preranked file
	my $comparisonFile=$diffExp_outDir.$comparison.".rnk";
	
	#checks how many files there are to process	
	my $counter=@files;
	my $memory=$counter;
	
	#in case that the number of allowed threads exceeds the number of samples to process,
	#the number of allowed threads is adjusted to the number of samples
	if($nThreads>$memory){$nThreads=$memory}
	
	#creates an array of threads
	my @setOfThreads=();
	my @hasThisThreadFinishedExecution=();
	for(my $i=0;$i<$counter;$i++){
		undef $setOfThreads[$i];
		$hasThisThreadFinishedExecution[$i]=0;
	}

	$counter=0; #counter that accounts for the number of threads that are currently active
	my $fileNumber=0; #takes care of the current file that must be processed (threaded)
	my $currentThread=0;
	
	while(1){ #while there are files to be processed (threaded)
		my $inputFile=$files[$fileNumber]; #read files		
		
		#check files
		if(!-e $inputFile){
				print STDERR "\n[ERROR GSEA.pl]: $inputFile file does not exist\n\n";
				exit(-1);
		}				

		if(!-e $comparisonFile){
				print STDERR "\n[ERROR GSEA.pl]: $comparisonFile file does not exist\n\n";
				exit(-1);
		}
						
		$setOfThreads[$currentThread]=threads->create(\&GSEA::runGSEA,
		$collapse,$mode,$norm,$nperm,$scoring_scheme,$include_only_symbols,$make_sets,
		$plot_top_x,$rnd_seed,$set_max,$set_min,$zip_report,$GSEAoutDir,$inputFile,$comparisonFile,$gseaPath,$gseaChip,$gseamaxMemory,$comparison,$tmpDir);		
		
		$currentThread++;
		$fileNumber++;		
		$counter++;
		
		#checks if the number of active threads reached the maximum of threads allowed
		if($counter>=$nThreads){
			$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
		}					

		#finishes if all the files were threaded (sent to be processed)		
		if($fileNumber>=@files){
			last;
		}

	}#while(1)
	
	#finally, joins the threads
	$counter=0;
	my $joinedThreads=0;
	while(1){
		if($setOfThreads[$counter]->is_joinable()){
			$setOfThreads[$counter]->join();
			print "joined thread $counter\n";
			$joinedThreads++;
		}
			
		$counter++;
			
		if($counter>=@setOfThreads){$counter=0}
		if($joinedThreads>=@setOfThreads){last}			
	}		
}

sub checkThreads($$$){
#when the maximun number of allowed threads was reached (or run), it controls when one
#of them has finished. In that case it returns that a new thread can be launched
#(by decreasing the counter, allowing one more thread to be launched)
	
	my($threads,$counter,$finishedExecution)=@_;
	
	while(1){
		for(my $i=0;$i<@$threads;$i++){
			if(!$finishedExecution->[$i] && defined $threads->[$i] && !$threads->[$i]->is_running()){
				$finishedExecution->[$i]=1;
				return ($counter-1);
			}
		}
	}
}

sub help(){
	my $usage = qq{
--quality.pl--
			
	Execution error. Check out the execution modes in the following examples:
			
	Example for trimming a single end experiment
		perl preprocess.pl --nThreads 6 --mode 1 --seqtkPath /path --positions 0:3,0:3 --samples /path/file1.fq:/path/file1.trimmed.fq
				
	Example for trimming a paired end experiment
		perl preprocess.pl --nThreads 6 --mode 1 --seqtkPath /path --positions 0:3,0:3 --samples /path/file1.fq:/path/file1.trimmed.fq,/path/file2.fq:/path/file2.trimmed.fq
				
	Example for downsampling a single end experiment
		perl preprocess.pl --nThreads 6 --mode 2 --seqtkPath /path --seed 123L --nReads 10000 --samples /path/file1.fq,/path/file2.fq
				
	Example for downsampling a paired end experiment
		perl preprocess.pl --nThreads 6 --mode 2 --seqtkPath /path --seed 123L --nReads 10000 --samples /path/file1_R1.fq:/path/file1_R2.fq,/path/file2_R1.fq:/path/file2_R2.fq
				


};
 
	print STDERR $usage;
	exit(1);
}

