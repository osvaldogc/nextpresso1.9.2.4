#!/usr/bin/perl -w

# Author: Osvaldo Grana
# Description: runs tophatFusion to predict gene fusions
# v0.1		oct2014

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use fusions;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options



sub main();
sub runTopHatFusion($$$$$$$$$$$$$$$$$$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;	

	my($numFusionReads,$numFusionPairs,$numFusionBoth,$fusionReadMismatches,$fusionMultireads);
	my($nonHuman,$nTophatFusionThreads,$pathToAnnotationFiles,$samples,$alignmentsDir);
	my($experimentName,$workspace,$tophatPath,$bowtiePath,$samtoolsPath,$pathToBlastAll,$pathToBlastn,$referenceSequence);

	undef($numFusionReads);
	undef($numFusionPairs);
	undef($numFusionBoth);
	undef($fusionReadMismatches);
	undef($fusionMultireads);
	undef($nonHuman);
	undef($nTophatFusionThreads);
	undef($pathToAnnotationFiles);
	undef($samples);
	undef($alignmentsDir);
	undef($experimentName);
	undef($workspace);
	undef($tophatPath);
	undef($bowtiePath);
	undef($samtoolsPath);
	undef($pathToBlastAll);
	undef($pathToBlastn);
	undef($referenceSequence);

	GetOptions(
		"numFusionReads=i"=>\$numFusionReads,#integer
		"numFusionPairs=i"=>\$numFusionPairs,	#integer
		"numFusionBoth=i"=>\$numFusionBoth, #integer
		"fusionReadMismatches=i"=>\$fusionReadMismatches, #integer
		"nTophatFusionThreads=i"=>\$nTophatFusionThreads,	#integer
		"fusionMultireads=i"=>\$fusionMultireads,	#integer
		"nonHuman=s"=>\$nonHuman,	#string
		"samples=s"=>\$samples,	#string
		"pathToAnnotationFiles=s"=>\$pathToAnnotationFiles,	#string
		"alignmentsDir=s"=>\$alignmentsDir,	#string
		"experimentName=s"=>\$experimentName,
		"workspace=s"=>\$workspace,
		"tophatPath=s"=>\$tophatPath,
		"bowtiePath=s"=>\$bowtiePath,
		"samtoolsPath=s"=>\$samtoolsPath,
		"pathToBlastAll=s"=>\$pathToBlastAll,
		"pathToBlastn=s"=>\$pathToBlastn,
		"referenceSequence=s"=>\$referenceSequence,
		
	);

#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 

	#*****IMPORTANT:
	#here $nThreads is no used as tophatFusion does it for all samples at once, hence
	#you don't specifiy sample names, and so there is no possibility of sending several threads at once
	runTopHatFusion($numFusionReads,$numFusionPairs,$numFusionBoth,$fusionReadMismatches,$fusionMultireads,$nonHuman,
	$nTophatFusionThreads,$pathToAnnotationFiles,$samples,$alignmentsDir,$experimentName,$workspace,
	$tophatPath,$bowtiePath,$samtoolsPath,$pathToBlastn,$pathToBlastAll,$referenceSequence);
		
}

sub runTopHatFusion($$$$$$$$$$$$$$$$$$){
	my($numFusionReads,$numFusionPairs,$numFusionBoth,$fusionReadMismatches,$fusionMultireads,$nonHuman,
	$nTophatFusionThreads,$pathToAnnotationFiles,$samples,$alignmentsDir,$experimentName,$workspace,
	$tophatPath,$bowtiePath,$samtoolsPath,$pathToBlastn,$pathToBlastAll,$referenceSequence)=@_;

	#1.creates symbolic links to the files required by tophatFusion	
	if($pathToAnnotationFiles!~ /\/$/){
		$pathToAnnotationFiles.=$pathToAnnotationFiles."/";
	}
	
	if(-e $pathToAnnotationFiles){
		if(-e $pathToAnnotationFiles."ensGene.txt"){
			system("cd ".$alignmentsDir."; ln -sf ".$pathToAnnotationFiles."ensGene.txt");
		}
		
		if(-e $pathToAnnotationFiles."mcl"){
			system("cd ".$alignmentsDir."; ln -sf ".$pathToAnnotationFiles."mcl");
		}
		
		if(-e $pathToAnnotationFiles."refGene.txt"){
			system("cd ".$alignmentsDir."; ln -sf ".$pathToAnnotationFiles."refGene.txt");
		}
		
		if(-e $pathToAnnotationFiles."blast"){
			system("cd ".$alignmentsDir."; ln -sf ".$pathToAnnotationFiles."blast");
		}
	
	
	}else{
		print STDERR "\n[ERROR fusions.pl]: $pathToAnnotationFiles link does not exist\n\n";
		exit(-1);
	}
	
	#2.creates symbolic links to the input files
	my @files=split('\,',$samples);
	my $fileNumber=0;	
	while(1){ #while there are files to be processed (threaded)
		my $inputFile=$files[$fileNumber]; #read files		

		my $symbolicLink=$alignmentsDir."tophat_".$inputFile;

		if(!-e $symbolicLink){
				system("cd ".$alignmentsDir."; ln -s ".$alignmentsDir.$inputFile." ".$symbolicLink);
				print STDERR "\n[fusions.pl]: $symbolicLink created\n";
		}
					
		$fileNumber++;		

		#finishes if all the files were threaded (sent to be processed)		
		if($fileNumber>=@files){
			last;
		}

	}#while(1)
	
	#3.launches tophatFusion
	&fusions::runTophatFusion($numFusionReads,$numFusionPairs,$numFusionBoth,$fusionReadMismatches,$fusionMultireads,$nonHuman,
	$nTophatFusionThreads,$alignmentsDir,$tophatPath,$bowtiePath,$samtoolsPath,$pathToBlastn,$pathToBlastAll,$referenceSequence);
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

