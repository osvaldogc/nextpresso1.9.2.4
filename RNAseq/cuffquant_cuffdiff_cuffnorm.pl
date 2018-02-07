#!/usr/bin/perl -w

# Author: Osvaldo Grana
# Description: performs transcripts quantification and assembly
# v0.3		oct2016
# v0.4		dic2016 - adds colors to sample names in PCA plots

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use cuffquant_cuffdiff_cuffnorm;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;


sub main();
sub runCuffquant($$$$$$$$$$$$$$$$$);
sub runCuffdiff($$$$$$$$$$$$$$$$$$$$$$$);
sub runCuffnorm($$$$$$$$$$$$$$);
sub calculateCorrelationsAndPCA($$$$$);
sub checkThreads($$$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;
		
	my($alignmentsDir,$cufflinksPath,$samtoolsPath,$libraryType,$referenceSequence,$cuffnormLibraryNormalizationMethod,$extraPathsRequired);
	my($cuffnormNormalization,$cuffdiffLibraryNormalizationMethod,$cuffnormLibraryType,$cuffnormLabels,$cuffnormInputFiles);
	my($GTF,$samples,$cuffquantOutDir,$cuffdiffOutDir,$cuffnormOutDir,$cuffdiffLibraryType,$cuffdiffComparisonLabels,$cuffdiffInputFiles);
	my($cuffquant_noEffectiveLengthCorrection,$cuffquant_noLengthCorrection);
	my($cuffdiff_noEffectiveLengthCorrection,$cuffdiff_noLengthCorrection,$cuffdiff_dispersionMethod,$cuffnormColors,$executionCreatedTempDir);
	
	my $mode="123";
	my $nThreads=1;
	my $cuffquantNThreads=1;
	my $cuffnormNThreads=1;
	my $cuffnormOutputFormat="simple-table";
	my $cuffdiffNThreads=1;
	my $cuffdiffFDR=0.05;
	my $cuffdiffMinAlignmentCount=10;
	my $cuffnormSeed=123;
	my $cuffquantSeed=123;
	my $cuffdiffSeed=123;
	my $cuffquantFragBiasCorrect="true";
	my $cuffdiffFragBiasCorrect="true";
	my $cuffquantMultiReadCorrect="true";
	my $cuffdiffMultiReadCorrect="true";
	my $FPKMthreshold="0.05";
	my $cuffquantMaxBundleFrags=500000;
	my $cuffdiffMaxBundleFrags=500000;
	
	undef($cufflinksPath);
	undef($mode);
	undef($referenceSequence);
	undef($samples);
	undef($GTF);
	undef($cuffquantOutDir);
	undef($cuffdiffOutDir);
	undef($cuffnormOutDir);
	undef($alignmentsDir);
	undef($extraPathsRequired);
	undef($cuffquant_noEffectiveLengthCorrection);
	undef($cuffquant_noLengthCorrection);
	undef($cuffdiff_noEffectiveLengthCorrection);
	undef($cuffdiff_noLengthCorrection);
	undef($cuffdiff_dispersionMethod);
	undef($cuffnormColors);
	undef($executionCreatedTempDir);
	
	GetOptions(
		"mode=i"=>\$mode, #integer
		"nThreads=i"=>\$nThreads,	#integer
		"cufflinksPath=s"=>\$cufflinksPath, #string
		"samtoolsPath=s"=>\$samtoolsPath, #string
		"cuffquantNThreads=i"=>\$cuffquantNThreads,	#integer
		"cuffnormNThreads=i"=>\$cuffnormNThreads,	#integer
		"cuffdiffNThreads=i"=>\$cuffdiffNThreads,	#integer
		"cuffnormOutputFormat=s"=>\$cuffnormOutputFormat,	#string
		"cuffnormlibraryNormalizationMethod=s"=>\$cuffnormLibraryNormalizationMethod,	#string
		"cuffnormNormalization=s"=>\$cuffnormNormalization,	#string
		"cuffnormSeed=s"=>\$cuffnormSeed,	#string
		"samples=s"=>\$samples,	#string
		"cuffdiffFDR=s"=>\$cuffdiffFDR,	#string
		"cuffdiffMinAlignmentCount=s"=>\$cuffdiffMinAlignmentCount,	#string
		"cuffdiffFragBiasCorrect=s"=>\$cuffdiffFragBiasCorrect,	#string
		"cuffdiffMultiReadCorrect=s"=>\$cuffdiffMultiReadCorrect,	#string
		"cuffdiffLibraryNormalizationMethod=s"=>\$cuffdiffLibraryNormalizationMethod,	#string
		"libraryType=s"=>\$libraryType,	#string
		"GTF=s"=>\$GTF, #string
		"referenceSequence=s"=>\$referenceSequence, #string
		"alignmentsDir=s"=>\$alignmentsDir, #string
		"cuffquantOutDir=s"=>\$cuffquantOutDir, #string
		"cuffdiffOutDir=s"=>\$cuffdiffOutDir, #string
		"cuffnormOutDir=s"=>\$cuffnormOutDir, #string
		"cuffquantFragBiasCorrect=s"=>\$cuffquantFragBiasCorrect, #string
		"cuffquantMultiReadCorrect=s"=>\$cuffquantMultiReadCorrect, #string
		"cuffquantSeed=s"=>\$cuffquantSeed, #string
		"cuffdiffSeed=s"=>\$cuffdiffSeed,	#string
		"cuffdiffInputFiles=s"=>\$cuffdiffInputFiles,	#string
		"cuffnormInputFiles=s"=>\$cuffnormInputFiles,	#string
		"cuffdiffComparisonLabels=s"=>\$cuffdiffComparisonLabels,	#string
		"cuffnormLabels=s"=>\$cuffnormLabels,	#string
		"cuffdiffLibraryType=s"=>\$cuffdiffLibraryType,	#string
		"cuffnormLibraryType=s"=>\$cuffnormLibraryType,	#string
		"FPKMthreshold=s"=>\$FPKMthreshold,	#string
		"extraPathsRequired=s"=>\$extraPathsRequired,
		"cuffquantMaxBundleFrags=i"=>\$cuffquantMaxBundleFrags,
		"cuffdiffMaxBundleFrags=i"=>\$cuffdiffMaxBundleFrags,		
		"cuffquant_noEffectiveLengthCorrection=s"=>\$cuffquant_noEffectiveLengthCorrection,
		"cuffquant_noLengthCorrection=s"=>\$cuffquant_noLengthCorrection,		
		"cuffdiff_noEffectiveLengthCorrection=s"=>\$cuffdiff_noEffectiveLengthCorrection,
		"cuffdiff_noLengthCorrection=s"=>\$cuffdiff_noLengthCorrection,
		"cuffdiff_dispersionMethod=s"=>\$cuffdiff_dispersionMethod,
		"colors=s"=>\$cuffnormColors,
		"tmpDir=s"=>\$executionCreatedTempDir
	);

#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 
	
	if($mode=~ "1"){		
		runCuffquant($extraPathsRequired,$cufflinksPath,$samtoolsPath,$alignmentsDir,$referenceSequence,$samples,$GTF,$libraryType,$nThreads,
		$cuffquantNThreads,$cuffquantSeed,$cuffquantFragBiasCorrect,$cuffquantMultiReadCorrect,$cuffquantOutDir,$cuffquantMaxBundleFrags,
		$cuffquant_noEffectiveLengthCorrection,$cuffquant_noLengthCorrection);
	}
	if($mode=~ "2"){
		# replaces the initial separation between the samples from different conditions (represented by ':'), by a space (required in this way by cuffdiff)
		$cuffdiffInputFiles=~ s/\:/ /g;
		
		runCuffdiff($extraPathsRequired,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$cuffdiffOutDir,$GTF,$cuffdiffLibraryType,$cuffdiffNThreads,
		$cuffdiffFragBiasCorrect,$cuffdiffMultiReadCorrect,$cuffdiffMinAlignmentCount,$cuffdiffSeed,
		$cuffdiffInputFiles,$referenceSequence,$cuffdiffFDR,$cuffdiffLibraryNormalizationMethod,$cuffdiffComparisonLabels,$FPKMthreshold,$cuffdiffMaxBundleFrags,
		$cuffdiff_noEffectiveLengthCorrection,$cuffdiff_noLengthCorrection,$cuffdiff_dispersionMethod,$executionCreatedTempDir);		
	}
	if($mode=~ "3"){
		# replaces the initial separation between the samples from the different conditions (represented by ':'), by a space (required in this way by cuffnorm)
		$cuffnormInputFiles=~ s/\:/ /g;
		$cuffnormColors=~ s/,$//;
		$cuffnormColors=~ s/,/","/g;
		$cuffnormColors="\"".$cuffnormColors."\"";

		runCuffnorm($extraPathsRequired,$cuffnormNThreads,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$GTF,$cuffnormOutputFormat,
		$cuffnormLibraryNormalizationMethod,$cuffnormSeed,$cuffnormNormalization,$cuffnormLibraryType,$cuffnormOutDir,$cuffnormLabels,$cuffnormInputFiles);

		calculateCorrelationsAndPCA($extraPathsRequired,$cuffnormOutDir,$cuffnormInputFiles,$cuffnormColors,$executionCreatedTempDir);
	}	
		
}

sub runCuffquant($$$$$$$$$$$$$$$$$){	
	my ($extraPathsRequired,$cufflinksPath,$samtoolsPath,$alignmentsDir,$referenceSequence,$samples,$GTF,$libraryType,$nThreads,
		$cuffquantNThreads,$cuffquantSeed,$cuffquantFragBiasCorrect,$cuffquantMultiReadCorrect,$cuffquantOutDir,$cuffquantMaxBundleFrags,
		$cuffquant_noEffectiveLengthCorrection,$cuffquant_noLengthCorrection)=@_;

	my @files=split(',',$samples);
	my @libraries=split(',',$libraryType);
	#checks how many files there are to process
	#in case of paired end files, they are counted as being only one file
	my $counter=@files;		
	my $memory=$counter;
	
	#in case that the number of allowed threads exceeds the number of samples to process,
	#the number of allowed threads is adjusted to the number of samples
	if($nThreads>$memory){$nThreads=$memory};
	
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

		#sample library type (it can be different for each sample)
		#posibilities: unstranded, firststrand, secondstrand
		my $library=$libraries[$fileNumber];		
		
		#prepare files
		my $bamFile=$alignmentsDir.$inputFile."/accepted_hits.bam";
		if(!-e $bamFile){
				print STDERR "\n[ERROR cuffquant]: $bamFile file does not exist\n\n";
				exit(-1);
		}				
				
		$setOfThreads[$currentThread]=threads->create(\&cuffquant_cuffdiff_cuffnorm::runCuffquant,
		$extraPathsRequired,$cufflinksPath,$samtoolsPath,$alignmentsDir,$referenceSequence,$samples,$GTF,$libraryType,$nThreads,
		$cuffquantNThreads,$cuffquantSeed,$cuffquantFragBiasCorrect,$cuffquantMultiReadCorrect,$cuffquantOutDir,$cuffquantMaxBundleFrags,
		$cuffquant_noEffectiveLengthCorrection,$cuffquant_noLengthCorrection,$inputFile);		
		
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

sub runCuffdiff($$$$$$$$$$$$$$$$$$$$$$$){
	my($extraPathsRequired,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$cuffdiffOutDir,$GTF,$cuffdiffLibraryType,$cuffdiffNThreads,
		$cuffdiffFragBiasCorrect,$cuffdiffMultiReadCorrect,$cuffdiffMinAlignmentCount,$cuffdiffSeed,
		$cuffdiffInputFiles,$referenceSequence,$cuffdiffFDR,$cuffdiffLibraryNormalizationMethod,$cuffdiffComparisonLabels,$FPKMthreshold,$cuffdiffMaxBundleFrags,
		$cuffdiff_noEffectiveLengthCorrection,$cuffdiff_noLengthCorrection,$cuffdiff_dispersionMethod,$executionCreatedTempDir)=@_;
	
	&cuffquant_cuffdiff_cuffnorm::runCuffdiff($extraPathsRequired,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$cuffdiffOutDir,$GTF,$cuffdiffLibraryType,$cuffdiffNThreads,
		$cuffdiffFragBiasCorrect,$cuffdiffMultiReadCorrect,$cuffdiffMinAlignmentCount,$cuffdiffSeed,
		$cuffdiffInputFiles,$referenceSequence,$cuffdiffFDR,$cuffdiffLibraryNormalizationMethod,$cuffdiffComparisonLabels,$cuffdiffMaxBundleFrags,
		$cuffdiff_noEffectiveLengthCorrection,$cuffdiff_noLengthCorrection,$cuffdiff_dispersionMethod);
		
	&cuffquant_cuffdiff_cuffnorm::createExcel($cuffdiffOutDir,$cuffdiffFDR,$cuffdiffComparisonLabels,$FPKMthreshold,$executionCreatedTempDir);
	
	&cuffquant_cuffdiff_cuffnorm::createGSEArnkFile($cuffdiffOutDir);
			
}

sub runCuffnorm($$$$$$$$$$$$$$){
	my($extraPathsRequired,$cuffnormNThreads,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$GTF,$cuffnormOutputFormat,
		$cuffnormLibraryNormalizationMethod,$cuffnormSeed,$cuffnormNormalization,$cuffnormLibraryType,$cuffnormOutDir,$cuffnormLabels,$cuffnormInputFiles)=@_;

	&cuffquant_cuffdiff_cuffnorm::runCuffnorm($extraPathsRequired,$cuffnormNThreads,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$GTF,$cuffnormOutputFormat,
		$cuffnormLibraryNormalizationMethod,$cuffnormSeed,$cuffnormNormalization,$cuffnormLibraryType,$cuffnormOutDir,$cuffnormLabels,$cuffnormInputFiles);	
}

sub calculateCorrelationsAndPCA($$$$$){
	my ($extraPathsRequired,$cuffnormOutDir,$cuffnormInputFiles,$cuffnormColors,$executionCreatedTempDir)=@_;
	
	&cuffquant_cuffdiff_cuffnorm::calculateCorrelationsAndPCA_GeneLevel($extraPathsRequired,$cuffnormOutDir,$cuffnormInputFiles,$cuffnormColors,$executionCreatedTempDir);
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

