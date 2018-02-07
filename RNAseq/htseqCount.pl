#!/usr/bin/perl -w

# Author: Osvaldo Grana
# Description: runs htseqcountPath to count reads in transcripts
# v0.2		nov2016

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use htseqCount;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;


sub main();
sub runSamtools($$$$$);
sub runHtseqCount($$$$$$$$$$$$$$$);
sub checkThreads($$$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;


	my($mode,$minaqual,$featuretype,$idattr,$htseqcountPath,$htseqCountPythonpath,$samtoolsPath,$samples,
	$GTF,$nThreads,$alignmentsDir,$htseqcountOutDir,$libraryType,$extraPathsRequired,$perl5lib,$executionCreatedTempDir);

	undef($mode);
	undef($minaqual);
	undef($featuretype);
	undef($idattr);
	undef($htseqcountPath);
	undef($htseqCountPythonpath);
	undef($samtoolsPath);
	undef($samples);
	undef($GTF);
	undef($nThreads);
	undef($alignmentsDir);
	undef($htseqcountOutDir);
	undef($nThreads);
	undef($libraryType);
	undef($extraPathsRequired);
	undef($perl5lib);
	undef($executionCreatedTempDir);

	GetOptions(
		"mode=s"=>\$mode,
		"minaqual=i"=>\$minaqual,	#integer
		"featuretype=s"=>\$featuretype, #string
		"idattr=s"=>\$idattr, #string
		"numberOfThreads=i"=>\$nThreads,	#integer
		"samtoolsPath=s"=>\$samtoolsPath,	#string
		"htseqcountPath=s"=>\$htseqcountPath,	#string
		"htseqCountPythonpath=s"=>\$htseqCountPythonpath,	#string
		"samples=s"=>\$samples,	#string
		"GTF=s"=>\$GTF,	#string
		"alignmentsDir=s"=>\$alignmentsDir,	#string
		"htseqcountOutDir=s"=>\$htseqcountOutDir,	#string
		"libraryType=s"=>\$libraryType,
		"extraPathsRequired=s"=>\$extraPathsRequired,
		"executionCreatedTempDir=s"=>\$executionCreatedTempDir,
		"perl5lib=s"=>\$perl5lib,
	);

#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 

	runSamtools($samtoolsPath,$samples,$alignmentsDir,$nThreads,$executionCreatedTempDir);				
	runHtseqCount($perl5lib,$extraPathsRequired,$mode,$minaqual,$featuretype,$idattr,$nThreads,$samtoolsPath,
		$htseqcountPath,$htseqCountPythonpath,$samples,$libraryType,$GTF,$alignmentsDir,$htseqcountOutDir);
	
	
		
}

sub runSamtools($$$$$){	
#to create sam files derived from bam files

	my ($samtoolsPath,$samples,$alignmentsDir,$nThreads,$executionCreatedTempDir)=@_;

	my @files=split(',',$samples);
	
	#checks how many files there are to process
	#in case of paired end files, they are counted as being only one file
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
		
		#prepare files
		my $bamFile=$alignmentsDir.$inputFile."/accepted_hits.bam";
		my $samFile=$alignmentsDir.$inputFile."/accepted_hits.sam";
		if(!-e $bamFile){
				print STDERR "\n[ERROR htseqcountPath.pl]: $bamFile file does not exist\n\n";
				exit(-1);
		}				
				
		$setOfThreads[$currentThread]=threads->create(\&htseqCount::runSamtools,
		$samtoolsPath,$bamFile,$samFile,$executionCreatedTempDir);		
		
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

sub runHtseqCount($$$$$$$$$$$$$$$){
	my($perl5lib,$extraPathsRequired,$mode,$minaqual,$featuretype,$idattr,$nThreads,$samtoolsPath,
	$htseqcountPath,$htseqCountPythonpath,$samples,$libraryType,$GTF,$alignmentsDir,$htseqcountOutDir)=@_;

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
		my $samFile=$alignmentsDir.$inputFile."/accepted_hits.sam";
		my $htseqcountPathOutputFile=$htseqcountOutDir.$inputFile.".xls";
		if(!-e $samFile){
				print STDERR "\n[ERROR htseqcountPath.pl]: $samFile file does not exist\n\n";
				exit(-1);
		}				
				
		$setOfThreads[$currentThread]=threads->create(\&htseqCount::runHtseqCount,
		$perl5lib,$extraPathsRequired,$mode,$minaqual,$featuretype,$idattr,$htseqcountPath,$htseqCountPythonpath,$GTF,
		$library,$samFile,$htseqcountPathOutputFile);

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

