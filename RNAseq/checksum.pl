#!/usr/bin/perl -w

# Author: Osvaldo Grana
# Description: perform sample validation doing checksum
# v0.1		may2016

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use checksum;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;


sub main();
sub runChecksum($$$);
sub validateChecksumFile($$$);
sub checkThreads($$$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;

	my($samples,$nThreads,$outputFile,$fileWithChecksumCodesToValidate);

	undef($samples);
	undef($nThreads);
	undef($outputFile);
	undef($fileWithChecksumCodesToValidate);
	
	GetOptions(
		"nThreads=i"=>\$nThreads, # integer	
		"samples=s"=>\$samples,
		"outputFile=s"=>\$outputFile,
		"fileWithChecksumCodesToValidate=s"=>\$fileWithChecksumCodesToValidate
	);

#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 

	runChecksum($samples,$nThreads,$outputFile);
	validateChecksumFile($outputFile,$samples,$fileWithChecksumCodesToValidate);
}

sub runChecksum($$$){	
#creates sam files derived from bam files

	my ($samples,$nThreads,$outputFile)=@_;
	my @files=split(',',$samples);
	
	#checks how many files there are to process
	#in case of paired end files, they are counted as being only one file (only one alignment file is created anyway)
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
		if(!-e $inputFile){
				print STDERR "\n[ERROR checksum.pl]: $inputFile file does not exist\n\n";
				exit(-1);
		}				
				
		$setOfThreads[$currentThread]=threads->create(\&checksum::performChecksum,$inputFile,$outputFile);		
		
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

sub validateChecksumFile($$$){
# validation is done at two different levels:
# 1. checks if the SAME md5 code was given for TWO different files (basically reading all the md5 codes in the $outputFile, 
#    and comparing them in a pairwise fashion)
#	Why? although it is not frecuent, it could happen that two sample files (fastq/bam) had the same content due to an error (e.g., same file with two different file names)
#	This case wouldn't be detected by 'md5sum -c'

	my ($outputFile,$samples,$fileWithChecksumCodesToValidate)=@_;
	
	open(IN,$outputFile);
	my @filesToCheck=<IN>;
	close(IN);
	
	
	### FIRST CHECK
	for(my $i=0;$i<@filesToCheck-1; $i++){
		my $file1=$filesToCheck[$i];
		chomp($file1);
		
		for(my $j=$i+1;$j<@filesToCheck; $j++){
			my $file2=$filesToCheck[$j];
			chomp($file2);
			
			my $checksum1=(split(/  /,$file1))[0];
			my $fileName1=(split(/  /,$file1))[1];
			
			my $checksum2=(split(/  /,$file2))[0];
			my $fileName2=(split(/  /,$file2))[1];			
			
			if($checksum1 eq $checksum2){
				print STDERR "\n\n\n\n[ERROR]\n";
				print STDERR $fileName1."\n\tand\n".$fileName2."\n\thave the same checksum value!!!!\n";
				print STDERR "[Execution stopped]\n\n\n\n";
				exit(-1);						
			}
		}
	}
	
	### SECOND CHECK: validation of received md5 codes (with md5sum)
	if(-e $fileWithChecksumCodesToValidate){&checksum::validate($fileWithChecksumCodesToValidate,$samples)}
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
--checksum.pl--
			
	Execution error. Check out the execution modes in the following examples:
			
	Example for trimming a single end experiment
		perl preprocess.pl --nThreads 6 --mode 1 --seqtkPath /path --positions 0:3,0:3 --samples /path/file1.fq:/path/file1.trimmed.fq
				
	Example for trimming a paired end experiment
		perl preprocess.pl --nThreads 6 --mode 1 --seqtkPath /path --positions 0:3,0:3 --samples /path/file1.fq:/path/file1.trimmed.fq,/path/file2.fq:/path/file2.trimmed.fq


};
 
	print STDERR $usage;
	exit(1);
}

