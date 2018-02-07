#!/usr/bin/perl -w

# preprocess.pl
# Author: Osvaldo Grana
# Description: performs trimming and downsampling of reads
# v0.1		apr2014

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use preprocess;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;





sub main();
sub doTrimming($$$$);
sub doFastqScreen($$$$);
sub fromBamToFastq($$$);
sub checkThreads($$$);
sub help();

main();




# function: main() 
# arguments:
# $mode can be 1 (fastqc), 2 (fastqScreen) or both
# $fastqcPath: path for fastqc
# $fastqcOutDir: output directory for fastqc
# $fastqScreenPath: path for fastqscreen
# $fastqScreenOutDir: output directory for fastqscreen
# $inputFile: samples list. Sample file names must be separated with ",". 
# In case of a paired-end experiment, sample file names (left and right) have to be separated with ":"

sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;

	my($mode,$seqtkPath,$positions,$samples,$conditions,$bedtoolsPath);
	my $nThreads=1;
	
	undef($mode);
	undef($seqtkPath);
	undef($positions);
	undef($samples);
	undef($conditions);
	undef($bedtoolsPath);
	
	GetOptions(
		"mode=s"=>\$mode, #string
		"seqtkPath=s"=>\$seqtkPath,	#string
		"positions=s"=>\$positions,	#string			
		"conditions=s"=>\$conditions,	#string
		"samples=s"=>\$samples,	#string
		"nThreads=i"=>\$nThreads,	#integer
		"bedtoolsPath=s"=>\$bedtoolsPath,	#string
	);

#	if(!defined($seqtkPath) || !defined($mode) || !defined($samples) || (!defined($conditions) && !defined($positions)))
#	{
#		help();
#	} 
	
	if($mode eq "1"){		
		doTrimming($seqtkPath,$positions,$nThreads,$samples);		
	}elsif($mode eq "2"){
		#fastqc only
		doDownsampling($seqtkPath,$conditions,$samples,$nThreads);
	}elsif($mode eq "3"){
		#converts bam files to fastq files
		fromBamToFastq($bedtoolsPath,$samples,$nThreads);
	}
	
		
}

sub doTrimming($$$$){
	my ($seqtkPath,$positions,$nThreads,$samples)=@_;

	my @files=split(',',$samples);
	my @cuts=split(',',$positions);	
	
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
		my $pos=$files[$fileNumber];
		my $cut=$cuts[$fileNumber];
		my($inputFile,$outputFile);
		undef($inputFile);
		undef($outputFile);
		my ($nNucleotidesLeftEnd,$nNucleotidesRightEnd);
		undef($nNucleotidesLeftEnd);
		undef($nNucleotidesRightEnd);				

		#prepare files
		if($pos=~ /:/){
			$inputFile=(split(':',$pos))[0];
			$outputFile=(split(':',$pos))[1];
			$nNucleotidesLeftEnd=(split(':',$cut))[0];
			$nNucleotidesRightEnd=(split(':',$cut))[1];

			#checks that both files exist
			if(!-e $inputFile){
				print STDERR "\n[ERROR preprocess.pl->doTrimming]: $inputFile file does not exist\n\n";
				exit(-1);
			}				
		}

		
		$setOfThreads[$currentThread]=threads->create(\&preprocess::trimming,$seqtkPath,$nNucleotidesLeftEnd,$nNucleotidesRightEnd,$inputFile,$outputFile);
		$currentThread++;
		$fileNumber++;
		$counter++;
		
		#checks if the number of active threads reached the maximum of threads allowed
		if($counter>=$nThreads){
			$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
		}					

#		#right sample
#		if(defined($outputFile)){
#			$setOfThreads[$currentThread]=threads->create(\&preprocess::trimming,$seqtkPath,$nNucleotidesLeftEnd,$nNucleotidesRightEnd,$inputFile,$outputFile);			
#			$currentThread++;
#			
#			#****IMPORTANT: in the case of a paired end experiment $fileNumber must not be increased
#			#here, because the two files are originally counted as only one (the split is done just after
#			#within the while loop)
#
#
#			#checks if the number of actived threads reached the maximum of threads allowed
#			$counter++;
#
#			#checks if the number of active threads reached the maximum of threads allowed
#			if($counter>=$nThreads){
#				$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
#			}					
#		}			

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

sub doDownsampling($$$$){
	my ($seqtkPath,$conditions,$samples,$nThreads)=@_;

	my @files=split(',',$samples);
	my @cuts=split(',',$conditions);
	
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
		my $pos=$files[$fileNumber];
		my $cut=$cuts[$fileNumber];
		my($inputFile,$outputFile);
		undef($inputFile);
		undef($outputFile);
		my ($seed,$nReads);
		undef($seed);
		undef($nReads);
					

		#prepare files
		if($pos=~ /:/){
			$inputFile=(split(':',$pos))[0];
			$outputFile=(split(':',$pos))[1];
			$seed=(split(':',$cut))[0];
			$nReads=(split(':',$cut))[1];

			#checks that both files exist
			if(!-e $inputFile){
				print STDERR "\n[ERROR preprocess.pl->doDownsampling]: $inputFile file does not exist\n\n";
				exit(-1);
			}				
		}

		
		$setOfThreads[$currentThread]=threads->create(\&preprocess::downsampling,$seqtkPath,$seed,$nReads,$inputFile,$outputFile);
		$currentThread++;
		$fileNumber++;
		$counter++;
		
		#checks if the number of active threads reached the maximum of threads allowed
		if($counter>=$nThreads){
			$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
		}					

#		#right sample
#		if(defined($outputFile)){
#			$setOfThreads[$currentThread]=threads->create(\&preprocess::downsampling,$seqtkPath,$seed,$nReads,$inputFile,$outputFile);
#			$currentThread++;
#			
#			#****IMPORTANT: in the case of a paired end experiment $fileNumber must not be increased
#			#here, because the two files are originally counted as only one (the split is done just after
#			#within the while loop)
#
#
#			#checks if the number of actived threads reached the maximum of threads allowed
#			$counter++;
#		
#			#checks if the number of active threads reached the maximum of threads allowed
#			if($counter>=$nThreads){
#				$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
#			}					
#		}			

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

sub fromBamToFastq($$$){
	my($bedtoolsPath,$samples,$nThreads)=@_;

	my @files=split(',',$samples);
		
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
		my $pos=$files[$fileNumber];
		my($inputFile,$outputFile);
		undef($inputFile);
		undef($outputFile);				

		#prepare files
		if($pos=~ /:/){
			$inputFile=(split(':',$pos))[0];
			$outputFile=(split(':',$pos))[1];

			#checks that both files exist
			if(!-e $inputFile){
				print STDERR "\n[ERROR preprocess.pl->fromBamToFastq]: $inputFile file does not exist\n\n";
				exit(-1);
			}				
		}

		
		$setOfThreads[$currentThread]=threads->create(\&preprocess::fromBamToFastq,$bedtoolsPath,$inputFile,$outputFile);
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
--preprocess.pl--
			
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

