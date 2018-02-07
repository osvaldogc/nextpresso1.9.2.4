#!/usr/bin/perl -w

# quality.pl
# Author: Osvaldo Grana
# Description: quality checks for the samples with fastqc and/or fastqscreen
# v0.2	ene2017 - adds PERL5LIB as argument for fastqscreen
# v0.1	mar2014

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use quality;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;



sub main();
sub doFastqc($$$$);
sub doFastqScreen($$$$$$$$);
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
# $samples: samples list. Sample file names must be separated with ",". 
# In case of a paired-end experiment, sample file names (left and right) have to be separated with ":"

sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;

	my($mode,$fastqcPath,$fastqcOutDir,$fastqScreenPath,$fastqScreenOutDir,$samples,$fastqScreenConf,$fastqFileIlluminaQualityEncodingForFastqScreen,$subset,$perl5lib);
	my $nThreads=1;
	$fastqFileIlluminaQualityEncodingForFastqScreen=""; #no explicit specification for illumina quality encoding
	undef($mode);
	undef($fastqcPath);
	undef($fastqcOutDir);
	undef($fastqScreenPath);
	undef($fastqScreenOutDir);
	undef($samples);
	undef($fastqScreenConf);
	undef($subset);
	undef($perl5lib);
	
	GetOptions(
		"mode=s"=>\$mode, #string
		"fastqcPath=s"=>\$fastqcPath,	#string
		"fastqcOutDir=s"=>\$fastqcOutDir,	#string				
		"fastqScreenPath=s"=>\$fastqScreenPath,		#string
		"fastqScreenOutDir=s"=>\$fastqScreenOutDir,	#string
		"samples=s"=>\$samples,	#string
		"nThreads=i"=>\$nThreads,	#integer
		"fastqScreenConfFile=s"=>\$fastqScreenConf, #string
		"subset=i"=>\$subset, #integer
		"fastqFileQualityEncodingForFastqScreen=s"=>\$fastqFileIlluminaQualityEncodingForFastqScreen,
		"perl5lib=s"=>\$perl5lib
	);

#	if(!defined($mode) || !defined($samples) || ((!defined($fastqcPath) || !defined($fastqcOutDir)) && (!defined($fastqScreenPath) || !defined($fastqScreenOutDir) || !defined($fastqScreenConf) || !defined($subset))))
#	{
#		help();
#	} 
	
	if($mode eq "1"){
		#fastqc only
		doFastqc($fastqcPath,$fastqcOutDir,$samples,$nThreads);		
	}elsif($mode eq "2"){
		#fastqc only
		doFastqScreen($perl5lib,$fastqScreenPath,$fastqScreenOutDir,$samples,$nThreads,$fastqScreenConf,$subset,$fastqFileIlluminaQualityEncodingForFastqScreen);		
	}elsif($mode eq "both"){
		#do both
		doFastqc($fastqcPath,$fastqcOutDir,$samples,$nThreads);	
		doFastqScreen($perl5lib,$fastqScreenPath,$fastqScreenOutDir,$samples,$nThreads,$fastqScreenConf,$subset,$fastqFileIlluminaQualityEncodingForFastqScreen);		
	}
	
		
}

sub doFastqc($$$$){
	my ($fastqcPath,$fastqcOutDir,$samples,$nThreads)=@_;
	
	#looks for fastqc
	if($fastqcPath!~ /\/$/){$fastqcPath.="/"}
	if(!-e $fastqcPath."fastqc"){
			print STDERR "\n[ERROR]: The program ".$fastqcPath."fastqc doesn't exist\n\n";
			exit(-1);
	}
	
	#creates fastqc output directory
	if(!-d $fastqcOutDir){
		File::Path::make_path($fastqcOutDir);
	}
	
	#finishes execution if the output directory was not created
	if(-d $fastqcOutDir){
		my @files=split(',',$samples);
		
		#checks how many files there are to process.
		#in case of paired end files, they are counted twice because they are processed independently
		my $counter=0;	
		foreach my $pos(@files){
			if($pos=~ /:/){
				$counter+=2;
			}else{
				$counter++;
			}
		}
		
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
			my($left,$right);
			undef($left);
			undef($right);			
			
			#prepare files
			if($pos=~ /:/){
				$left=(split(':',$pos))[0];
				$right=(split(':',$pos))[1];
				
				#checks that both files exist
				if(!-e $left){
					print STDERR "\n[ERROR]: $left file does not exist\n\n";
					exit(-1);
				}
				if(!-e $right){
					print STDERR "\n[ERROR]: $right file does not exist\n\n";
					exit(-1);
				}				
			}else{
				$left=$pos;
				
				#checks that the file exists
				if(!-e $left){
					print STDERR "\n[ERROR]: $left file does not exist\n\n";
					exit(-1);
				}
			}

			#left sample
			$setOfThreads[$currentThread]=threads->create(\&quality::singleFastqc,$fastqcPath,$left,$fastqcOutDir);
			$currentThread++;
			$fileNumber++;
			$counter++;

			#checks if the number of active threads reached the maximum of threads allowed
			if($counter>=$nThreads){
				$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
			}					
			
			#right sample
			if(defined($right)){
				$setOfThreads[$currentThread]=threads->create(\&quality::singleFastqc,$fastqcPath,$right,$fastqcOutDir);
				$currentThread++;								
					
				#****IMPORTANT: in the case of a paired end experiment $fileNumber must not be increased
				#here, because the two files are originally counted as only one (the split is done just after
				#within the while loop)
				
				
				#checks if the number of actived threads reached the maximum of threads allowed
				$counter++;
		
				#checks if the number of active threads reached the maximum of threads allowed
				if($counter>=$nThreads){
					$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
				}					
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
			if(defined $setOfThreads[$counter] && $setOfThreads[$counter]->is_joinable()){
				$setOfThreads[$counter]->join();
				#print "joined thread $counter\n";
				$joinedThreads++;
			}
			
			$counter++;
			
			if($counter>=@setOfThreads){$counter=0}
			if($joinedThreads>=@setOfThreads){last}			
		}
	
	}else{
		print STDERR "\n[ERROR]: $fastqcOutDir directory does not exist\n\n";
		exit(-1);
	}
	
	
}

sub doFastqScreen($$$$$$$$){
	my ($perl5lib,$fastqScreenPath,$fastqScreenOutDir,$samples,$nThreads,$fastqScreenConf,$subset,$fastqFileIlluminaQualityEncodingForFastqScreen)=@_;
	
	#looks for fastqScreen
	if($fastqScreenPath!~ /\/$/){$fastqScreenPath.="/"}
	if(!-e $fastqScreenPath."fastq_screen"){
			print STDERR "\n[ERROR]: The program ".$fastqScreenPath."fastq_screen doesn't exist\n\n";
			exit(-1);
	}
	
	#deletes and creates fastqScreen output directory
	#File::Path::remove_tree($fastqScreenOutDir);
	File::Path::make_path($fastqScreenOutDir);
	
	
	#finishes execution if the output directory was not created
	if(-d $fastqScreenOutDir){
		my @files=split(',',$samples);
		
		#checks how many files there are to process
		# (paired end files go together in just one step, so their counted only once)
		my $counter=@files;
				
		#creates an array of threads
		my @setOfThreads=();
		my @hasThisThreadFinishedExecution=();
		for(my $i=0;$i<$counter;$i++){
			undef $setOfThreads[$i];
			$hasThisThreadFinishedExecution[$i]=0;
		}
		
		$counter=0; #resets the counter to account for the position to include the current thread in the array
		#in case that the number of allowed threads exceeds the number of samples to process,
		#the number of allowed threads is adjusted to the number of samples
		if($nThreads>@files){$nThreads=@files};
		
	
		for(my $i=0;$i<@files;$i++){
			my $pos=$files[$i];
			chomp($pos);
			
			my($left,$right);
			undef($left);
			undef($right);			
			
			#prepare files
			if($pos=~ /:/){
				$left=(split(':',$pos))[0];
				$right=(split(':',$pos))[1];
				
				#checks that both files exist
				if(!-e $left){
					print STDERR "\n[ERROR]: $left file does not exist\n\n";
					exit(-1);
				}
				if(!-e $right){
					print STDERR "\n[ERROR]: $right file does not exist\n\n";
					exit(-1);
				}				
			}else{
				$left=$pos;
				
				#checks that the file exists
				if(!-e $left){
					print STDERR "\n[ERROR]: $left file does not exist\n\n";
					exit(-1);
				}
			}

			#paired end ($left should be defined in any case)
			my $inputFile="";
			#if paired end experiment
			if(defined($right)){
				$inputFile=$left." ".$right;
			}else{ #if single end experiment, only left file
				$inputFile=$left;			
			}
			
			$setOfThreads[$i]=threads->create(\&quality::fastqScreen,$perl5lib,$fastqScreenPath,$inputFile,$fastqScreenOutDir,$fastqScreenConf,$subset,$fastqFileIlluminaQualityEncodingForFastqScreen);			
			
			$counter++;

			#checks if the number of active threads reached the maximum of threads allowed
			if($counter>=$nThreads){
				$counter=checkThreads(\@setOfThreads,$counter,\@hasThisThreadFinishedExecution);
			}					
					
		}#for(my $i=0;$i<@files;$i++)
			
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
	
	}else{
		print STDERR "\n[ERROR]: $fastqScreenOutDir directory does not exist\n\n";
		exit(-1);
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
			
	Example for fastqc in a single end experiment
		perl quality.pl --nThreads 6 --mode 1 --fastqcPath /path --fastqcOutDir /path --samples /path/file1.fq,/path/file2.fq
				
	Example for fastqc in a paired end experiment
		perl quality.pl --nThreads 6 --mode 1 --fastqcPath /path --fastqcOutDir /path --samples /path/file1_R1.fq:/path/file1_R2.fq,/path/file2_R1.fq:/path/file2_R2.fq
				
	Example for fastqScreen in a single end experiment
		perl quality.pl [--fastqFileQualityEncodingForFastqScreen illumina1_3] --nThreads 6 --mode 2 --fastqScreenConfFile /path/fastq_screen.conf --subset 1000000 --fastqScreenPath /path --fastqScreenOutDir /path --samples /path/file1.fq,/path/file2.fq
				
	Example for fastqScreen in a paired end experiment
		perl quality.pl [--fastqFileQualityEncodingForFastqScreen illumina1_3] --nThreads 6 --mode 2 --fastqScreenConfFile /path/fastq_screen.conf --subset 1000000 --fastqScreenPath /path --fastqScreenOutDir /path --samples /path/file1_R1.fq:/path/file1_R2.fq,/path/file2_R1.fq:/path/file2_R2.fq
				
	Example for fastqc and fastqScreen in a single end experiment
		perl quality.pl [--fastqFileQualityEncodingForFastqScreen illumina1_3] --nThreads 6 --mode both --fastqScreenConfFile /path/fastq_screen.conf --subset 1000000 --fastqcPath /path --fastqcOutDir /path --fastqScreenPath /path --fastqScreenOutDir /path --samples /path/file1.fq,/path/file2.fq

};
 
	print STDERR $usage;
	exit(1);
}

