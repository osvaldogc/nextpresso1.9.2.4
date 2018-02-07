#!/usr/bin/perl -w

# Author: Osvaldo Grana
# Description: creates bigWigFiles from bam files
# v0.1		aug2014

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use createbigwigFiles;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;


sub main();
sub runSamtools($$$$);
sub runGenomeCoverageBed($$$$$$);
sub runBedGraphToBigWig($$$$$);
sub createBigWigTrackHeadersForUCSCgenomeBrowser($$$$);
sub checkThreads($$$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;

	my($chromosomeSizesFile,$bigDataUrlPrefix,$samtoolsPath,$bedtoolsPath,$bedGraphToBigWigPath,$samples,$nThreads,$alignmentsDir,$bigwigFilesDir,$experimentName);

	undef($chromosomeSizesFile);
	undef($bigDataUrlPrefix);
	undef($samtoolsPath);
	undef($bedtoolsPath);
	undef($bedGraphToBigWigPath);
	undef($samples);
	undef($nThreads);
	undef($alignmentsDir);
	undef($bigwigFilesDir);
	undef($experimentName);

	GetOptions(
		"chromosomeSizesFile=s"=>\$chromosomeSizesFile,
		"bigDataUrlPrefix=s"=>\$bigDataUrlPrefix,
		"numberOfThreads=i"=>\$nThreads, # integer	
		"samtoolsPath=s"=>\$samtoolsPath,	
		"bedtoolsPath=s"=>\$bedtoolsPath,
		"bedGraphToBigWigPath=s"=>\$bedGraphToBigWigPath,
		"samples=s"=>\$samples,	
		"alignmentsDir=s"=>\$alignmentsDir,	
		"bigwigFilesDir=s"=>\$bigwigFilesDir,	
		"experimentName=s"=>\$experimentName
	);

#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 

	runSamtools($samtoolsPath,$samples,$alignmentsDir,$nThreads);				
	runGenomeCoverageBed($chromosomeSizesFile,$bedtoolsPath,$samples,$nThreads,$alignmentsDir,$bigwigFilesDir);
	runBedGraphToBigWig($chromosomeSizesFile,$bedGraphToBigWigPath,$samples,$nThreads,$bigwigFilesDir);
	createBigWigTrackHeadersForUCSCgenomeBrowser($experimentName,$bigDataUrlPrefix,$samples,$bigwigFilesDir);
}

sub runSamtools($$$$){	
#creates sam files derived from bam files

	my ($samtoolsPath,$samples,$alignmentsDir,$nThreads)=@_;

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
		my $bamFile=$alignmentsDir.$inputFile."/accepted_hits.bam";
		my $sortedBamFile=$alignmentsDir.$inputFile."/accepted_hits.sorted"; # **** only adds the PREFIX for the output file (the 'bam' extension is added by samtools)
		if(!-e $bamFile){
				print STDERR "\n[ERROR createbigwigFiles.pl]: $bamFile file does not exist\n\n";
				exit(-1);
		}				
				
		$setOfThreads[$currentThread]=threads->create(\&createbigwigFiles::runSamtools,
		$samtoolsPath,$bamFile,$sortedBamFile);		
		
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

sub runGenomeCoverageBed($$$$$$){
	my($chromosomeSizesFile,$bedtoolsPath,$samples,$nThreads,$alignmentsDir,$bigwigFilesDir)=@_;

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
	
	my @colors=("0,0,139","255,64,64","0,0,0","152,245,255","0,100,0","255,140,0","205,192,176");
	my $colorCounter=0;
	
	while(1){ #while there are files to be processed (threaded)
		my $inputFile=$files[$fileNumber]; #read files
		
		my $bedGraphFileHeader="'name=\"".$inputFile."\" description=\"".$inputFile."\" visibility=full color=".$colors[$colorCounter]."'";
		$colorCounter++;
		if($colorCounter>=@colors){
			$colorCounter=0;
		}		
		
		#prepare files
		my $sortedBamFile=$alignmentsDir.$inputFile."/accepted_hits.sorted.bam";
		my $bedGraphFile=$bigwigFilesDir.$inputFile.".bg";
		if(!-e $sortedBamFile){
				print STDERR "\n[ERROR createbigwigFiles.pl]: $sortedBamFile file does not exist\n\n";
				exit(-1);
		}				
				
		$setOfThreads[$currentThread]=threads->create(\&createbigwigFiles::runGenomeCoverageBed,
		$chromosomeSizesFile,$bedtoolsPath,$sortedBamFile,$bedGraphFile,$bedGraphFileHeader);

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

sub runBedGraphToBigWig($$$$$){
	my($chromosomeSizesFile,$bedGraphToBigWigPath,$samples,$nThreads,$bigwigFilesDir)=@_;

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
		my $inputFile=$files[$fileNumber]; #read files
		
		#prepare files
		my $bedGraphFile=$bigwigFilesDir.$inputFile.".bg";
		my $bigWigFile=$bigwigFilesDir.$inputFile.".bw";
		if(!-e $bedGraphFile){
				print STDERR "\n[ERROR createbigwigFiles.pl]: $bedGraphFile file does not exist\n\n";
				exit(-1);
		}				
				
		$setOfThreads[$currentThread]=threads->create(\&createbigwigFiles::runBedGraphToBigWig,
		$chromosomeSizesFile,$bedGraphToBigWigPath,$bedGraphFile,$bigWigFile);

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

sub createBigWigTrackHeadersForUCSCgenomeBrowser($$$$){
	my($experimentName,$bigDataUrlPrefix,$samples,$bigwigFilesDir)=@_;

	my @colors=("205,192,176","0,0,139","255,64,64","0,0,0","152,245,255","0,100,0","255,140,0");
	
	my @files=split(',',$samples);
	
	open(OUT,">",$bigwigFilesDir."trackHeaders.txt");
	
	if($bigDataUrlPrefix!~ /\/$/){
		$bigDataUrlPrefix.="/";
	}
	
	my $counter=0;
	foreach my $sample(@files){
		my $trackURL=$bigDataUrlPrefix.$sample.".bw";
		my $track="track type=bigWig name=\"".$experimentName."_".$sample."\" description=\"".$experimentName."_".$sample."\" color=".$colors[$counter]." bigDataUrl=".$trackURL."\n";
		
		print OUT $track;
		
		$counter++;
		if($counter>=@colors){
			$counter=0;
		}
	}
	
	close(OUT);
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

