#!/usr/bin/perl -w

# align.pl
# Author: Osvaldo Grana
# Description: performs aligning of samples
# v0.1		apr2014

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use align;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;





sub main();
sub buildIndex($$$);
sub align($$$$$$$$$$$$$$$$$$$$$$$$$$);
sub checkThreads($$$);
sub getAligningStatistics($$);
sub annotateAlignmentFiles($$$$$);
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

	my($mode,$tophatPath,$bowtiePath,$samtoolsPath,$bedtoolsPath,$peakAnnotatorPath,$referenceSequence);
	my($indexPrefixForReferenceSequence,$samples,$GTF);
	my($outputFileNames,$libraryType);
	my $nThreads=1;
	my $nTophatThreads=1;
	my $maxMultihits=5;
	my $readMismatches=2;
	my $segmentLength=20;
	my $segmentMismatches=1;
	my $spliceMismatches=0;
	my $reportSecondaryAlignments="false";
	my $bowtieVersion="--bowtie1";
	my $readEditDist=2;
	my $readGapLength=2;
	my $pairedEnd="false";
	my $solexaQualityEncoding="--solexa-quals";
	my $coverageSearch="";
	my $performFusionSearch="false";
	my $innerDist_StdDev="";
	my $alignmentsOutputDirectory="";
	my $useGTF="true";
	
	undef($mode);
	undef($tophatPath);
	undef($bowtiePath);
	undef($samtoolsPath);
	undef($referenceSequence);
	undef($samples);
	undef($GTF);
	undef($indexPrefixForReferenceSequence);
	undef($outputFileNames);
	undef($bedtoolsPath);
	undef($peakAnnotatorPath);
	
	GetOptions(
		"mode=s"=>\$mode, #string
		"tophatPath=s"=>\$tophatPath,	#string
		"bowtiePath=s"=>\$bowtiePath,	#string
		"samtoolsPath=s"=>\$samtoolsPath,	#string
		"bedtoolsPath=s"=>\$bedtoolsPath,	#string
		"peakAnnotatorPath=s"=>\$peakAnnotatorPath,	#string
		"nTophatThreads=i"=>\$nTophatThreads,	#integer
		"samples=s"=>\$samples,	#string
		"nThreads=i"=>\$nThreads,	#integer
		"maxMultihits=i"=>\$maxMultihits,	#integer
		"readMismatches=i"=>\$readMismatches,	#integer
		"segmentLength=i"=>\$segmentLength,	#integer
		"segmentMismatches=i"=>\$segmentMismatches,	#integer
		"spliceMismatches=i"=>\$spliceMismatches,	#integer
		"reportSecondaryAlignments=s"=>\$reportSecondaryAlignments,	#string
		"bowtieVersion=s"=>\$bowtieVersion,	#string
		"readEditDist=i"=>\$readEditDist,	#integer
		"readGapLength=i"=>\$readGapLength,	#integer
		"pairedEnd=s"=>\$pairedEnd,	#string
		"solexaQualityEncoding=s"=>\$solexaQualityEncoding,	#string
		"libraryType=s"=>\$libraryType,	#string
		"coverageSearch=s"=>\$coverageSearch,	#string
		"performFusionSearch=s"=>\$performFusionSearch,	#string
		"referenceSequence=s"=>\$referenceSequence,	#string
		"GTF=s"=>\$GTF, #string
		"indexPrefixForReferenceSequence=s"=>\$indexPrefixForReferenceSequence, #string
		"outputFileNames=s"=>\$outputFileNames, #string
		"innerDist_StdDev=s"=>\$innerDist_StdDev, #string
		"alignmentsOutputDirectory=s"=>\$alignmentsOutputDirectory,
		"useGTF=s"=>\$useGTF
	);
	
#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 
	
	if($mode=~ "1"){		
		buildIndex($bowtiePath,$referenceSequence,$indexPrefixForReferenceSequence);	       
	}
	if($mode=~ "2"){
		align($tophatPath,$bowtiePath,$samtoolsPath,$referenceSequence,$indexPrefixForReferenceSequence,$samples,$GTF,
		$nTophatThreads,$nThreads,$maxMultihits,$readMismatches,$segmentLength,$segmentMismatches,$spliceMismatches,$reportSecondaryAlignments,
		$bowtieVersion,$readEditDist,$readGapLength,$innerDist_StdDev,$pairedEnd,$solexaQualityEncoding,$libraryType,
		$coverageSearch,$performFusionSearch,$outputFileNames,$useGTF);		
	}
	if($mode=~ "3"){
		getAligningStatistics($outputFileNames,$alignmentsOutputDirectory);
	}
	if($mode=~ "4"){
		annotateAlignmentFiles($bedtoolsPath,$outputFileNames,$alignmentsOutputDirectory,$GTF,$peakAnnotatorPath);
	}	
		
}

sub buildIndex($$$){
	my ($bowtiePath,$referenceFasta,$indexPrefixForReferenceSequence)=@_;
	&align::buildIndex($bowtiePath,$referenceFasta,$indexPrefixForReferenceSequence);  	     
}

sub align($$$$$$$$$$$$$$$$$$$$$$$$$$){
	my ($tophatPath,$bowtiePath,$samtoolsPath,$referenceSequence,$indexPrefixForReferenceSequence,$samples,$GTF,
	$nTophatThreads,$nThreads,$maxMultihits,$readMismatches,$segmentLength,$segmentMismatches,$spliceMismatches,$reportSecondaryAlignments,
	$bowtieVersion,$readEditDist,$readGapLength,$innerDist_StdDev,$pairedEnd,$solexaQualityEncoding,$libraryType,
	$coverageSearch,$performFusionSearch,$outputFileNames,$useGTF)=@_;

	my @files=split(',',$samples);
	my @outputFiles=split(',',$outputFileNames);
	my @pairs=split(',',$innerDist_StdDev);# one pair for each sample
	my @libraries=split(',',$libraryType);# one value for each sample
	my @quals=split(',',$solexaQualityEncoding); # one value for each sample
	
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
		my $pos=$files[$fileNumber]; #read files
		my $pair=$pairs[$fileNumber];
		my $output=$outputFiles[$fileNumber]; #output file names
		my $qual=$quals[$fileNumber]; #samples qualities
		
		#sample library type (it can be different for each sample)
		#posibilities: unstranded, firststrand, secondstrand
		my $library=$libraries[$fileNumber]; 
		my($inputFile);			
		my $mateInnerDist="";
		my $mateStdDev="";
		
		#prepare files
		if($pos=~ /:/){ #paired-end experiment
			my $leftInputFile=(split(':',$pos))[0];
			my $rightInputFile=(split(':',$pos))[1];
			$mateInnerDist=(split(':',$pair))[0];
			$mateStdDev=(split(':',$pair))[1];
			
			#checks that both files exist
			if(!-e $leftInputFile){
				print STDERR "\n[ERROR]: $leftInputFile file does not exist\n\n";
				exit(-1);
			}
			if(!-e $rightInputFile){
				print STDERR "\n[ERROR]: $rightInputFile file does not exist\n\n";
				exit(-1);
			}			
			
			#paired-end, both input files together
			$inputFile=$leftInputFile." ".$rightInputFile;			
							
		}else{
			$inputFile=$pos;	
			
			if(!-e $inputFile){
				print STDERR "\n[ERROR]: $inputFile file does not exist\n\n";
				exit(-1);
			}				
		}		
		
		$setOfThreads[$currentThread]=threads->create(\&align::doAligning,$tophatPath,$bowtiePath,$samtoolsPath,
		$referenceSequence,$indexPrefixForReferenceSequence,$GTF,$nTophatThreads,$maxMultihits,
		$readMismatches,$segmentLength,$segmentMismatches,$spliceMismatches,$reportSecondaryAlignments,
		$bowtieVersion,$readEditDist,$readGapLength,$mateInnerDist,$mateStdDev,$pairedEnd,$qual,
		$library,$coverageSearch,$performFusionSearch,$inputFile,$output,$useGTF);		
		
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

sub getAligningStatistics($$){
#reads the alignment files generated by tophat and concatenate them in just a single file
#with alignment statistics
	
	my($outputFileNames,$alignmentsOutputDirectory)=@_;
	chomp($alignmentsOutputDirectory);
	
	my @outputFiles=split(',',$outputFileNames);
	
	my $statistics="";
	
	if($alignmentsOutputDirectory!~ /\/$/){
		$alignmentsOutputDirectory.="/";
	}
	
	foreach my $file(@outputFiles){
		chomp($file);
		
		$statistics.=&align::getTophatAligningStatistics($file);	
	}

	if(open(OUT,">",$alignmentsOutputDirectory."tophatAligningStatistics.xls")){
	
		print OUT $statistics;
		close(OUT);	
	}else{
		print STDERR "\n[ERROR]: Cannot create ".$alignmentsOutputDirectory."tophatAligningStatistics.xls\n\n";
	}
}

sub annotateAlignmentFiles($$$$$){
#converts bam output files to bed files and after annotations with PeakAnnotator

	my($bedtoolsPath,$outputFileNames,$alignmentsOutputDirectory,$GTF,$peakAnnotatorPath)=@_;
	chomp($alignmentsOutputDirectory);
	
	my @outputFiles=split(',',$outputFileNames);

	if($alignmentsOutputDirectory!~ /\/$/){
		$alignmentsOutputDirectory.="/";
	}
	
	my $annotationsDir=$alignmentsOutputDirectory."annotations/";
	File::Path::make_path($annotationsDir);
	
	if(-d $annotationsDir){	
		
		foreach my $file(@outputFiles){
			chomp($file);

			my @justTheFileName=split('/',$file);	

			&align::convertBAMtoBED($bedtoolsPath,$alignmentsOutputDirectory,$annotationsDir,$justTheFileName[@justTheFileName-1]);	
			&align::runPeakAnnotator($peakAnnotatorPath,$annotationsDir.$justTheFileName[@justTheFileName-1].".bed",$GTF,$annotationsDir,$justTheFileName[@justTheFileName-1]);
		} 	
	}else{
		print STDERR "\n[ERROR]: cannot create '$annotationsDir' directory\n\n";
		#exit(-1);
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
	exit(-1);
}

