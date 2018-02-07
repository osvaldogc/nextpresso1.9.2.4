#!/usr/bin/perl -w


# Author: Osvaldo Grana
# Description: performs transcripts quantification and assembly
# v0.3		oct2016

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use cufflinks;
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options
use threads;





sub main();
sub runCufflinks($$$$$$$$$$$$$$$$$$$);
sub calculateCorrelationsAndPCA_geneLevel($$$$);
sub calculateCorrelationsAndPCA_IsoformLevel($$$$);
sub runCuffmerge($$$$$$$$$$);
sub doSpikeInCorrection($$);
sub checkThreads($$$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;

	my($experimentName,$mode,$cufflinksPath,$samtoolsPath,$outDir,$alignmentsDir,$referenceSequence,$indexPrefixForReferenceSequence);
	my($samples,$GTF,$libraryType,$libraryNormalizationMethod,$cuffmergeOutDir,$extraPathsRequired,$inputSpikesRef,$workspace);
	my($noEffectiveLengthCorrection,$noLengthCorrection,$normalization,$executionCreatedTempDir);
	my $nThreads=1;
	my $nCufflinksThreads=1;
	my $fragBiasCorrect="false";
	my $multiReadCorrect="false";
	my $useGTFwithCufflinks="true";
	my $maxBundleFrags=500000;
	
	undef($mode);
	undef($referenceSequence);
	undef($samples);
	undef($GTF);
	undef($outDir);
	undef($alignmentsDir);
	undef($libraryNormalizationMethod);
	undef($experimentName);
	undef($extraPathsRequired);
	undef($inputSpikesRef);
	undef($noEffectiveLengthCorrection);
	undef($noLengthCorrection);
	undef($normalization);
	undef($executionCreatedTempDir);
	
	GetOptions(
		"mode=s"=>\$mode, #string
		"cufflinksPath=s"=>\$cufflinksPath, #string
		"samtoolsPath=s"=>\$samtoolsPath, #string
		"nCufflinksThreads=i"=>\$nCufflinksThreads,	#integer
		"samples=s"=>\$samples,	#string
		"nThreads=i"=>\$nThreads,	#integer
		"libraryType=s"=>\$libraryType,	#string
		"GTF=s"=>\$GTF, #string
		"indexPrefixForReferenceSequence=s"=>\$indexPrefixForReferenceSequence, #string
		"referenceSequence=s"=>\$referenceSequence, #string
		"outDir=s"=>\$outDir, #string
		"alignmentsDir=s"=>\$alignmentsDir, #string
		"cuffmergeOutDir=s"=>\$cuffmergeOutDir, #string
		"fragBiasCorrect=s"=>\$fragBiasCorrect, #string
		"multiReadCorrect=s"=>\$multiReadCorrect, #string
		"libraryNormalizationMethod=s"=>\$libraryNormalizationMethod, #string
		"experimentName=s"=>\$experimentName,
		"extraPathsRequired=s"=>\$extraPathsRequired,
		"inputSpikesRef=s"=>\$inputSpikesRef,
		"useGTFwithCufflinks=s"=>\$useGTFwithCufflinks,
		"maxBundleFrags=i"=>\$maxBundleFrags,
		"workspace=s"=>\$workspace,		
		"noEffectiveLengthCorrection=s"=>\$noEffectiveLengthCorrection,
		"noLengthCorrection=s"=>\$noLengthCorrection,
		"normalization=s"=>\$normalization,
		"executionCreatedTempDir=s"=>\$executionCreatedTempDir
	);
	
#	if(!defined($experimentName) && !defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 
	
	if($mode=~ "1"){		
		runCufflinks($extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$alignmentsDir,$referenceSequence,$samples,$GTF,
		$libraryType,$nThreads,$nCufflinksThreads,$fragBiasCorrect,$multiReadCorrect,$libraryNormalizationMethod,$useGTFwithCufflinks,$maxBundleFrags,
		$noEffectiveLengthCorrection,$noLengthCorrection,$normalization);
		calculateCorrelationsAndPCA_GeneLevel($experimentName,$outDir,$samples,$executionCreatedTempDir);
		calculateCorrelationsAndPCA_IsoformLevel($experimentName,$outDir,$samples,$executionCreatedTempDir);
	}
	if($mode=~ "2"){
		runCuffmerge($extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$referenceSequence,$samples,$GTF,$nThreads,$nCufflinksThreads,$cuffmergeOutDir);
	}
	if($mode=~ "3"){
		#gene level
		doSpikeInCorrection($workspace."cufflinks/samplesFPKMs_geneLevel.xls",$inputSpikesRef);
		#isoform level
		doSpikeInCorrection($workspace."cufflinks/samplesFPKMs_isoformLevel.xls",$inputSpikesRef);		
	}	
		
}

sub calculateCorrelationsAndPCA_GeneLevel($$$$){
	my ($experimentName,$outDir,$samples,$executionCreatedTempDir)=@_;
	
	&cufflinks::calculateCorrelationsAndPCA_GeneLevel($experimentName,$outDir,$samples,$executionCreatedTempDir);
}

sub calculateCorrelationsAndPCA_IsoformLevel($$$$){
	my ($experimentName,$outDir,$samples,$executionCreatedTempDir)=@_;
	
	&cufflinks::calculateCorrelationsAndPCA_IsoformLevel($experimentName,$outDir,$samples,$executionCreatedTempDir);
}

sub runCufflinks($$$$$$$$$$$$$$$$$$$){
	my ($extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$alignmentsDir,$referenceSequence,$samples,$GTF,
		$libraryType,$nThreads,$nCufflinksThreads,$fragBiasCorrect,$multiReadCorrect,$libraryNormalizationMethod,$useGTFwithCufflinks,$maxBundleFrags,
		$noEffectiveLengthCorrection,$noLengthCorrection,$normalization)=@_;

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
				print STDERR "\n[ERROR]: $bamFile file does not exist\n\n";
				exit(-1);
		}				
		
		$setOfThreads[$currentThread]=threads->create(\&cufflinks::runCufflinks,
		$extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$alignmentsDir,$referenceSequence,$samples,$GTF,
		$libraryType,$nThreads,$nCufflinksThreads,$fragBiasCorrect,$multiReadCorrect,$libraryNormalizationMethod,$useGTFwithCufflinks,$maxBundleFrags,
		$noEffectiveLengthCorrection,$noLengthCorrection,$normalization,$inputFile);		
		
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

sub runCuffmerge($$$$$$$$$$){
	my($extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$referenceSequence,$samples,$GTF,$nThreads,$nCufflinksThreads,$cuffmergeOutDir)=@_;

	my @files=split(',',$samples);
	
	my $assemblies="";
	
	foreach my $file(@files){
		chomp($file);
		
		$assemblies.=$outDir.$file."/transcripts.gtf\n";	
	}
	
	#all the transcript file names are put together in a common file for cuffmerge
	my $assembliesFile=$outDir."assemblies.txt";
	open(OUT,">",$assembliesFile);
	print OUT $assemblies;
	close(OUT);
	
	&cufflinks::runCuffmerge($extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$GTF,$referenceSequence,$nCufflinksThreads,$assembliesFile,$cuffmergeOutDir);	
}

sub doSpikeInCorrection($$){
	my($inputFPKMmatrix,$inputSpikesRef)=@_;
					
	# name prefixes (of spike-in control mixes) need to be inferred from the spike-in control REF file
	#to afterwards find their positions inside the $inputFPKMmatrix (for the R script)

	# reads the file as a list 'aggregate => 0' to infer name prefixes
	open(SPIKESREF,$inputSpikesRef);
	my($prefix);
	undef($prefix);	
	foreach my $line(<SPIKESREF>){
		chomp($line);
		
		if($line=~ /\>/){
			$line=~ s/\>//;
			$prefix=(split('-',$line))[0]; 
			last;
		}
	}		
	close(SPIKESREF);
	
	#now it finds the first and last position of the prefix in the matrix
	open(MATRIX,$inputFPKMmatrix);
	my @matrix =<MATRIX>;
	close(MATRIX);
	my($firstPosition,$lastPosition);
	undef($firstPosition);
	undef($lastPosition);
	my $dimMatrix="";
	
	for(my $i=0;$i<@matrix;$i++){
		my $line=$matrix[$i];
		chomp($line);
		
		if(!defined($firstPosition) && $line=~ /^$prefix\-/){
			$firstPosition=$i+1; # R is 1-based
			
			my @aux=split('\t',$line);
			$dimMatrix=@aux; # and not @aux-1 because R is one based
		}
		
		if($line=~ /^$prefix\-/){
			$lastPosition=$i+1; # R is 1-based
		}		
	}	
	undef(@matrix);

	my $outputFPKMmatrix=$inputFPKMmatrix;
	$outputFPKMmatrix=~ s/xls$/correctedWithSpikeIncontrols\.xls/;
	&cufflinks::doSpikeInCorrection($inputFPKMmatrix,$outputFPKMmatrix,$firstPosition,$lastPosition,$dimMatrix);
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

