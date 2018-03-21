#!/usr/bin/perl -W

# ExecutionLevels.pm 
# Author: Osvaldo Grana
# Description : executes the different levels
# v1.9		ene2017 - adds PCAs with all samples together, for both Cuffnorm and DESeq2
#
# v1.9.2	ene2018 - removes genes with expression levels below background + removes genes with flat pattern expression.
#			  After this removal a new GTF is created that contains only those genes that passed the filtering: cuffquant+cuffdiff+cuffnorm
#			  or htseqcount+deseq is run again using this reduced gene annotation (new GTF): added removeBackgroundLevelGenesANDFlatPatternGenes for
#			  cuffdiff and deseq2
# v1.9.2.1	mar2018 - the criteria to avoid filtering out a gene is now $percentageOfSamplesOverBackgroundLevel>=$required_PercentageOfSamplesOverBackgroundLevel
#			  for both cuffdiff and deseq2 (before the condition was set as strictly > , which filtered out genes where, for instance, its expression in three
#			  out of six samples were over background expression level. In that case, the gene should be kept and not filtered out).
#
#			  When there are no flat pattern genes to filter out, the number of flat pattern genes is now written as '0'. In v1.9.2 this value was the result of 
#			  doing 'wc -l ' of the number of lines of the file with the list of flat pattern filtered genes, but when no genes are filtered out due to this reason,
#			  the expected file does not exist, and the 'wc -l ' command results in a bug, now solved this way.




package ExecutionLevels;
use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use File::Copy;

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use queue;#loads queue module from Utils
use Miscellaneous;

if($Bin!~ /\/$/){$Bin.="/"}




########## level 0: converts raw read bam files to fastq files if needed or/and prepares reference and GTF index files in case of having spike-in control mixes
sub level_0{
	my($fileWithChecksumCodesToValidate,$workspace,$experimentName,$logfh,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$spikeInControlMixesParams,$referenceSequence,$GTF,
	$samples,$bedtoolsPath,$pairedEnd,$bowtiePath,$indexPrefixForReferenceSequence,
	$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;

	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait to perform this level execution

	
	# <----------------- only for testing ----------------------->	 
	#use Data::Dumper;
	#print STDERR "FILE:\n".Dumper($samples)."\n";
	
		
	# <----------------- checksum for sample validation ----------------------->
	my $ALLsamplesToCheck="";
	if($fileWithChecksumCodesToValidate eq ""){$fileWithChecksumCodesToValidate="NO_FILE"}
	foreach my $key (keys %$samples){
		$ALLsamplesToCheck.=$samples->{$key}{leftFile}.",";		
		
		if($pairedEnd eq "true"){
			$ALLsamplesToCheck.=$samples->{$key}{rightFile}[0].",";
		}
	}

	$ALLsamplesToCheck=~ s/,$//;
	my $checksumDir=$workspace."checksum/";
	File::Path::make_path($checksumDir);
	if(!-d $checksumDir){
				print STDERR "\n[ERROR]: $checksumDir directory does not exist\n";
				exit(-1);
	}
	
	my $checksumOutputFile=$checksumDir."checksum.md5.txt";
	#if checksum was already done in a previous execution, then this step is skipped
	#basically it checks for the existence of the $checksumOutputFile, if it exists this step is skipped
	if(!-e $checksumOutputFile){
		print STDOUT "\n[DOING] sample validation with checksum codes\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[DOING] checksum for sample validation\n";
		

		my $command="perl ".$Bin."checksum.pl --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep;
		$command.=" --samples ".$ALLsamplesToCheck." --outputFile ".$checksumOutputFile." --fileWithChecksumCodesToValidate ".$fileWithChecksumCodesToValidate." > ".$workspace."checksum.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";					
		
		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"checksum".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/waitingUntilChecksumHasFinished", $command,$wait,$multiCFlag);
		
		if($queueSystem ne "none"){
			copy($executionCreatedTempDir."/waitingUntilChecksumHasFinished",$executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/waitingUntilChecksumHasFinished to ".$executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting\n\n";		
			copy($executionCreatedTempDir."/waitingUntilChecksumHasFinished",$executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/waitingUntilChecksumHasFinished to ".$executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting\n\n";
			copy($executionCreatedTempDir."/waitingUntilChecksumHasFinished",$executionCreatedTempDir."/tophatIsWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/waitingUntilChecksumHasFinished to ".$executionCreatedTempDir."/tophatIsWaiting\n\n";	
		}
		
		open(FILE_WITH_CHECKSUM_SUMMARY,$workspace."checksum.log");
		my @checksumRESULT=<FILE_WITH_CHECKSUM_SUMMARY>;		
		close(FILE_WITH_CHECKSUM_SUMMARY);
		
		foreach my $sentence(@checksumRESULT){
			chomp($sentence);
			
			if($sentence=~ /computed checksums did NOT match/){
				print STDERR "[ERROR] ".$sentence."\n";
				print STDERR "\tplease check: ".$workspace."checksum.log\n";
				print STDERR "[Execution stopped]\n\n";
				exit(-1);
			}		
		}
		
	}else{
		print STDOUT "\n[CHECKSUM SKIPPED] This step was already done\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[CHECKSUM SKIPPED] This step was already done\n";	
	}
	
	#wait step until checksum is done
	if(-e $executionCreatedTempDir."/waitingUntilChecksumHasFinished"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/waitingUntilChecksumHasFinished",$queueSystem);	
	}
	
	# <----------------- bam file converstion to fastq file ----------------------->
	
	#checks if raw read files are in bam format:
	#if that was the case, then they are converted to fastq
	my $ALLsamples="";
	foreach my $key (keys %$samples){
		my $type=$samples->{$key}{type}[0];			

		if($type eq "bam"){

			#new directory for the converted fastq files				
			my $fastqFilesDir=$workspace."fastqfiles/";
			File::Path::make_path($fastqFilesDir);				

			if(!-d $fastqFilesDir){
				print STDERR "\n[ERROR]: $fastqFilesDir directory does not exist\n";
				exit(-1);
			}

			#input and output files
			my $bam=$samples->{$key}{leftFile};
			my @tokens=split('\/',$bam);
			my $fastq=$tokens[@tokens-1];
			$fastq=~ s/bam$/fastq/;
			$fastq=$fastqFilesDir.$fastq;

			#in case that the fastq file was created before (in a previous execution of the analysis),
			#then the conversion is not done again
			if(!-e $fastq){
				$ALLsamples.=$bam.":".$fastq.",";
				print STDOUT "\n[DOING] file conversion: ".$bam." to ".$fastq."\n";
				print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[DOING] conversion: ".$bam." to ".$fastq."\n";
			}
			else{
				print STDOUT "\n[SKIPPED] file conversion: ".$bam." to ".$fastq."\n";
				print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[SKIPPED] conversion: ".$bam." to ".$fastq."\n";
			}
			$samples->{$key}{leftFile}=$fastq;
			$samples->{$key}{type}[0]="fastq";

			#if paired-end experiment
			if($pairedEnd eq "true"){
				$bam=$samples->{$key}{rightFile}[0];
				@tokens=split('\/',$bam);
				$fastq=$tokens[@tokens-1];
				$fastq=~ s/bam$/fastq/;
				$fastq=$fastqFilesDir.$fastq;
				if(!-e $fastq){
					$ALLsamples.=$bam.":".$fastq.",";
					print STDOUT "\n[DOING] file conversion: ".$bam." to ".$fastq."\n";
					print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[DOING] conversion: ".$bam." to ".$fastq."\n";			
				}
				else{
					print STDOUT "\n[SKIPPED] file conversion: ".$bam." to ".$fastq."\n";
					print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[SKIPPED] conversion: ".$bam." to ".$fastq."\n";
				}
				$samples->{$key}{rightFile}[0]=$fastq;							
			}							
		}#if($type eq "bam")
	}
	
	$ALLsamples=~ s/,$//;
	
	#if there are bam files to convert to fastq
	if($ALLsamples ne ""){
		my $command="perl ".$Bin."preprocess.pl --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --mode 3 --bedtoolsPath ".$bedtoolsPath;
		$command.=" --samples ".$ALLsamples." > ".$workspace."bamToFastq.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";				

		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"b2Fq".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting", $command,$wait,$multiCFlag);
		
		if($queueSystem ne "none"){
			 copy($executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting",$executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting to ".$executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting\n\n";
			 copy($executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting",$executionCreatedTempDir."/tophatIsWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting to ".$executionCreatedTempDir."/tophatIsWaiting\n\n";	
		}
	}	
	
	# <----------------- reference indexing including spike-in control mixes referencing indexing and gtf adding ----------------------->
	
	# in case of having spike-in control mixes, the spike-in control reference have to be indexed together with the genomic reference,
	# and at the same time, the spike-in control annotation have to be added to the genomic annotation

	my $doSpikesAndGenomeRefIndexing=lc($spikeInControlMixesParams->[0]->{do});

	#decides to perform or not the reference indexing of both the genomic reference and the spike-in control mixes reference
	#and also the creation of a new combined annotation file (GTF)
	if($doSpikesAndGenomeRefIndexing eq "true"){
		print STDOUT "\n[DOING] spikes and genome references indexing\n";
		#new directory for the indexed reference
		my $refDir=$workspace."reference/";
		my $gtfDir=$workspace."gtf/";
		File::Path::make_path($refDir);				
		File::Path::make_path($gtfDir);

		if(!-d $refDir || !-d $gtfDir){
			if(!-d $refDir){print STDERR "\n[ERROR level 0]: $refDir directory does not exist\n"}
			if(!-d $gtfDir){print STDERR "\n[ERROR level 0]: $gtfDir directory does not exist\n"}
			exit(-1);
		}

		#concatenation of both references
		my $spikesRef=$spikeInControlMixesParams->[0]->{ref};
		my $bothRef=$refDir."reference.fa";
		my $bothRefFileOneBowtie=$refDir."reference.1.ebwt";
		my $bothRefPrefix=$refDir."reference";
		
		if(!-e $spikesRef || !-e $referenceSequence){
			if(!-e $spikesRef){
				print STDERR "\n[ERROR level 0]: You've selected to correct transcripts expression with spike-in control mixes, but there ";
				print "is no reference file for spike-in controls. Check the xml 'spikeInControlMixes' element definition.\n\n";
			}
			if(!-e $referenceSequence){print STDERR "\n[ERROR level 0]: $referenceSequence file does not exist\n"}
			exit(-1);
		}
		
		#if mixed reference file has not been created yet
		if(!-e $bothRef && !-e $bothRefFileOneBowtie){
			my $command="cat ".$referenceSequence." ".$spikesRef." > ".$bothRef;
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
			system($command);
			$referenceSequence=$bothRef;
			$indexPrefixForReferenceSequence=$bothRefPrefix;
			print $logfh (Miscellaneous->getCurrentDateAndTime())."new reference file: ".$referenceSequence."\n";
			print STDOUT "new reference file: ".$referenceSequence."\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())."new index file: ".$indexPrefixForReferenceSequence."\n";
			print STDOUT "new index file: ".$indexPrefixForReferenceSequence."\n";
			
			#indexing of both references
			print $logfh (Miscellaneous->getCurrentDateAndTime())."indexing both references\n";
			print STDOUT "indexing both references\n";
			$command="perl ".$Bin."align.pl --mode 1 --bowtiePath ".$bowtiePath." --indexPrefixForReferenceSequence ".$bothRefPrefix;
			$command.=" --referenceSequence ".$bothRef." >> ".$workspace."tophat.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"SpIdx".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">>".$executionCreatedTempDir."/tophatIsWaiting",$command,$wait,$multiCFlag);

		}else{
			print STDOUT "\n[DONE level 0]: $bothRef file already exists, skipping reference concatenation\n";
			print STDOUT "\n[DONE level 0]: $bothRef already exists and seems to be indexed for bowtie1, skipping bowtie1 indexing\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[DONE level 0]: $bothRef file already exists, skipping reference concatenation\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[DONE level 0]: $bothRef already exists and seems to be indexed for bowtie1, skipping bowtie1 indexing\n";			
			
		}
		
		#concatenation of both gtfs
		my $spikesGTF=$spikeInControlMixesParams->[0]->{gtf};		
		my $bothGTF=$gtfDir."annotation.gtf";
		#if gtf for bot annotations has not been created yet
		if(!-e $bothGTF){
			if(!-e $GTF || !-e $spikesGTF){
				if(!-e $GTF){print STDERR "\n[ERROR level 0]: $GTF file does not exist\n"}
				if(!-e $spikesGTF){print STDERR "\n[ERROR level 0]: You've selected to correct transcripts expression with spike-in control mixes, but there ";
				print "is no annotation GTF file for spike-in controls. Check the xml 'spikeInControlMixes' element definition.\n\n"}
				exit(-1);
			}


			my $command="cat ".$GTF." ".$spikesGTF." > ".$bothGTF;
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
			system($command);		
			$GTF=$bothGTF;
			print $logfh (Miscellaneous->getCurrentDateAndTime())."new GTF file: ".$GTF."\n";
			print STDOUT "new GTF file: ".$GTF."\n";
		}else{
			print STDOUT "\n[DONE level 0]: $bothGTF file already exists, skipping GTFs concatenation\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[DONE level 0]: $bothGTF file already exists, skipping GTFs concatenation\n";
		}

		return($bothRef,$bothGTF,$bothRefPrefix);
		
	}#if($performReferenceIndexing eq "true")
	else{
		#in order to avoid loosing the old files, they are also returned if no spike-in control mixes
		return($referenceSequence,$GTF,$indexPrefixForReferenceSequence);
	}
}

#level 1: analyses sequencing qualities and possible contamination sources
sub level_1{
    	my($perl5lib,$fastQCpath,$fastQScreenPath,$fastQScreenConf,$bowtiePath,$experimentName,$workspace,$referenceSequence,$GTF,
	$samples,$logfh,$executionCreatedTempDir,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$subset,$pairedEnd,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;


	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	if(-e $executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting",$queueSystem);	
	}
	
	# <----------------- only for testing ----------------------->	 
	#use Data::Dumper;
	#print STDERR "FILE:\n".Dumper($samples)."\n";
	

	# <----------------- FastQC ----------------------->
	my $fastQCoutputDir=$workspace."fastQC";
	File::Path::make_path($fastQCoutputDir);
	
	my $ALLsamples="";
	
	
	if(-d $fastQCoutputDir){		
		
		foreach my $key (keys %$samples){
			my $type=$samples->{$key}{type}[0];			

			if($type eq "fastq"){
				$ALLsamples.=$samples->{$key}{leftFile};					

				#if paired-end experiment it should also check the quality on the right end reads
				if($pairedEnd eq "true"){
					$ALLsamples.=":".$samples->{$key}{rightFile}[0];							
				}
				
				$ALLsamples.=",";			
			}
		}
		
		$ALLsamples=~ s/,$//;
		
		print STDOUT "\n[DOING] sequencing quality check\n";
		my $command="perl ".$Bin."quality.pl --mode 1 --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --fastqcPath ".$fastQCpath;
		$command.=" --fastqcOutDir ".$fastQCoutputDir." --samples ".$ALLsamples." >> ".$workspace."fastQC.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
			
		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"fQC".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);
		
	}else{
		print STDERR "\n[ERROR]: $fastQCoutputDir directory does not exist\n\n";
		exit(-1);
	}
	
	# <----------------- FastQScreen ----------------------->
	my $fastQScreenOutputDir=$workspace."fastQScreen";
	#the fastqscreen directory must be deleted in case it already exists,
	#otherwise it is not able to overwrite the new files
	File::Path::remove_tree($fastQScreenOutputDir);
	File::Path::make_path($fastQScreenOutputDir);
	
	if(-d $fastQScreenOutputDir){
		print STDOUT "\n[DOING] sample contamination check\n";
		my $command="";
		
		my $fastqFileQualityEncodingForFastqScreen=""; #="--fastqFileQualityEncodingForFastqScreen illumina1_3";
		if($fastqFileQualityEncodingForFastqScreen eq ""){
			$command="perl ".$Bin."quality.pl --mode 2 --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --fastqScreenConfFile ".$fastQScreenConf." --subset ".$subset." --perl5lib ".$perl5lib;
			$command.=" --fastqScreenPath ".$fastQScreenPath." --fastqScreenOutDir ".$fastQScreenOutputDir." --samples ".$ALLsamples." >> ".$workspace."fastQScreen.log 2>&1";	 
		}else{
			$command="perl ".$Bin."quality.pl ".$fastqFileQualityEncodingForFastqScreen." --mode 2 --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --conf ".$fastQScreenConf." --subset ".$subset." --perl5lib ".$perl5lib;
			$command.=" --fastqScreenPath ".$fastQScreenPath." --fastqScreenOutDir ".$fastQScreenOutputDir." --samples ".$ALLsamples." >> ".$workspace."fastQScreen.log 2>&1";
		}		
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
		
		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"fQS".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);												        	       
		
	}else{
		print STDERR "\n[ERROR]: $fastQScreenOutputDir directory does not exist\n\n";
		exit(-1);
	}	
}


#level 2: trimming && downsampling
sub level_2{
    	my($fastQCpath,$fastQScreenPath,$fastQScreenConf,$bowtiePath,$experimentName,$workspace,$referenceSequence,
	$GTF,$samples,$logfh,$executionCreatedTempDir,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$fastQScreenSubset,$pairedEnd,$seqtkPath,
	$seqtk_maximunNumberOfInstancesForDownSampling,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;

	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore);	
	my $wait=""; #empty=don't wait

	#first time for trimming
	if(-e $executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting",$queueSystem);	
	}

	# <----------------- only for testing ----------------------->	 
	#use Data::Dumper;
	#print #&ERR "FILE:\n".Dumper($samples)."\n";

	
	# <----------------- trimming ----------------------->	 
	my $trimmingDir=$workspace."trimmedSamples/";
	
	if(!-d $trimmingDir){
		File::Path::make_path($trimmingDir);
	}
	
	my $ALLsamples="";
	my $ALLnucleotidePositions="";
	
	my $doTrimming=0;
	
	#checks that the directory was built, then does the rest
	if(-d $trimmingDir){
			
		foreach my $key (keys %$samples){
	
			my $type=$samples->{$key}{type}[0];	
			#is this sample going to be trimmed?
			my $trimming=$samples->{$key}{trimming}[0]{do};

			if($type eq "fastq" && $trimming eq "true"){
				$doTrimming=1;
				$ALLnucleotidePositions.=$samples->{$key}{trimming}[0]{nNucleotidesLeftEnd}[0];
				$ALLnucleotidePositions.=":".$samples->{$key}{trimming}[0]{nNucleotidesRightEnd}[0].",";
				
				
				my $leftFastqFile=$samples->{$key}{leftFile};
				$ALLsamples.=$leftFastqFile;
				my @aux=split('\/',$leftFastqFile);
				my $trimmedLeftFastqFile=$aux[@aux-1]; #my $trimmedLeftFastqFile=$aux[length(@aux)-2];
				my $trimmedFile=$trimmingDir.$trimmedLeftFastqFile;
				$trimmedFile.=".trimmed";
				$ALLsamples.=":".$trimmedFile.",";

				#if paired-end experiment it should also check the quality on the right end reads
				if($pairedEnd eq "true"){
					$ALLnucleotidePositions.=$samples->{$key}{trimming}[0]{nNucleotidesLeftEnd}[0];
					$ALLnucleotidePositions.=":".$samples->{$key}{trimming}[0]{nNucleotidesRightEnd}[0].",";
				
					my $rightFastqFile=$samples->{$key}{rightFile}[0];
					$ALLsamples.=$rightFastqFile;
					@aux=split('\/',$rightFastqFile);
					my $trimmedRightFastqFile=$aux[@aux-1]; #my $trimmedRightFastqFile=$aux[length(@aux)-2];
					$trimmedFile=$trimmingDir.$trimmedRightFastqFile;
					$trimmedFile.=".trimmed";
					$ALLsamples.=":".$trimmedFile.",";
	
				}			
			
			}#if($type eq "fastq" && $trimming eq "true")

		}#foreach my $key (keys %$samples)
		
		$ALLsamples=~ s/,$//;
		$ALLnucleotidePositions=~ s/,$//;

		if($doTrimming){
			print STDOUT "\n[DOING] read trimming\n";
			my $command="perl ".$Bin."preprocess.pl --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --mode 1 --seqtkPath ".$seqtkPath;
			$command.=" --positions ".$ALLnucleotidePositions." --samples ".$ALLsamples." > ".$workspace."trimming.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"tri".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting",$command,$wait,$multiCFlag);
			
			#in case of doing trimming and downsampling:
			#1.downsampling has to wait for trimming, that's way trimmingANDdownsamplingAreWaiting it is read again below
			#2.tophat has to wait for trimming as well
			system("cat ".$executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting >> ".$executionCreatedTempDir."/tophatIsWaiting");
			system("cat ".$executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting >> ".$executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting");			
		}
	
	}else{
		print #&ERR "\n[ERROR]: $trimmingDir directory does not exist\n\n";
		exit(-1);
	}
	
	
		
	# <----------------- downsampling ----------------------->	 
	my $downsamplingDir=$workspace."downsampledSamples/";
	#attempts to build the directory
	
	#second time for downsampling
	if(-e $executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/trimmingANDdownsamplingAreWaiting",$queueSystem);	
	}

	if(!-d $downsamplingDir){
		File::Path::make_path($downsamplingDir);
	}
	
	$ALLsamples="";
	my $conditions="";
	
	my $doDownsampling=0;
	
	#checks that the directory was built, then do the rest
	if(-d $downsamplingDir){
			
		foreach my $key (keys %$samples){
		
			my $type=$samples->{$key}{type}[0];	
			#is this sample going to be downsampled?
			my $downsampling=$samples->{$key}{downsampling}[0]{do};

			if($type eq "fastq" && $downsampling eq "true"){
				$doDownsampling=1;
				
				$conditions.=$samples->{$key}{downsampling}[0]{seed}[0];
				$conditions.=":".$samples->{$key}{downsampling}[0]{nReads}[0].",";
				
				my $outDownsampledLeftFile=$downsamplingDir;
				my $outDownsampledRightFile=$downsamplingDir;
				my $leftFastqFile=$samples->{$key}{leftFile};
				my $rightFastqFile=$samples->{$key}{rightFile}[0];
				
				#if the sample was trimmed, it must use the trimmed file instead
				my $trimming=$samples->{$key}{trimming}[0]{do};				
				if($trimming eq	"true"){
					my @aux=split('\/',$leftFastqFile);
					my $trimmedLeftFastqFile=$aux[@aux-1]; #my $trimmedLeftFastqFile=$aux[length(@aux)-2];
					my $trimmedFile=$trimmingDir.$trimmedLeftFastqFile;
					$trimmedFile.=".trimmed";
					$leftFastqFile=$trimmedFile;
					$ALLsamples.=$leftFastqFile;
					
					#the output path is different for downsampling
					# that's why it is redefined here
					$outDownsampledLeftFile.=$trimmedLeftFastqFile.".trimmed.downsampled";
					$ALLsamples.=":".$outDownsampledLeftFile.",";
					
					if($pairedEnd eq "true"){
						#it's necessary to write the conditions once more if there is right file
						$conditions.=$samples->{$key}{downsampling}[0]{seed}[0];
						$conditions.=":".$samples->{$key}{downsampling}[0]{nReads}[0].",";
						
						@aux=split('\/',$rightFastqFile);
						my $trimmedRightFastqFile=$aux[@aux-1]; #my $trimmedRightFastqFile=$aux[length(@aux)-2];
						$trimmedFile=$trimmingDir.$trimmedRightFastqFile;
						$trimmedFile.=".trimmed";
						$rightFastqFile=$trimmedFile;	
						$ALLsamples.=$rightFastqFile;
						
						#the output path is different for downsampling
						# that's why it is redefined here
						$outDownsampledRightFile.=$trimmedRightFastqFile.".trimmed.downsampled";
						$ALLsamples.=":".$outDownsampledRightFile.",";				
				
					}
				}else{#the sample was not trimmed
					$ALLsamples.=$leftFastqFile;
					
					#the output path is different for downsampling
					# that's why it is redefined here
					my @aux=split('\/',$leftFastqFile);
					my $fileName=$aux[@aux-1]; #my $fileName=$aux[length(@aux)-2];
					$outDownsampledLeftFile.=$fileName.".downsampled";
					$ALLsamples.=":".$outDownsampledLeftFile.",";
					
					if($pairedEnd eq "true"){
						#it's necessary to write the conditions once more if there is right file
						$conditions.=$samples->{$key}{downsampling}[0]{seed}[0];
						$conditions.=":".$samples->{$key}{downsampling}[0]{nReads}[0].",";
						
						$ALLsamples.=$rightFastqFile;
						#the output path is different for downsampling
						# that's why it is redefined here
						@aux=split('\/',$rightFastqFile);
						$fileName=$aux[@aux-1]; #$fileName=$aux[length(@aux)-2];
						$outDownsampledRightFile.=$fileName.".downsampled";
						$ALLsamples.=":".$outDownsampledRightFile.",";											
					}
				}			
			
			}#if($type eq "fastq" && $trimming eq "true")

		}#foreach my $key (keys %$samples)
		
		$ALLsamples=~ s/,$//;
		$conditions=~ s/,$//;
		
		if($doDownsampling){
			print STDOUT "\n[DOING] read downsampling\n";
			my $command="perl ".$Bin."preprocess.pl --nThreads ".$seqtk_maximunNumberOfInstancesForDownSampling." --mode 2 --seqtkPath ".$seqtkPath;
			$command.=" --conditions ".$conditions." --samples ".$ALLsamples." > ".$workspace."downsampling.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

			#tophat has to wait for downsampling
			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"down".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">>".$executionCreatedTempDir."/tophatIsWaiting",$command,$wait,$multiCFlag);		
		
			system("cat ".$executionCreatedTempDir."/tophatIsWaiting >> ".$executionCreatedTempDir."/fastQCANDfastQScreenAreWaiting");
		}
		#parsing de ficheros $executionCreatedTempFile con funcion de queue.pm
		#el wait devuelto por esta funcion entra al nivel 3 como parametro $wait

	}else{
		print #&ERR "\n[ERROR]: $downsamplingDir directory does not exist\n\n";
		exit(-1);
	}

	#FINALLY:
	#checks wich one is going to be the fastq file for each sample that is going to be aligned by tophat,
	#depending on if it must be the original fastq file (meaning that nothing changed),
	#or if it was trimmed, downsampled or both trimmed and downsampled
	changeSampleNamesInCaseTheyWereTrimmedAndOrDownsampledBefore($samples,$workspace,$pairedEnd);
}

#level 3: reference indexing and aligning of samples
sub level_3{
    	my($tophatPath,$bowtiePath,$samtoolsPath,$bedtoolsPath,$peakAnnotatorPath,$referenceSequence,$indexPrefixForReferenceSequence,$samples,$GTF,$tophatParams,
	$nThreads,$workspace,$experimentName,$logfh,$executionCreatedTempDir,$pairedEnd,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;

	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	#first reading for the indexing step (in case of having to do indexing)
	if(-e $executionCreatedTempDir."/tophatIsWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/tophatIsWaiting",$queueSystem);	
	}
	
	# <----------------- only for testing ----------------------->	 
	#use Data::Dumper;
	#print STDERR "FILE:\n".Dumper($samples)."\n";
	

	# <----------------- INDEXING and ALIGNING ----------------------->
	my $alignmentsDir=$workspace."alignments/";
	File::Path::make_path($alignmentsDir);
	
	if(-d $alignmentsDir){	
				
		my $ALLsamples="";
		my $samplesOutputFileName="";
		my $libraryType="";
		my $innerDist_StdDev="";
		my $solexaQual="";		
		

		foreach my $key (keys %$samples){
			my $type=$samples->{$key}{type}[0];			
			if($type eq "fastq"){
				$ALLsamples.=$samples->{$key}{leftFile};					
				$samplesOutputFileName.=$alignmentsDir.${\$key}.",";	
				$libraryType.=$samples->{$key}{libraryType}[0].",";
				$innerDist_StdDev.=$samples->{$key}{mateInnerDist}[0].":".$samples->{$key}{mateStdDev}[0].",";
				
				#we have to control when $samples->{$key}{solexaQualityEncoding}[0] is empty
				if(!ref($samples->{$key}{solexaQualityEncoding}[0])){$solexaQual.=$samples->{$key}{solexaQualityEncoding}[0]}
				$solexaQual.=",";
				
				#if paired-end experiment it should also check the quality on the right end reads
				if(lc($pairedEnd) eq "true"){
					$ALLsamples.=":".$samples->{$key}{rightFile}[0];							
				}
				
				$ALLsamples.=",";			
			}
		}
		
		$ALLsamples=~ s/,$//;
		$samplesOutputFileName=~ s/,$//;
		$libraryType=~ s/,$//;
		$innerDist_StdDev=~ s/,$//;
		$solexaQual=~ s/,$//;		

		my $nTophatThreads=$tophatParams->[0]->{nTophatThreads};		
		my $maxMultihits=$tophatParams->[0]->{maxMultihits};		
		my $readMismatches=$tophatParams->[0]->{readMismatches};	
		my $segmentLength=$tophatParams->[0]->{segmentLength};		
		my $segmentMismatches=$tophatParams->[0]->{segmentMismatches};	
		my $spliceMismatches=$tophatParams->[0]->{spliceMismatches};		
		my $reportSecondaryAlignments=$tophatParams->[0]->{reportSecondaryAlignments};		
		my $bowtieVersion=$tophatParams->[0]->{bowtie};		
		my $readEditDist=$tophatParams->[0]->{readEditDist};		
		my $readGapLength=$tophatParams->[0]->{readGapLength};		
		my $referenceIndexing=$tophatParams->[0]->{referenceIndexing};		
		my $coverageSearch=$tophatParams->[0]->{coverageSearch}[0];
		my $performFusionSearch=$tophatParams->[0]->{fusionSearchExperiment}[0]->{performFusionSearch};
		my $useGTF=$tophatParams->[0]->{useGTF};

		my $command="";

		if(lc($referenceIndexing) eq "true"){
			print STDOUT "\n[DOING] reference indexing\n";
			$command="perl ".$Bin."align.pl --mode 1 --bowtiePath ".$bowtiePath." --indexPrefixForReferenceSequence ".$indexPrefixForReferenceSequence;
			$command.=" --referenceSequence ".$referenceSequence." >> ".$workspace."tophat.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
				
			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"Idx".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/tophatIsWaiting",$command,$wait,$multiCFlag);
		}

		#second reading, for the aligning step (it could have to wait for the indexing step also here - 2 lines above)
		if(-e $executionCreatedTempDir."/tophatIsWaiting"){
			$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/tophatIsWaiting",$queueSystem);	
		}

		print STDOUT "\n[DOING] read alignment\n";
		$command="perl ".$Bin."align.pl --mode 234 --nThreads ".$nThreads." --tophatPath ".$tophatPath." --bowtiePath ".$bowtiePath." --samtoolsPath ".$samtoolsPath;
		$command.=" --bedtoolsPath ".$bedtoolsPath." --peakAnnotatorPath ".$peakAnnotatorPath." --nTophatThreads ".$nTophatThreads." --bowtieVersion ".$bowtieVersion;
		$command.=" --solexaQualityEncoding ".$solexaQual." --libraryType ".$libraryType." --indexPrefixForReferenceSequence ".$indexPrefixForReferenceSequence;
		$command.=" --referenceSequence ".$referenceSequence." --alignmentsOutputDirectory ".$alignmentsDir." --useGTF ".$useGTF;
		if($GTF ne ""){
			if(!-e $GTF){
				print STDERR "\n[ERROR level 3]: $GTF file does not exist\n\n";
				exit(-1);
			}
			
			$command.=" --GTF ".$GTF;	
		}
		
		$command.=" --maxMultihits ".$maxMultihits." --readMismatches ".$readMismatches." --segmentLength ".$segmentLength." --segmentMismatches ".$segmentMismatches;
		$command.=" --spliceMismatches ".$spliceMismatches." --reportSecondaryAlignments ".$reportSecondaryAlignments." --readEditDist ".$readEditDist;
		$command.=" --readGapLength ".$readGapLength." --coverageSearch ".$coverageSearch." --performFusionSearch ".$performFusionSearch;
		$command.=" --pairedEnd ".$pairedEnd." --innerDist_StdDev ".$innerDist_StdDev." --samples ".$ALLsamples." --outputFileNames ".$samplesOutputFileName." >> ".$workspace."tophat.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
				
		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"aln".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">>".$executionCreatedTempDir."/cufflinksIsWaiting",$command,$wait,$multiCFlag);
		
		if($queueSystem ne "none"){
			copy($executionCreatedTempDir."/cufflinksIsWaiting",$executionCreatedTempDir."/cuffquantIsWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/cufflinksIsWaiting to ".$executionCreatedTempDir."/cuffquantIsWaiting\n";
			copy($executionCreatedTempDir."/cufflinksIsWaiting",$executionCreatedTempDir."/htseqCountIsWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/cufflinksIsWaiting to ".$executionCreatedTempDir."/htseqCountIsWaiting\n";
			copy($executionCreatedTempDir."/cufflinksIsWaiting",$executionCreatedTempDir."/bedGraphToBigWigIsWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/cufflinksIsWaiting to ".$executionCreatedTempDir."/bedGraphToBigWigIsWaiting\n";
			copy($executionCreatedTempDir."/cufflinksIsWaiting",$executionCreatedTempDir."/GSEAIsWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/cufflinksIsWaiting to ".$executionCreatedTempDir."/GSEAIsWaiting\n";
			copy($executionCreatedTempDir."/cufflinksIsWaiting",$executionCreatedTempDir."/tophatFusionIsWaiting") || print "\n[WARNING] Cannot copy ".$executionCreatedTempDir."/cufflinksIsWaiting to ".$executionCreatedTempDir."/tophatFusionIsWaiting\n";
		}		
	}else{
		print STDERR "\n[ERROR]: $alignmentsDir directory does not exist\n\n";
		exit(-1);
	}
}

#level 4: transcripts assembly and quantification (cufflinks and cuffmerge)
sub level_4{
    	my($extraPathsRequired,$spikeInControlMixesParams,$cufflinksPath,$samtoolsPath,$bedtoolsPath,$referenceSequence,$indexPrefixForReferenceSequence,$samples,$GTF,$cufflinksParams,
	$cuffmergeParams,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,$experimentName,$logfh,$executionCreatedTempDir,
	$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;

	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	if(-e $executionCreatedTempDir."/cufflinksIsWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/cufflinksIsWaiting",$queueSystem);	
	}

	
	# <----------------- only for testing ----------------------->	 
	#use Data::Dumper;
	#print STDERR "FILE:\n".Dumper($samples)."\n";
	

	# <----------------- Transcripts assembly and quantification ----------------------->
	my $cufflinksOutDir=$workspace."cufflinks/";
	my $alignmentsDir=$workspace."alignments/";
	my $cuffmergeOutDir=$workspace."cuffmerge/";
	File::Path::make_path($cufflinksOutDir);
	File::Path::make_path($cuffmergeOutDir);
	
	if(-d $cufflinksOutDir && -d $alignmentsDir && -d $cuffmergeOutDir){				
		my $ALLsamples="";
		my $libraryType="";		
		

		foreach my $key (keys %$samples){
			$ALLsamples.=${\$key}.",";
			$libraryType.=$samples->{$key}{libraryType}[0].",";
		}
		
		$ALLsamples=~ s/,$//; 
		$libraryType=~ s/,$//;		

		my $nCufflinksThreads=$cufflinksParams->[0]->{nThreads};		
		my $fragBiasCorrect=$cufflinksParams->[0]->{fragBiasCorrect};		
		my $multiReadCorrect=$cufflinksParams->[0]->{multiReadCorrect};	
		my $libraryNormalizationMethod=$cufflinksParams->[0]->{libraryNormalizationMethod};
		my $useGTFwithCufflinks=$cufflinksParams->[0]->{useGTF};
		my $maxBundleFrags=$cufflinksParams->[0]->{maxBundleFrags};
		my $noEffectiveLengthCorrection=$cufflinksParams->[0]->{noEffectiveLengthCorrection};
		my $noLengthCorrection=$cufflinksParams->[0]->{noLengthCorrection};
		my $normalization=$cufflinksParams->[0]->{normalization};
		
		print STDOUT "\n[DOING] transcript assembly\n";
		my $command="perl ".$Bin."cufflinks.pl --experimentName ".$experimentName;
		$command.=" --mode 12 --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --cufflinksPath ".$cufflinksPath." --nCufflinksThreads ".$nCufflinksThreads;
		$command.=" --libraryType ".$libraryType." --indexPrefixForReferenceSequence ".$indexPrefixForReferenceSequence." --referenceSequence ".$referenceSequence;
		$command.=" --outDir ".$cufflinksOutDir." --alignmentsDir ".$alignmentsDir." --cuffmergeOutDir ".$cuffmergeOutDir." --fragBiasCorrect ".$fragBiasCorrect." --multiReadCorrect ".$multiReadCorrect;
		$command.=" --libraryNormalizationMethod ".$libraryNormalizationMethod." --samtoolsPath ".$samtoolsPath." --extraPathsRequired ".$extraPathsRequired." --useGTF ".$useGTFwithCufflinks." --maxBundleFrags ".$maxBundleFrags;
		$command.=" --noEffectiveLengthCorrection ".$noEffectiveLengthCorrection." --noLengthCorrection ".$noLengthCorrection." --normalization ".$normalization;
		if($GTF ne ""){
			if(!-e $GTF){
				print STDERR "\n[ERROR level 4]: $GTF file does not exist\n\n";
				exit(-1);
			}
			
			$command.=" --GTF ".$GTF;
		}
		$command.=" --executionCreatedTempDir ".$executionCreatedTempDir;
		$command.=" --samples ".$ALLsamples." >> ".$workspace."cufflinks.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
				
		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"clink".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/spikesAreWaiting",$command,$wait,$multiCFlag);
		
		
		
		my $doSpikes=lc($spikeInControlMixesParams->[0]->{do});

		if(-e $executionCreatedTempDir."/spikesAreWaiting"){
			$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/spikesAreWaiting",$queueSystem);	
		}

		#Gene expression CORRECTION WITH SPIKE-in control mixes
		if($doSpikes eq "true"){
		
			my $inputSpikesRef=$spikeInControlMixesParams->[0]->{ref};			
			
			#executes cufflinks.pl for spike-in correction (--mode 3)
			my $command="perl ".$Bin."cufflinks.pl";
			$command.=" --mode 3 --workspace ".$workspace;
			$command.=" --inputSpikesRef ".$inputSpikesRef." >> ".$workspace."cufflinks.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"spkcor".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);

		}#if($doSpikes eq "true")
			
	}else{
		print STDERR "\n[ERROR]: Cufflinks cannot run:\n";
		if(!-d $cufflinksOutDir){print STDERR "\n[ERROR]: $cufflinksOutDir directory does not exist\n"}
		if(!-d $alignmentsDir){print STDERR "\n[ERROR]: $alignmentsDir directory does not exist\n"}
		if(!-d $cuffmergeOutDir){print STDERR "\n[ERROR]: $cuffmergeOutDir directory does not exist\n"}
		exit(-1);
	}
}

########## level 5: differential expression (cuffquant, cuffdiff and cuffnorm) ##########
sub level_5{
    	my($doSpikesAndGenomeRefIndexing,$initialGTF,$extraPathsRequired,$comparisons,$cufflinksPath,$samtoolsPath,$bedtoolsPath,$referenceSequence,$indexPrefixForReferenceSequence,
	$samples,$GTF,$cuffquantParams,$cuffnormParams,$cuffdiffParams,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,
	$experimentName,$logfh,$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;

	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	if(-e $executionCreatedTempDir."/cuffquantIsWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/cuffquantIsWaiting",$queueSystem);	
	}

	# <----------------- only for testing ----------------------->	 
#	use Data::Dumper;
#	print STDERR "FILE:\n".Dumper($comparisons)."\n";
	

	# <----------------- cuffquant, cuffdiff, cuffnorm ----------------------->
	my $cuffquantOutDir=$workspace."cuffquant/";
	my $cuffdiffOutDir=$workspace."cuffdiff/";
	my $cuffnormOutDir=$workspace."cuffnorm/";
	my $alignmentsDir=$workspace."alignments/";
	my $cuffmergeOutDir=$workspace."cuffmerge/";
	File::Path::make_path($cuffquantOutDir);
	File::Path::make_path($cuffdiffOutDir);
	File::Path::make_path($cuffnormOutDir);
	
	if(-d $cuffquantOutDir && -d $cuffdiffOutDir && -d $cuffnormOutDir && -d $alignmentsDir){				
		my $ALLsamples="";
		my $libraryType="";		
		
		foreach my $key (keys %$samples){
			#sample names
			$ALLsamples.=${\$key}.",";
			#library type
			$libraryType.=$samples->{$key}{libraryType}[0].",";
		}
		
		$ALLsamples=~ s/,$//;
		$libraryType=~ s/,$//;		
		
		my $cuffquantUseCuffmergeAssembly=$cuffquantParams->[0]->{useCuffmergeAssembly};
		my $cuffquantNThreads=$cuffquantParams->[0]->{nThreads};
		my $cuffquantFragBiasCorrect=$cuffquantParams->[0]->{fragBiasCorrect};
		my $cuffquantMultiReadCorrect=$cuffquantParams->[0]->{multiReadCorrect};
		my $cuffquantSeed=$cuffquantParams->[0]->{seed};
		my $cuffquant_noEffectiveLengthCorrection=$cuffquantParams->[0]->{noEffectiveLengthCorrection};
		my $cuffquant_noLengthCorrection=$cuffquantParams->[0]->{noLengthCorrection};
		
		my $cuffnormUseCuffmergeAssembly=$cuffnormParams->[0]->{useCuffmergeAssembly};
		my $cuffnormNThreads=$cuffnormParams->[0]->{nThreads};
		my $cuffnormOutputFormat=$cuffnormParams->[0]->{outputFormat};
		my $cuffnormlibraryNormalizationMethod=$cuffnormParams->[0]->{libraryNormalizationMethod};
		my $cuffnormSeed=$cuffnormParams->[0]->{seed};
		my $cuffnormNormalization=$cuffnormParams->[0]->{normalization};
		my $cuffdiffUseCuffmergeAssembly=$cuffdiffParams->[0]->{useCuffmergeAssembly};		
		my $cuffdiffNThreads=$cuffdiffParams->[0]->{nThreads};
		my $cuffdiffFragBiasCorrect=$cuffdiffParams->[0]->{fragBiasCorrect};		
		my $cuffdiffMultiReadCorrect=$cuffdiffParams->[0]->{multiReadCorrect};	
		my $cuffdiffLibraryNormalizationMethod=$cuffdiffParams->[0]->{libraryNormalizationMethod};
		my $cuffdiffFDR=$cuffdiffParams->[0]->{FDR};
		my $cuffdiffMinAlignmentCount=$cuffdiffParams->[0]->{minAlignmentCount};
		my $cuffdiffSeed=$cuffdiffParams->[0]->{seed};
		my $FPKMthreshold=$cuffdiffParams->[0]->{FPKMthreshold};
		my $cuffnormLabels="";
		my $cuffquantMaxBundleFrags=$cuffquantParams->[0]->{maxBundleFrags};
		my $cuffdiffMaxBundleFrags=$cuffdiffParams->[0]->{maxBundleFrags};		
		my $cuffdiff_noEffectiveLengthCorrection=$cuffdiffParams->[0]->{noEffectiveLengthCorrection};
		my $cuffdiff_noLengthCorrection=$cuffdiffParams->[0]->{noLengthCorrection};
		my $cuffdiff_dispersionMethod=$cuffdiffParams->[0]->{dispersionMethod};

		#to keep track of the original GTF
		my $originalGTF=$GTF;
		
		# REGULAR cuffquant EXECUTION (ALWAYS DONE)
		# When there are NO spikes, this execution is valid for both cuffdiff and cuffnorm (the same GTF file is valid for both
		# cuffdiff and cuffnorm, so cuffquant is excecuted ONLY HERE) 
		# When there are spikes available, this is the FIRST execution of cuffquant, that produces the proper files for cuffnorm,
		# using a GTF file with genes+spikes
		print STDOUT "\n[DOING] transcript quantification and differential expression\n";
		my $command="perl ".$Bin."cuffquant_cuffdiff_cuffnorm.pl --mode 1 --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --cuffquantNThreads ".$cuffquantNThreads;
		$command.=" --cuffquantFragBiasCorrect ".$cuffquantFragBiasCorrect." --cuffquantMultiReadCorrect ".$cuffquantMultiReadCorrect;
		$command.=" --alignmentsDir ".$alignmentsDir." --samtoolsPath ".$samtoolsPath." --cufflinksPath ".$cufflinksPath;
		$command.=" --libraryType ".$libraryType." --referenceSequence ".$referenceSequence." --extraPathsRequired ".$extraPathsRequired." --cuffquantMaxBundleFrags ".$cuffquantMaxBundleFrags;
		$command.=" --cuffquant_noLengthCorrection ".$cuffquant_noLengthCorrection;
		$command.=" --cuffquant_noEffectiveLengthCorrection ".$cuffquant_noEffectiveLengthCorrection;
		$command.=" --cuffdiffMaxBundleFrags ".$cuffdiffMaxBundleFrags;
		
		#if a GTF file is provided or if the user would like to use the GTF file created by cuffmerge
		if(($GTF ne "") || ($cuffquantUseCuffmergeAssembly eq "true")){
			
			#if neither the GTF file is provided nor the user wants to use the cuffmerge derived GTF
			if((!-e $GTF) && ($cuffquantUseCuffmergeAssembly eq "false")){
				print STDERR "\n[ERROR cuffquant]: $GTF file does not exist\n\n";
				exit(-1);				
			}
			
			#in case that the user wants to use the cuffmerge assembly			
			if($cuffquantUseCuffmergeAssembly eq "true"){
				if(-e $cuffmergeOutDir."merged.gtf"){
					$GTF=$cuffmergeOutDir."merged.gtf";
				}else{
					print STDERR "\n[ERROR cuffquant]: ".$cuffmergeOutDir."merged.gtf file does not exist\n\n";
					exit(-1);
				}
			}			
			
			$command.=" --GTF ".$GTF;
		}
			
		$command.=" --cuffquantOutDir ".$cuffquantOutDir;		
		$command.=" --cuffquantSeed ".$cuffquantSeed;
		$command.=" --samples ".$ALLsamples." >> ".$workspace."cuffquant.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"cquan".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/cuffdiffAndCuffnormAreWaiting",$command,$wait,$multiCFlag);

		if(-e $executionCreatedTempDir."/cuffdiffAndCuffnormAreWaiting"){
			$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/cuffdiffAndCuffnormAreWaiting",$queueSystem);	
		}		


		#**************************************************************************************************************************
		# 2nd CUFFquant EXECUTION (ONLY when having SPIKES, done to get rid of them in the annotation GTF file):
		# in case of having spikes, cuffquant must be run again, but NOW with the INITIAL GTF file that does NOT contain SPIKES
		# the cuffquant files produced here without considering spikes, are going to be used by CUFFDIFF (but not by CUFFNORM)
		$GTF=$originalGTF;
		my $cuffquantOutDirWithoutConsideringSpikes=$workspace."cuffquantOutDirWithoutConsideringSpikes_ForCuffdiffOnly/";
		if($doSpikesAndGenomeRefIndexing eq "true"){
			#if true => EXECUTES cuffquant again, but without spikes (to be used with cuffdiff ONLY)			
			File::Path::make_path($cuffquantOutDirWithoutConsideringSpikes);
			
			$command="perl ".$Bin."cuffquant_cuffdiff_cuffnorm.pl --mode 1 --nThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep." --cuffquantNThreads ".$cuffquantNThreads;
			$command.=" --cuffquantFragBiasCorrect ".$cuffquantFragBiasCorrect." --cuffquantMultiReadCorrect ".$cuffquantMultiReadCorrect;
			$command.=" --alignmentsDir ".$alignmentsDir." --samtoolsPath ".$samtoolsPath." --cufflinksPath ".$cufflinksPath;
			$command.=" --libraryType ".$libraryType." --referenceSequence ".$referenceSequence." --extraPathsRequired ".$extraPathsRequired." --cuffquantMaxBundleFrags ".$cuffquantMaxBundleFrags;
			$command.=" --cuffquant_noEffectiveLengthCorrection ".$cuffquant_noEffectiveLengthCorrection." --cuffdiffMaxBundleFrags ".$cuffdiffMaxBundleFrags;
			
			#if a GTF file is provided or if the user would like to use the GTF file created by cuffmerge
			if(($GTF ne "") || ($cuffquantUseCuffmergeAssembly eq "true")){

				#if neither the GTF file is provided nor the user wants to use the cuffmerge derived GTF
				if((!-e $GTF) && ($cuffquantUseCuffmergeAssembly eq "false")){
					print STDERR "\n[ERROR cuffquant]: $GTF file does not exist\n\n";
					exit(-1);				
				}

				#in case that the user wants to use the cuffmerge assembly			
				if($cuffquantUseCuffmergeAssembly eq "true"){
					if(-e $cuffmergeOutDir."merged.gtf"){
						$GTF=$cuffmergeOutDir."merged.gtf";
					}else{
						print STDERR "\n[ERROR cuffquant]: ".$cuffmergeOutDir."merged.gtf file does not exist\n\n";
						exit(-1);
					}
				}

				$command.=" --GTF ".$initialGTF;					
			}

			$command.=" --cuffquantOutDir ".$cuffquantOutDirWithoutConsideringSpikes;		
			$command.=" --cuffquantSeed ".$cuffquantSeed;
			$command.=" --samples ".$ALLsamples." >> ".$workspace."cuffquant.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"cquan".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/cuffdiffAndCuffnormAreWaiting",$command,$wait,$multiCFlag);

			if(-e $executionCreatedTempDir."/cuffdiffAndCuffnormAreWaiting"){
				$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/cuffdiffAndCuffnormAreWaiting",$queueSystem);	
			}
		}#if($doSpikesAndGenomeRefIndexing eq "true")
		#************************************************************************************************************************

		
		#<---- cuffdiff && cuffnorm ---->
		#<---- READS COMPARISONS and inserts them in the proper order for cuffdiff ---->
		
		# it loops for 100 posible positions max, more than required, just in case...

		
		#the loop is a way to set the order of the conditions for cuffdiff, according to what is
		#pointed out in the xml file
		
		$GTF=$originalGTF;
		my($cuffnormInputFiles);
		my @colors=("lightgreen","lightpink","lightblue","lightyellow","ligthred",
		"brightgreen","brightpink","brightblue","brightyellow","ligthred",
		"black","red",",green","yellow","blue","magenta","cyan");
		
		#these variables below are required to assign similar colors to all the samples of the same group/condition,
		#when performing the last cuffnorm comparison (with all the samples at once)
		my %sample_class=();
		my $sample_class_counter=0;
		my $samples_already_included="";
					
		foreach my $comparison (keys %$comparisons){
			my $cuffdiffComparisonLabels="";
			my $cuffdiffInputFiles="";
			$cuffnormInputFiles="";							
			my $comparisonName=${\$comparison};

			#print @$comparisons."\n";
			my $conditions=$comparisons->{$comparison}->{condition};		

			my $rememberOneLibraryName="";
			my $i=1;
			
			#finds out the number of conditions
			my @howManyConditions=keys %$conditions; 
			
			my $cuffnormColors="";
			
			#while it hasn't checked all the conditions
			while($i<=@howManyConditions){
				foreach my $condition(sort keys %$conditions){
					#condition name					
					#print "cond: ".${\$condition}."\n";

					my $cuffdiffPosition=$conditions->{$condition}->{cuffdiffPosition};					
					#print $cuffdiffPosition."\n";

					if($cuffdiffPosition==$i){ #includes this condition in the proper order
						$cuffdiffComparisonLabels.=${\$condition}.","; #condition name
						my $libraries=$conditions->{$condition}->{libraryName};
						
						my $are_new_samples_added=0;
						
						foreach my $library(@$libraries){
							if($doSpikesAndGenomeRefIndexing eq "false"){ #uses regular cuffquant files
								$cuffdiffInputFiles.=$cuffquantOutDir.$library."/abundances.cxb,"; #library (replicate)
							}else{ #uses additional cuffquant files without spikes info
								$cuffdiffInputFiles.=$cuffquantOutDirWithoutConsideringSpikes.$library."/abundances.cxb,"; #library (replicate)
							}
							
							$cuffnormInputFiles.=$cuffquantOutDir.$library."/abundances.cxb,"; #library (replicate)
							$cuffnormColors.=$colors[$i-1].",";
							$cuffnormLabels.=$library.",";
							$rememberOneLibraryName=$library;
						
						}
						
						if($are_new_samples_added){
							$sample_class_counter++;						
						}
						

						#moves to the next sample: deletes last "," and inserts an extra ':' character
						#to separate files from different samples
						$cuffdiffInputFiles=~ s/\,$//;
						$cuffdiffInputFiles.=":";
						$cuffnormInputFiles=~ s/\,$//;
						$cuffnormInputFiles.=":";
						
						
						$i++;
					}
				}
				
				
			}

			$cuffdiffComparisonLabels=~ s/\,$//;
			$cuffdiffInputFiles=~ s/\,$//;
			$cuffdiffInputFiles=~ s/\:$//;
			$cuffnormLabels=~ s/\,$//;
			$cuffnormInputFiles=~ s/\,$//;
			$cuffnormInputFiles=~ s/\:$//;
			
			#print $cuffdiffComparisonLabels."\n";
			#print $cuffdiffInputFiles."\n";

			#CUFFDIFF executes one comparison each time it runs
			$command="perl ".$Bin."cuffquant_cuffdiff_cuffnorm.pl --mode 2 --cuffdiffNThreads ".$cuffdiffNThreads;
			$command.=" --samtoolsPath ".$samtoolsPath." --cufflinksPath ".$cufflinksPath;
			$command.=" --cuffdiffMinAlignmentCount ".$cuffdiffMinAlignmentCount;
			$command.=" --cuffdiffFragBiasCorrect ".$cuffdiffFragBiasCorrect;
			$command.=" --referenceSequence ".$referenceSequence;
			$command.=" --cuffdiffMultiReadCorrect ".$cuffdiffMultiReadCorrect;
			$command.=" --cuffdiffMaxBundleFrags ".$cuffdiffMaxBundleFrags;
			$command.=" --cuffdiff_noEffectiveLengthCorrection ".$cuffdiff_noEffectiveLengthCorrection;
			$command.=" --cuffdiff_noLengthCorrection ".$cuffdiff_noLengthCorrection;
			$command.=" --cuffdiff_dispersionMethod ".$cuffdiff_dispersionMethod;

			#if a GTF file is provided or if the user would like to use the GTF file created by cuffmerge
			if(($GTF ne "") || ($cuffdiffUseCuffmergeAssembly eq "true")){

				#if neither the GTF file is provided nor the user wants to use the cuffmerge derived GTF
				if((!-e $GTF) && ($cuffdiffUseCuffmergeAssembly eq "false")){
					print STDERR "\n[ERROR cuffdiff]: $GTF file does not exist\n\n";
					exit(-1);				
				}

				#in case that the user wants to use the cuffmerge assembly			
				if($cuffdiffUseCuffmergeAssembly eq "true"){
					if(-e $cuffmergeOutDir."merged.gtf"){
						$GTF=$cuffmergeOutDir."merged.gtf";
					}else{
						print STDERR "\n[ERROR cuffdiff]: ".$cuffmergeOutDir."merged.gtf file does not exist\n\n";
						exit(-1);
					}
				}

				if($doSpikesAndGenomeRefIndexing eq "false"){
					$command.=" --GTF ".$GTF;
				}else
				{
					# when having spikes, it is more appropriate to not consider them for cuffdiff as they could affect
					# normalization values (FPKMs) for regular genes. So in this case, the original GTF is given instead of the
					# one with the combined annotation (genes+spikes)
					$command.=" --GTF ".$initialGTF;
				}
			}			


			#regarding the library type, it has to assume that all the libraries in this comparison are of the same library type.
			#so it only checks the type of one of them, it should be enough (different types cannot be mixed with cuffdiff)
			my @libs=split(',',$ALLsamples);		
			my @libType=split(',',$libraryType);

			my $cuffdiffLibraryType="";
			for(my $k=0;$k<@libs;$k++){			
				if($libs[$k] eq $rememberOneLibraryName){
					$cuffdiffLibraryType=$libType[$k];
					last;
				}
			}

			$cuffdiffLibraryType=~ s/ $//;
			
			$command.=" --cuffdiffLibraryNormalizationMethod ".$cuffdiffLibraryNormalizationMethod;
			$command.=" --cuffdiffLibraryType ".$cuffdiffLibraryType;
			$command.=" --cuffdiffSeed ".$cuffdiffSeed;
			$command.=" --cuffdiffFDR ".$cuffdiffFDR;
			$command.=" --FPKMthreshold ".$FPKMthreshold;
			$command.=" --cuffdiffOutDir ".$cuffdiffOutDir.$comparisonName."/";
			$command.=" --cuffdiffComparisonLabels ".$cuffdiffComparisonLabels;				
			$command.=" --cuffdiffInputFiles ".$cuffdiffInputFiles;
			$command.=" --extraPathsRequired ".$extraPathsRequired;
			$command.=" >> ".$workspace."cuffdiff.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";				

			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"cdif".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">>".$executionCreatedTempDir."/GSEAIsWaiting",$command,$wait,$multiCFlag);

		
			#CUFFNORM executes the samples that belong to one comparison each time it runs
			$GTF=$originalGTF;
			
			$command="perl ".$Bin."cuffquant_cuffdiff_cuffnorm.pl --mode 3 --cuffnormNThreads ".$cuffnormNThreads;
			$command.=" --samtoolsPath ".$samtoolsPath." --cufflinksPath ".$cufflinksPath;
			$command.=" --cuffnormOutputFormat ".$cuffnormOutputFormat;
			$command.=" --cuffnormlibraryNormalizationMethod ".$cuffnormlibraryNormalizationMethod;
			$command.=" --cuffnormSeed ".$cuffnormSeed;
			$command.=" --cuffnormNormalization ".$cuffnormNormalization;
			$command.=" --colors ".$cuffnormColors;

			#if a GTF file is provided or if the user would like to use the GTF file created by cuffmerge
			if(($GTF ne "") || ($cuffnormUseCuffmergeAssembly eq "true")){

				#if neither the GTF file is provided nor the user wants to use the cuffmerge derived GTF
				if((!-e $GTF) && ($cuffnormUseCuffmergeAssembly eq "false")){
					print STDERR "\n[ERROR cuffnorm]: $GTF file does not exist\n\n";
					exit(-1);				
				}

				#in case that the user wants to use the cuffmerge assembly			
				if($cuffnormUseCuffmergeAssembly eq "true"){
					if(-e $cuffmergeOutDir."merged.gtf"){
						$GTF=$cuffmergeOutDir."merged.gtf";
					}else{
						print STDERR "\n[ERROR cuffnorm]: ".$cuffmergeOutDir."merged.gtf file does not exist\n\n";
						exit(-1);
					}
				}

				$command.=" --GTF ".$GTF;
			}
			
			$command.=" --cuffnormLibraryType ".$cuffdiffLibraryType;
			$command.=" --cuffnormOutDir ".$cuffnormOutDir.$comparisonName."/";
			#$command.=" --cuffnormLabels ".$cuffnormLabels; #uses values form cuffdiff
			$command.=" --cuffnormLabels ".$cuffdiffComparisonLabels; #uses values form cuffdiff
			$command.=" --cuffnormInputFiles ".$cuffnormInputFiles;			
			$command.=" --extraPathsRequired ".$extraPathsRequired;
			$command.=" --tmpDir ".$executionCreatedTempDir;
			$command.=" >> ".$workspace."cuffnorm.log 2>&1";
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";				

			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"cnor".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);


		}#foreach my $comparison (keys %$comparisons)
		
		
		#EXECUTES the last CUFFNORM for ALL SAMPLES together		
		$GTF=$originalGTF;		
		
		my @unsorted_Samples=split(',',$ALLsamples);
		
		#samples are sorted here (they are recevied in some kind of random order)
		my @sorted_samples = sort { lc($a) cmp lc($b) } @unsorted_Samples;
		
		my $ALL_LABELS="";
		foreach my $label(@sorted_samples){
			$ALL_LABELS.=$label.",";
		}
		$ALL_LABELS=~ s/\,$//;
		
		$cuffnormInputFiles="";
		my $colorsForAllSamples="";
		#my $counter=0;
		foreach my $lib(@sorted_samples){
			$cuffnormInputFiles.=$cuffquantOutDir.$lib."/abundances.cxb:";
			
			#COLOR patterns for all sample PCAs
			
			#OPTION 1: assigns a color from @colors corresponding to sample number
			#$colorsForAllSamples.=$colors[$counter].",";
			#$counter++;
			
			#OPTION 2: assigns a RANDOM color for each sample from @colors: 17 colors in total
			#$colorsForAllSamples.=$colors[int(rand(17))].",";
			
			#OPTION 3: all samples in black color
			$colorsForAllSamples.="black,"; # assigns black color to all samples
		}
		
		$cuffnormInputFiles=~ s/\,$//;
		$cuffnormInputFiles=~ s/\:$//;
		
		#cuffnorm has to assume that all the samples are of the same type, because they are all normalized together
		my $onelibType=(split(',',$libraryType))[0];
				
		$command="perl ".$Bin."cuffquant_cuffdiff_cuffnorm.pl --mode 3 --cuffnormNThreads ".$cuffnormNThreads;
		$command.=" --samtoolsPath ".$samtoolsPath." --cufflinksPath ".$cufflinksPath;
		$command.=" --cuffnormOutputFormat ".$cuffnormOutputFormat;
		$command.=" --cuffnormlibraryNormalizationMethod ".$cuffnormlibraryNormalizationMethod;
		$command.=" --cuffnormSeed ".$cuffnormSeed;
		$command.=" --cuffnormNormalization ".$cuffnormNormalization;
		$command.=" --colors ".$colorsForAllSamples;

		#if a GTF file is provided or if the user would like to use the GTF file created by cuffmerge
		if(($GTF ne "") || ($cuffnormUseCuffmergeAssembly eq "true")){

			#if neither the GTF file is provided nor the user wants to use the cuffmerge derived GTF
			if((!-e $GTF) && ($cuffnormUseCuffmergeAssembly eq "false")){
				print STDERR "\n[ERROR cuffnorm]: $GTF file does not exist\n\n";
				exit(-1);				
			}

			#in case that the user wants to use the cuffmerge assembly			
			if($cuffnormUseCuffmergeAssembly eq "true"){
				if(-e $cuffmergeOutDir."merged.gtf"){
					$GTF=$cuffmergeOutDir."merged.gtf";
				}else{
					print STDERR "\n[ERROR cuffnorm]: ".$cuffmergeOutDir."merged.gtf file does not exist\n\n";
					exit(-1);
				}
			}

			$command.=" --GTF ".$GTF;
		}
		
		$command.=" --cuffnormLibraryType ".$onelibType;
		$command.=" --cuffnormOutDir ".$cuffnormOutDir."ALLsamples/";
		$command.=" --cuffnormLabels ".$ALL_LABELS;
		$command.=" --cuffnormInputFiles ".$cuffnormInputFiles; 
		$command.=" --extraPathsRequired ".$extraPathsRequired;
		$command.=" --tmpDir ".$executionCreatedTempDir;
		$command.=" >> ".$workspace."cuffnorm.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";				

		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"cnor".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);


	}else{
		print STDERR "\n[ERROR]: problem with directories:\n";
		if(!-d $cuffquantOutDir){print STDERR "\n[ERROR]: $cuffquantOutDir directory does not exist\n";}
		if(!-d $cuffdiffOutDir){print STDERR "\n[ERROR]: $cuffdiffOutDir directory does not exist\n";}
		if(!-d $cuffnormOutDir){print STDERR "\n[ERROR]: $cuffnormOutDir directory does not exist\n";}
		if(!-d $alignmentsDir){print STDERR "\n[ERROR]: $alignmentsDir directory does not exist\n";}
		exit(-1);
	}		
}

########## runs htseq-count (gets read counts for genes) + DESeq2 differential expression
sub level_6{
    	my($perl5lib,$comparisons,$deseqParams,$extraPathsRequired,$htseqcountPath,$htseqCountPythonpath,$htseqcountParams,$samtoolsPath,$samples,$GTF,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,
	$workspace,$experimentName,$logfh,$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;

	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	if(-e $executionCreatedTempDir."/htseqCountIsWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/htseqCountIsWaiting",$queueSystem);	
	}
	
	# <----------------- only for testing ----------------------->	 
#	use Data::Dumper;
#	print STDERR "FILE:\n".Dumper($comparisons)."\n";
	

	# <----------------- htseq-count ----------------------->
	
	my $alignmentsDir=$workspace."alignments/";
	my $htseqcountOutDir=$workspace."htseqCount/";
	my $deseqOutDir=$workspace."deseq/";
	File::Path::make_path($htseqcountOutDir);
	File::Path::make_path($deseqOutDir);
	
	
	
	if(-d $alignmentsDir && -d $htseqcountOutDir){				
		my $ALLsamples="";
		my $libraryType="";		
		
		foreach my $key (keys %$samples){
			#sample names
			$ALLsamples.=${\$key}.",";
			#library type
			$libraryType.=$samples->{$key}{libraryType}[0].",";
		}
		
		$ALLsamples=~ s/,$//;
		$libraryType=~ s/,$//;

		my $minaqual=$htseqcountParams->[0]->{minaqual};
		my $featuretype=$htseqcountParams->[0]->{featuretype};		
		my $idattr=$htseqcountParams->[0]->{idattr};
		my $mode=$htseqcountParams->[0]->{mode}[0];
		
		print STDOUT "\n[DOING] counting reads with Htseq-count\n";
		my $command="perl ".$Bin."htseqCount.pl";
		$command.=" --extraPathsRequired ".$extraPathsRequired;
		$command.=" --perl5lib ".$perl5lib;
		$command.=" --mode ".$mode;
		$command.=" --minaqual ".$minaqual;
		$command.=" --featuretype ".$featuretype;
		$command.=" --idattr ".$idattr;
		$command.=" --htseqcountPath ".$htseqcountPath;
		$command.=" --htseqcountPythonpath ".$htseqCountPythonpath;
		$command.=" --samtoolsPath ".$samtoolsPath;
		$command.=" --samples ".$ALLsamples;
		$command.=" --libraryType ".$libraryType;
		$command.=" --GTF ".$GTF;
		$command.=" --numberOfThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep;
		$command.=" --alignmentsDir ".$alignmentsDir;
		$command.=" --htseqCountOutDir ".$htseqcountOutDir;
		$command.=" --executionCreatedTempDir ".$executionCreatedTempDir;
		$command.=" >> ".$workspace."htseqCount.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";				
		
		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"htsq".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/deseqIsWaiting",$command,$wait,$multiCFlag);
		
		#wait until HTseqcount is running
		if(-e $executionCreatedTempDir."/deseqIsWaiting"){
			$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/tophatIsWaiting",$queueSystem);	
		}
		
		
		my $listWithAllSamples="";
		my $allComparisonLables="all_1,all_2";
		my $samplesIncludedIn_listWithAllSamples="";
		
		foreach my $comparison (keys %$comparisons){
			my $deseqComparisonLabels="";
			my $deseqInputFiles="";				
			my $comparisonName=${\$comparison};

			#print @$comparisons."\n";
			my $conditions=$comparisons->{$comparison}->{condition};		

			my $rememberOneLibraryName="";
			
			my $i=1;
			
			#finds out the number of conditions
			my @howManyConditions=keys %$conditions; 
			#while it hasn't checked all the conditions
			while($i<=@howManyConditions){
				foreach my $condition(sort keys %$conditions){
					#condition name					
					#print "cond: ".${\$condition}."\n";


					my $cuffdiffPosition=$conditions->{$condition}->{cuffdiffPosition};					
					#print $cuffdiffPosition."\n";

					if($cuffdiffPosition==$i){ #includes this condition in the proper order
						$deseqComparisonLabels.=${\$condition}.","; #condition name
						my $libraries=$conditions->{$condition}->{libraryName};						

						foreach my $library(@$libraries){
							$deseqInputFiles.=$htseqcountOutDir.$library.".xls,"; #library (replicate)
							$rememberOneLibraryName=$library;
							
							if($samplesIncludedIn_listWithAllSamples!~ /=$library=/){
								$listWithAllSamples.=$htseqcountOutDir.$library.".xls,"; #library (replicate)
								$samplesIncludedIn_listWithAllSamples.="=".$library."=";								
							}
						}

						#moves to the next sample: deletes last "," and inserts an extra ':' character
						#to separate files from different samples
						$deseqInputFiles=~ s/\,$//;
						$deseqInputFiles.=":";
						
						#for adding ':' only once, to force two (artificial) groups,
						#no matter what groups they are						
						if($listWithAllSamples!~ /\:/){
							$listWithAllSamples=~ s/\,$//;
							$listWithAllSamples.=":";
						}

						$i++;
					}
				}
			}

			$deseqComparisonLabels=~ s/\,$//;
			$deseqInputFiles=~ s/\,$//;
			$deseqInputFiles=~ s/\:$//;
						
			
			my $nThreads=$deseqParams->[0]->{nThreads};
			my $alpha=$deseqParams->[0]->{alpha};
			my $pAdjustMethod=$deseqParams->[0]->{pAdjustMethod};			
			
			#executes each comparison in each iteration
			print STDOUT "[DOING] differential expression with DESeq2\n";
			my $command="perl ".$Bin."deseq.pl";
			$command.=" --nThreads ".$nThreads;
			$command.=" --GTF ".$GTF;			
			$command.=" --alpha ".$alpha;
			$command.=" --pAdjustMethod ".$pAdjustMethod;
			$command.=" --deseqOutDir ".$deseqOutDir;
			$command.=" --comparisonName ".$comparisonName;
			$command.=" --deseqInputFiles ".$deseqInputFiles;
			$command.=" --deseqComparisonLabels ".$deseqComparisonLabels; #uses values form cuffdiff
			$command.=" --extraPathsRequired ".$extraPathsRequired;
			$command.=" --tmpDir ".$executionCreatedTempDir;
			$command.=" >> ".$workspace."deseq.log 2>&1";
			
			print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

			queue::executeScript($queueSystem,$queueName,$queueSGEProject,"deseq2".substr($experimentName,0,5),
			$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);

			
		}#foreach my $comparison (keys %$comparisons)
		
		$listWithAllSamples=~ s/\,$//;
		
		#executes DESeq2 once more with all the samples together
		my $nThreads=$deseqParams->[0]->{nThreads};
		my $alpha=$deseqParams->[0]->{alpha};
		my $pAdjustMethod=$deseqParams->[0]->{pAdjustMethod};		
		$command="perl ".$Bin."deseq.pl";
		$command.=" --nThreads ".$nThreads;
		$command.=" --GTF ".$GTF;			
		$command.=" --alpha ".$alpha;
		$command.=" --pAdjustMethod ".$pAdjustMethod;
		$command.=" --deseqOutDir ".$deseqOutDir;
		$command.=" --comparisonName ALLsamples";
		$command.=" --deseqInputFiles ".$listWithAllSamples;
		$command.=" --deseqComparisonLabels ".$allComparisonLables;
		$command.=" --extraPathsRequired ".$extraPathsRequired;
		$command.=" --tmpDir ".$executionCreatedTempDir;
		$command.=" >> ".$workspace."deseq.log 2>&1";	
			
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"deseq2".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);

		#deletes differential expression files for ALL samples as it makes no sense (as more than two conditions could be included)
		#deletes the associated rnk file
		$command="rm -f ".$deseqOutDir."ALLsamples.differentialExpression* ".$deseqOutDir."ALLsamples.rnk";
		system($command);
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";
	
	}else{
		print STDERR "\n[ERROR]: problem with the following directories:\n";
		if(!-d $alignmentsDir){print STDERR "\n[ERROR]: $alignmentsDir directory does not exist\n";}
		if(!-d $htseqcountOutDir){print STDERR "\n[ERROR]: $htseqcountOutDir directory does not exist\n";}
		if(!-d $deseqOutDir){print STDERR "\n[ERROR]: $deseqOutDir directory does not exist\n";}
		exit(-1);
	}	
}

########## level 7: creates wiggle files from bam alignments
sub level_7{
	my($bedGraphToBigWigPath,$bedGraphToBigWigParams,$bedtoolsPath,$samtoolsPath,$samples,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,$experimentName,$logfh,
			$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;

	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	if(-e $executionCreatedTempDir."/bedGraphToBigWigIsWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/bedGraphToBigWigIsWaiting",$queueSystem);	
	}

	
	# <----------------- only for testing ----------------------->	 
#	use Data::Dumper;
#	print STDERR "FILE:\n".Dumper($comparisons)."\n";
	

	# <----------------- bigwig files are created ----------------------->
	
	my $alignmentsDir=$workspace."alignments/";
	my $bigwigFilesDir=$workspace."bigwigFilesDir/";
	File::Path::make_path($bigwigFilesDir);
	
	if(-d $alignmentsDir && -d $bigwigFilesDir){				
		my $ALLsamples="";
		my $libraryType="";		
		
		foreach my $key (keys %$samples){
			#sample names
			$ALLsamples.=${\$key}.",";
			#library type
			$libraryType.=$samples->{$key}{libraryType}[0].",";
		}
		
		$ALLsamples=~ s/,$//;
		$libraryType=~ s/,$//;

		my $chromosomeSizesFile=$bedGraphToBigWigParams->[0]->{chromosomeSizesFile};
		my $bigDataUrlPrefix=$bedGraphToBigWigParams->[0]->{bigDataUrlPrefix};
		
		print STDOUT "\n[DOING] bigWig files creation\n";
		my $command="perl ".$Bin."createbigwigFiles.pl";
		$command.=" --chromosomeSizesFile ".$chromosomeSizesFile;
		$command.=" --bigDataUrlPrefix ".$bigDataUrlPrefix;
		$command.=" --samtoolsPath ".$samtoolsPath;
		$command.=" --bedtoolsPath ".$bedtoolsPath;
		$command.=" --bedGraphToBigWigPath ".$bedGraphToBigWigPath;
		$command.=" --samples ".$ALLsamples;
		$command.=" --numberOfThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep;
		$command.=" --alignmentsDir ".$alignmentsDir;
		$command.=" --bigwigFilesDir ".$bigwigFilesDir;
		$command.=" --experimentName ".$experimentName;
		$command.=" >> ".$workspace."bigwigFilesDir.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";				

		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"bWig".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);

	
	}else{
		print STDERR "\n[ERROR]: problem with directories:\n";
		if(!-d $alignmentsDir){print STDERR "\n[ERROR]: $alignmentsDir directory does not exist\n";}
		if(!-d $bigwigFilesDir){print STDERR "\n[ERROR]: $bigwigFilesDir directory does not exist\n";}
		exit(-1);
	}		
}

########## level 8: runs GSEA for all rnk files with the provided geneset files
sub level_8{
	my($gseaChip,$gseamaxMemory,$gseaPath,$gseaParams,$workspace,$experimentName,$logfh,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$comparisons,
	$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;
    
	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	if(-e $executionCreatedTempDir."/GSEAIsWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/GSEAIsWaiting",$queueSystem);	
	}
	
	# <----------------- only for testing ----------------------->
#	use Data::Dumper;
#	print STDERR "FILE:\n".Dumper($gseaParams)."\n";	


    	#GSEA FOR CUFFDIFF RESULTS
	my $cuffdiffOutDir=$workspace."cuffdiff/";	
	
	if(-d $cuffdiffOutDir){

		my $GSEAoutDir_Cuffdiff=$workspace."GSEA_Cuffdiff/";
		File::Path::make_path($GSEAoutDir_Cuffdiff);
		
		if(-d $GSEAoutDir_Cuffdiff){
			my $genesets=$gseaParams->[0]->{geneset};
			my $genesetFiles="";

			#creates a string with al the gene set files
			foreach my $gene(@$genesets){
				#if /^HASH/ means that there is no value in the XML, ie. <geneset></geneset>
				if(($gene ne "") && ($gene!~ /^HASH/)) {$genesetFiles.=$gene.","}
			}

			#if gene sets are provided
			if($genesetFiles ne ""){
				$genesetFiles=~ s/,$//;

				my $collapse=$gseaParams->[0]->{collapse};
				my $mode=$gseaParams->[0]->{mode};
				my $norm=$gseaParams->[0]->{norm};
				my $nperm=$gseaParams->[0]->{nperm};
        			my $scoring_scheme=$gseaParams->[0]->{scoring_scheme};
        			my $include_only_symbols=$gseaParams->[0]->{include_only_symbols};
        			my $make_sets=$gseaParams->[0]->{make_sets};
        			my $plot_top_x=$gseaParams->[0]->{plot_top_x};
        			my $rnd_seed=$gseaParams->[0]->{rnd_seed};
        			my $set_max=$gseaParams->[0]->{set_max};
        			my $set_min=$gseaParams->[0]->{set_min};
        			my $zip_report=$gseaParams->[0]->{zip_report};
        			my $geneset=$gseaParams->[0]->{geneset}[0];

				print STDOUT "\n[DOING] Gene Set Enrichment Analysis (GSEA) on Cuffdiff results\n";
				#reads comparisons files
				my $comparisonName="";
				foreach my $comparison (keys %$comparisons){
					$comparisonName=${\$comparison};			

					my $command="perl ".$Bin."GSEA.pl";
					$command.=" --gseaPath ".$gseaPath;
					$command.=" --gseaChip ".$gseaChip;
					$command.=" --gseamaxMemory ".$gseamaxMemory;
					$command.=" --collapse ".$collapse;
					$command.=" --mode ".$mode;
					$command.=" --norm ".$norm;
					$command.=" --nperm ".$nperm;
					$command.=" --scoring_scheme ".$scoring_scheme;
					$command.=" --include_only_symbols ".$include_only_symbols;
					$command.=" --make_sets ".$make_sets;
					$command.=" --plot_top_x ".$plot_top_x;
					$command.=" --rnd_seed ".$rnd_seed;
					$command.=" --set_max ".$set_max;
					$command.=" --set_min ".$set_min;
        				$command.=" --zip_report ".$zip_report;
        				$command.=" --genesets ".$genesetFiles;
        				$command.=" --GSEAoutDir ".$GSEAoutDir_Cuffdiff;
					$command.=" --diffExp_outDir ".$cuffdiffOutDir;
					$command.=" --numberOfThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep;
					$command.=" --comparison ".$comparisonName;
					$command.=" --tmpDir ".$executionCreatedTempDir;					
					$command.=" >> ".$workspace."GSEA.log 2>&1";
					print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

					queue::executeScript($queueSystem,$queueName,$queueSGEProject,"GSEA".substr($experimentName,0,5),
					$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);
				}		
			}else{
				print STDERR "\n[ERROR]: GSEA for Cuffdiff results cannot be executed because no gene sets has been provided\n";
				print $logfh (Miscellaneous->getCurrentDateAndTime())."[ERROR]: GSEA for Cuffdiff results cannot be executed because no gene sets has been provided\n";
			}
		}else{
			print STDOUT "\n[ERROR]: $GSEAoutDir_Cuffdiff directory does not exist, skipping GSEA for Cuffdiff results\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())."[ERROR]: $GSEAoutDir_Cuffdiff directory does not exist, skipping GSEA for Cuffdiff results\n";			
		
		}
        
	}else{
		print STDOUT "\n[WARNING]: Skipping GSEA for Cuffdiff results as $cuffdiffOutDir directory does not exist. Probably you haven't run differential expression with Cuffdiff\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[WARNING]: Skipping GSEA for Cuffdiff results as $cuffdiffOutDir directory does not exist. Probably you haven't run differential expression with Cuffdiff\n";		
	}
	
	
	#GSEA FOR DESeq2 RESULTS	
	my $deseqOutDir=$workspace."deseq/";
	if(-d $deseqOutDir){

		my $GSEAoutDir_deseq=$workspace."GSEA_DESeq/";
		File::Path::make_path($GSEAoutDir_deseq);
		
		if(-d $GSEAoutDir_deseq){
			my $genesets=$gseaParams->[0]->{geneset};
			my $genesetFiles="";

			#creates a string with al the gene set files
			foreach my $gene(@$genesets){
				#if /^HASH/ means that there is no value in the XML, ie. <geneset></geneset>
				if(($gene ne "") && ($gene!~ /^HASH/)) {$genesetFiles.=$gene.","}
			}

			#if gene sets are provided
			if($genesetFiles ne ""){
				$genesetFiles=~ s/,$//;

				my $collapse=$gseaParams->[0]->{collapse};
				my $mode=$gseaParams->[0]->{mode};
				my $norm=$gseaParams->[0]->{norm};
				my $nperm=$gseaParams->[0]->{nperm};
        			my $scoring_scheme=$gseaParams->[0]->{scoring_scheme};
        			my $include_only_symbols=$gseaParams->[0]->{include_only_symbols};
        			my $make_sets=$gseaParams->[0]->{make_sets};
        			my $plot_top_x=$gseaParams->[0]->{plot_top_x};
        			my $rnd_seed=$gseaParams->[0]->{rnd_seed};
        			my $set_max=$gseaParams->[0]->{set_max};
        			my $set_min=$gseaParams->[0]->{set_min};
        			my $zip_report=$gseaParams->[0]->{zip_report};
        			my $geneset=$gseaParams->[0]->{geneset}[0];

				print STDOUT "\n[DOING] Gene Set Enrichment Analysis (GSEA) on DESeq2 results\n";
				#reads comparisons files
				my $comparisonName="";
				foreach my $comparison (keys %$comparisons){
					$comparisonName=${\$comparison};			

					my $command="perl ".$Bin."GSEA.pl";
					$command.=" --gseaPath ".$gseaPath;
					$command.=" --gseaChip ".$gseaChip;
					$command.=" --gseamaxMemory ".$gseamaxMemory;
					$command.=" --collapse ".$collapse;
					$command.=" --mode ".$mode;
					$command.=" --norm ".$norm;
					$command.=" --nperm ".$nperm;
					$command.=" --scoring_scheme ".$scoring_scheme;
					$command.=" --include_only_symbols ".$include_only_symbols;
					$command.=" --make_sets ".$make_sets;
					$command.=" --plot_top_x ".$plot_top_x;
					$command.=" --rnd_seed ".$rnd_seed;
					$command.=" --set_max ".$set_max;
					$command.=" --set_min ".$set_min;
        				$command.=" --zip_report ".$zip_report;
        				$command.=" --genesets ".$genesetFiles;
        				$command.=" --GSEAoutDir ".$GSEAoutDir_deseq;
					$command.=" --diffExp_outDir ".$deseqOutDir;
					$command.=" --numberOfThreads ".$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep;
					$command.=" --comparison ".$comparisonName;
					$command.=" --tmpDir ".$executionCreatedTempDir;
					$command.=" >> ".$workspace."GSEA.log 2>&1";
					print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";

					queue::executeScript($queueSystem,$queueName,$queueSGEProject,"GSEA".substr($experimentName,0,5),
					$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);
				}		
			}else{
				print STDERR "\n[ERROR]: GSEA for DESeq2 results cannot be executed because no gene sets has been provided\n";
				print $logfh (Miscellaneous->getCurrentDateAndTime())."[ERROR]: GSEA for DESeq2 results cannot be executed because no gene sets has been provided\n";
			}
		}else{
			print STDOUT "\n[ERROR]: $GSEAoutDir_deseq directory does not exist, skipping GSEA for DESeq2 results\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())."[ERROR]: $GSEAoutDir_deseq directory does not exist, skipping GSEA for DESeq2 results\n";			
		
		}
        
	}else{
		print STDOUT "\n[WARNING]: Skipping GSEA for DESeq2 results as $deseqOutDir directory does not exist. Probably you haven't run differential expression with DESeq2\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[WARNING]: Skipping GSEA for DESeq2 results as $deseqOutDir directory does not exist. Probably you haven't run differential expression with DESeq2\n";		
	}
	
    
    
}

########## level 9: gene fusion prediction with Tophat-fusion
sub level_9{
	my($tophatfusionParams,$tophatPath,$bowtiePath,$samtoolsPath,$workspace,$experimentName,$logfh,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$samples,$referenceSequence,
	$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueSGEProject)=@_;
    
    
	my $multiCFlag = queue::multicoreFlag($queueSystem, $multicore); # for us '-pe' => qsub -pe multicore 4 -P Experiment -q ngs
	my $wait=""; #empty=don't wait

	if(-e $executionCreatedTempDir."/tophatFusionIsWaiting"){
		$wait = queue::waitPidsQUEUE($executionCreatedTempDir."/tophatFusionIsWaiting",$queueSystem);	
	}
	
	# <----------------- only for testing ----------------------->
#	use Data::Dumper;
#	print STDERR "FILE:\n".Dumper($gseaParams)."\n";	
    
	my $alignmentsDir=$workspace."alignments/";
	
	if(-d $alignmentsDir){
		my $ALLsamples="";				
		
		foreach my $key (keys %$samples){
			#sample names
			$ALLsamples.=${\$key}.",";
			#$ALLsamples.=$key.",";
		}
		
		$ALLsamples=~ s/,$//;

		my $numFusionReads=$tophatfusionParams->[0]->{numFusionReads};
		my $numFusionPairs=$tophatfusionParams->[0]->{numFusionPairs};
		my $numFusionBoth=$tophatfusionParams->[0]->{numFusionBoth};
		my $fusionReadMismatches=$tophatfusionParams->[0]->{fusionReadMismatches};
		my $fusionMultireads=$tophatfusionParams->[0]->{fusionMultireads};
		my $nonHuman=$tophatfusionParams->[0]->{nonHuman};
		my $pathToAnnotationFiles=$tophatfusionParams->[0]->{pathToAnnotationFiles};
		my $nTophatFusionThreads=$tophatfusionParams->[0]->{nTophatFusionThreads};
		my $pathToBlastAll=$tophatfusionParams->[0]->{pathToBlastAll};
		my $pathToBlastn=$tophatfusionParams->[0]->{pathToBlastn};

		#here $nThreads is no used as tophatFusion does it for all samples at once, hence
		#you don't specifiy sample names, and so there is no possibility of sending several threads at once
		print STDOUT "\n[DOING] fusion prediction\n";
		my $command="perl ".$Bin."fusions.pl";
		$command.=" --numFusionReads ".$numFusionReads;
		$command.=" --numFusionPairs ".$numFusionPairs;
		$command.=" --numFusionBoth ".$numFusionBoth;
		$command.=" --fusionReadMismatches ".$fusionReadMismatches;
		$command.=" --fusionMultireads ".$fusionMultireads;
		$command.=" --nonHuman ".$nonHuman;
		$command.=" --pathToAnnotationFiles ".$pathToAnnotationFiles;
		$command.=" --samples ".$ALLsamples;
		$command.=" --nTophatFusionThreads ".$nTophatFusionThreads;
		$command.=" --alignmentsDir ".$alignmentsDir;
		$command.=" --experimentName ".$experimentName;
		$command.=" --tophatPath ".$tophatPath;
		$command.=" --bowtiePath ".$bowtiePath;
		$command.=" --samtoolsPath ".$samtoolsPath;
		$command.=" --pathToBlastAll ".$pathToBlastAll;
		$command.=" --pathToBlastn ".$pathToBlastn;
		$command.=" --referenceSequence ".$referenceSequence;
		$command.=" --workspace ".$workspace;
		$command.=" >> ".$workspace."fusions.log 2>&1";
		print $logfh (Miscellaneous->getCurrentDateAndTime()).$command."\n";				

		queue::executeScript($queueSystem,$queueName,$queueSGEProject,"fus".substr($experimentName,0,5),
		$workspace.$experimentName."_error",$workspace.$experimentName."_error",">".$executionCreatedTempDir."/NOTHING",$command,$wait,$multiCFlag);		

	}else{
		print STDERR "\n[ERROR]: problem with directories:\n";
		if(!-d $alignmentsDir){print STDERR "\n[ERROR]: $alignmentsDir directory does not exist\n";}
		exit(-1);
	}

}


sub changeSampleNamesInCaseTheyWereTrimmedAndOrDownsampledBefore{
	#in case that trimming and/or downsampling are done over the samples, it finds out
	#what is going to be the derived fastq file for each sample that is going to be aligned by tophat:
	#it should be the original sample (in case of no trimming nor downsampling), or a new fastq file in case that
	#the file was trimmed, downsampled or both trimmed and downsampled (the fastq file is different from the original one)
	
	#***changes in $samples are applied globally, no return is needed
	
	my ($samples,$workspace,$pairedEnd)=@_;
	
	my $downsamplingDir=$workspace."downsampledSamples/";	
	my $trimmingDir=$workspace."trimmedSamples/";

	#checks that the directories were built, then do the rest
	if(-d $downsamplingDir && -d $trimmingDir){
		foreach my $key (keys %$samples){																   	   

			my $type=$samples->{$key}{type}[0];	
			#is this sample going to be trimmed?
			my $trimming=$samples->{$key}{trimming}[0]{do};
			my $downsampling=$samples->{$key}{downsampling}[0]{do};

			#the SAMPLE was trimmed and downsampled
			if($type eq "fastq" && $trimming eq "true" && $downsampling eq "true"){
				my $leftFastqFile=$samples->{$key}{leftFile};
				my @aux=split('\/',$leftFastqFile);
				my $leftFile=$aux[@aux-1]; #my $leftFile=$aux[length(@aux)-2];
				my $file=$downsamplingDir.$leftFile.".trimmed.downsampled";
				$samples->{$key}{leftFile}=$file;

				#if paired-end experiment
				if($pairedEnd eq "true"){
					my $rightFastqFile=$samples->{$key}{rightFile}[0];
					@aux=split('\/',$rightFastqFile);
					my $rightFile=$aux[@aux-1]; #my $rightFile=$aux[length(@aux)-2];
					$file=$downsamplingDir.$rightFile.".trimmed.downsampled";
					$samples->{$key}{rightFile}[0]=$file;

				}				
			}elsif($type eq "fastq" && $trimming eq "false" && $downsampling eq "true"){
				#it was downsampled but not trimmed			
				my $leftFastqFile=$samples->{$key}{leftFile};
				my @aux=split('\/',$leftFastqFile);
				my $leftFile=$aux[@aux-1]; #my $leftFile=$aux[length(@aux)-2];
				my $file=$downsamplingDir.$leftFile.".downsampled";
				$samples->{$key}{leftFile}=$file;

				#if paired-end experiment
				if($pairedEnd eq "true"){
					my $rightFastqFile=$samples->{$key}{rightFile}[0];
					@aux=split('\/',$rightFastqFile);
					my $rightFile=$aux[@aux-1]; #my $rightFile=$aux[length(@aux)-2];
					$file=$downsamplingDir.$rightFile.".downsampled";
					$samples->{$key}{rightFile}[0]=$file;

				}				
			}elsif($type eq "fastq" && $trimming eq "true" && $downsampling eq "false"){
 				 #it was trimmed but not downsampled
 				 my $leftFastqFile=$samples->{$key}{leftFile};
 				 my @aux=split('\/',$leftFastqFile);
 				my $leftFile=$aux[@aux-1]; #my $leftFile=$aux[length(@aux)-2];
 				 my $file=$trimmingDir.$leftFile.".trimmed";
 				 $samples->{$key}{leftFile}=$file;

 				 #if paired-end experiment
 				 if($pairedEnd eq "true"){
 					 my $rightFastqFile=$samples->{$key}{rightFile}[0];
 					 @aux=split('\/',$rightFastqFile);
 					my $rightFile=$aux[@aux-1]; #my $rightFile=$aux[length(@aux)-2];
 					 my $file=$trimmingDir.$rightFile.".trimmed";
	 				 $samples->{$key}{rightFile}[0]=$file;
 				 }
			}				 
		}#foreach my $key (keys %$samples)
	}
#*** este else no tiene sentido ya que ambos subdir de downsampling y trimming no van a existir si en el experimento
# no se ha hecho trimming ni downsampling, sin embargo obligaria a parar la ejecucion si entrase en el else
#	else{
#		print "\n[ERROR]: One or both of the following directories do not exist:\n\n";
#		print STDERR "$trimmingDir\n$downsamplingDir\n\n";
#		exit(-1);	
#	}
}


sub removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch{
#After running the initial cuffnorm for ALLsamples, it removes from 'ALLsamples.genes.fpkm_table.xls' those genes with expression values at background levels + flat pattern genes.
#Finally, it creates a new annotation GTF, without the removed genes, and runs cuffquant+cuffdiff+cuffnorm again
	
	##if using $initialGTF => spikes annotation is not considered (modified initial GTF with spikes info is not used here)
	my ($cuffnormParams,$initialGTF,$workspace,$logfh)=@_;
	

	# <----------------- only for testing ----------------------->	 
#	use Data::Dumper;
#	print STDERR "FILE:\n".Dumper($comparisons)."\n";
	

	# <----------------- cuffquant, cuffdiff, cuffnorm ----------------------->
	my $cuffnormOutDir=$workspace."cuffnorm/";
	my $backgroundFiltered_AND_flatPatternFiltered_DIR=$workspace."cuffdiff_backgroundFiltered_AND_flatPatternFiltered/";
	File::Path::make_path($backgroundFiltered_AND_flatPatternFiltered_DIR);

	my $filteringLogFile=$backgroundFiltered_AND_flatPatternFiltered_DIR."filtering.log.txt";
	
	open(LOGFILEOUT,">>",$filteringLogFile);
	print STDOUT "[removing background level genes and flat pattern genes for Cuffdiff]\n";
	print $logfh (Miscellaneous->getCurrentDateAndTime())."[removing background level genes and flat pattern genes for Cuffdiff]\n";	
	print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[removing background level genes and flat pattern genes for Cuffdiff]\n";	
	
	if(-d $cuffnormOutDir && -d $backgroundFiltered_AND_flatPatternFiltered_DIR){
		my $backgroundLevel=$cuffnormParams->[0]->{backgroundLevel};
		my $required_PercentageOfSamplesOverBackgroundLevel=$cuffnormParams->[0]->{requiredPercentageOfSamplesOverBackgroundLevel};
		my $IQRthresdhold_forFlatPatternFiltering=$cuffnormParams->[0]->{IQRthresdhold_forFlatPatternFiltering};
	
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[background expression threshold] ".$backgroundLevel."\n";
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[required percentage of samples for background expression threshold] ".$required_PercentageOfSamplesOverBackgroundLevel."\n";
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[IQR threshold for flat pattern filtering] ".$IQRthresdhold_forFlatPatternFiltering."\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[background expression threshold] ".$backgroundLevel."\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[required percentage of samples for background expression threshold] ".$required_PercentageOfSamplesOverBackgroundLevel."\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[IQR threshold for flat pattern filtering] ".$IQRthresdhold_forFlatPatternFiltering."\n";
		print STDOUT "[background expression threshold] ".$backgroundLevel."\n";
		print STDOUT "[required percentage of samples for background expression threshold] ".$required_PercentageOfSamplesOverBackgroundLevel."\n";
		print STDOUT "[IQR threshold for flat pattern filtering] ".$IQRthresdhold_forFlatPatternFiltering."\n";

		my $inputFileForBackgroundFiltering=$cuffnormOutDir."ALLsamples.genes.fpkm_table.xls";
		my $newTable_background_genes_OUT=$backgroundFiltered_AND_flatPatternFiltered_DIR."ALLsamples.genes.fpkm_table.background_genes_removed.xls";
		my $newGTF_without_background_AND_flatpattern_genes=$backgroundFiltered_AND_flatPatternFiltered_DIR."GTF_without_background_AND_flatpatternGenes.gtf";
		my $list_with_discardedBackgroundGenes=$backgroundFiltered_AND_flatPatternFiltered_DIR."background_level_removed_genes.txt";
		my $nTotalGenes=0;
		my $nSelectedGenes=0;
		
		system("rm -f ".$newTable_background_genes_OUT);
		system("rm -f ".$newGTF_without_background_AND_flatpattern_genes);
		system("rm -f ".$list_with_discardedBackgroundGenes);
		
		if(open(IN,"<",$inputFileForBackgroundFiltering)){
			open(OUT,">",$newTable_background_genes_OUT);
			open(DISCARDED,">",$list_with_discardedBackgroundGenes);
			
			my $nSamples=0;
			
			############ (1) Discards genes under background expression level ##################
			
			foreach my $line(<IN>){
				chomp($line);
				
				#prints out header line				
				if($line=~ /^Name/){
					print OUT $line."\n";
					
					#counts the number of total samples
					my @tokens=split('\t',$line);	
					
					for (my $j=1;$j<@tokens;$j++){
						$nSamples++;						
					}
					
					print $logfh (Miscellaneous->getCurrentDateAndTime())." [Number of total samples] ".$nSamples."\n";
					print STDOUT "[Number of total samples] ".$nSamples."\n";
					print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [Number of total samples] ".$nSamples."\n";
						
				}else{ #lines with data values
					$nTotalGenes++;
					my $nSamplesOverBackgroundLevel=0;								
					my @tokens=split('\t',$line);
					my $gene=$tokens[0];
										
					#reads expression values in all the samples for the current gene					
					for (my $j=1;$j<@tokens;$j++){
						my $value=$tokens[$j];
						
						if($value>$backgroundLevel){
							$nSamplesOverBackgroundLevel++;
						}						
					}					
					
					my $percentageOfSamplesOverBackgroundLevel=($nSamplesOverBackgroundLevel/$nSamples)*100;					
					
					# if the current gene is over background expression level for the required number of samples
					if($percentageOfSamplesOverBackgroundLevel>=$required_PercentageOfSamplesOverBackgroundLevel){
						#gene is selected as it has passed background expression level filter
						print OUT $line."\n";
						$nSelectedGenes++;
					}
					else{
						print DISCARDED $gene."\n";
					}			
				}			
			}
			
			close(OUT);
			close(DISCARDED);
			close(IN);

			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of total genes] ".$nTotalGenes."\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to background level filtering] ".($nTotalGenes-$nSelectedGenes)."\n";
			print STDOUT "[number of total genes] ".$nTotalGenes."\n";
			print STDOUT "[number of genes discarded due to background level filtering] ".($nTotalGenes-$nSelectedGenes)."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of total genes] ".$nTotalGenes."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to background level filtering] ".($nTotalGenes-$nSelectedGenes)."\n";
			
			############ (2) Discards flat pattern genes ##################
			
			#done with an R script
			my $outputFileWithoutFlatPatterns=$backgroundFiltered_AND_flatPatternFiltered_DIR."ALLsamples.genes.fpkm_table.background_genes_AND_flat_patterns_removed.xls";
			my $listOfDiscardedGenes=$backgroundFiltered_AND_flatPatternFiltered_DIR."flat_pattern_removed_genes.txt";
			system("rm -f ".$outputFileWithoutFlatPatterns);
			system("rm -f ".$listOfDiscardedGenes);
			system("head -n 1 ".$newTable_background_genes_OUT." > ".$outputFileWithoutFlatPatterns);
			my $fiteringRScritpt="countTable=read.table(\"".$newTable_background_genes_OUT."\",header=TRUE,row.names=1)\n";
			$fiteringRScritpt.="nRows=dim(countTable)[1]\n";
			$fiteringRScritpt.="nColumns=dim(countTable)[2]\n";
			$fiteringRScritpt.="for(i in 1:nRows){\n";
			$fiteringRScritpt.="\tQ1<-quantile(as.numeric(countTable[i,]))[2] #25% quantile\n";
			$fiteringRScritpt.="\tmedian<-quantile(as.numeric(countTable[i,]))[3] #50% quantile\n";
			$fiteringRScritpt.="\tQ3<-quantile(as.numeric(countTable[i,]))[4] #75% quantile\n";
			$fiteringRScritpt.="\tIQR=Q3-Q1\n";			
			$fiteringRScritpt.="\tif(IQR>".$IQRthresdhold_forFlatPatternFiltering."){write.table(countTable[i,],file=\"".$outputFileWithoutFlatPatterns."\",append=TRUE,sep=\"\t\",row.names=TRUE,col.names=FALSE)}\n";
			$fiteringRScritpt.="\telse{write.table(row.names(countTable[i,]),file=\"".$listOfDiscardedGenes."\",append=TRUE,sep=\"\t\",row.names=FALSE,col.names=FALSE,quote=FALSE)}\n";
			$fiteringRScritpt.="}\n";
			#guarda script para filtrar outliers
			my $RfilteredScript=$backgroundFiltered_AND_flatPatternFiltered_DIR."flatPatterns_filter.R";
			open(RFILTERED,">",$RfilteredScript);
			print RFILTERED $fiteringRScritpt;
			close(RFILTERED);
			
			#ejecuto R
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [Running R for] ".$RfilteredScript."\n";
			print STDOUT "[Running R for] ".$RfilteredScript."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [Running R for] ".$RfilteredScript."\n";
			
			my $command="R --vanilla < ".$RfilteredScript." > ".$RfilteredScript.".Rlog 2>&1";
			system($command);
		
			#counts the number of discarded genes
			my $nDiscarded=0;
			if(-e $listOfDiscardedGenes){
				open(NUMBER,"wc -l $listOfDiscardedGenes | cut -d' ' -f 1 |");
				$nDiscarded=<NUMBER>;
				close(NUMBER);
				chomp($nDiscarded);
			}
			
			#counts the number of genes lost due to flat pattern filtering
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to flat pattern filtering)] ".$nDiscarded."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to flat pattern filtering)] ".$nDiscarded."\n";
			print STDOUT "[number of genes discarded due to flat pattern filtering)] ".$nDiscarded."\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [total number of discarded genes] ".(($nTotalGenes-$nSelectedGenes)+$nDiscarded)."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [total number of discarded genes] ".(($nTotalGenes-$nSelectedGenes)+$nDiscarded)."\n";
			print STDOUT "[total number of discarded genes] ".(($nTotalGenes-$nSelectedGenes)+$nDiscarded)."\n";			
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of remaining/selected genes] ".($nTotalGenes-(($nTotalGenes-$nSelectedGenes)+$nDiscarded))."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of remaining/selected genes] ".($nTotalGenes-(($nTotalGenes-$nSelectedGenes)+$nDiscarded))."\n";
			print STDOUT "[number of remaining/selected genes] ".($nTotalGenes-(($nTotalGenes-$nSelectedGenes)+$nDiscarded))."\n";
	
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [creating new GTF] ".$newGTF_without_background_AND_flatpattern_genes."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [creating new GTF] ".$newGTF_without_background_AND_flatpattern_genes."\n";
			print STDOUT "[creating new GTF] ".$newGTF_without_background_AND_flatpattern_genes."\n\n";
		
			#generates the new GTF file that excludes background level genes + flat pattern genes
			if(open(REDUCEDLIST,"<",$outputFileWithoutFlatPatterns)){			
				foreach my $line(<REDUCEDLIST>){
					chomp($line);
				
					if($line!~ /^Name/){
						my @tokens=split('\t',$line);					
						my $command="grep 'gene_id ".$tokens[0].";' ".$initialGTF." >> ".$newGTF_without_background_AND_flatpattern_genes;
						system($command);
					}
				}
			}			
			
		}#if(open(IN,"<",$inputFileForBackgroundFiltering))
		else{
			print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch function]: ".$inputFileForBackgroundFiltering." does not exist\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch function]: ".$inputFileForBackgroundFiltering." does not exist\n";
			print STDERR "\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch function]: ".$inputFileForBackgroundFiltering." does not exist\n\n";
		exit(-1);
		}
			
	}#if(-d $cuffnormOutDir && -d $backgroundFiltered_AND_flatPatternFiltered_DIR)
	else{
		print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch function]: ".$cuffnormOutDir." or ".$backgroundFiltered_AND_flatPatternFiltered_DIR." do not exist\n";
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch function]: ".$cuffnormOutDir." or ".$backgroundFiltered_AND_flatPatternFiltered_DIR." do not exist\n";
		print STDERR "\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch function]: ".$cuffnormOutDir." or ".$backgroundFiltered_AND_flatPatternFiltered_DIR." do not exist\n\n";
		exit(-1);		
	}
		
	close(LOGFILEOUT);

}

sub removeBackgroundLevelGenesANDFlatPatternGenes_for_deseq_branch{
#After running the initial deseq for ALLsamples, it removes from 'ALLsamples.normalizedCounts.xls' those genes with expression values at background levels + flat pattern genes.
#Finally, it creates a new annotation GTF, without the removed genes, and runs DESeq2 again
	
	##if using $initialGTF => spikes annotation is not considered (modified initial GTF with spikes info is not used here)
	my ($deseqParams,$initialGTF,$workspace,$logfh)=@_;
	

	# <----------------- only for testing ----------------------->	 
#	use Data::Dumper;
#	print STDERR "FILE:\n".Dumper($comparisons)."\n";
	

	# <----------------- DESeq2 ----------------------->
	my $deseqOutDir=$workspace."deseq/";
	my $backgroundFiltered_AND_flatPatternFiltered_DIR=$workspace."deseq_backgroundFiltered_AND_flatPatternFiltered/";
	File::Path::make_path($backgroundFiltered_AND_flatPatternFiltered_DIR);

	my $filteringLogFile=$backgroundFiltered_AND_flatPatternFiltered_DIR."filtering.log.txt";
	
	open(LOGFILEOUT,">>",$filteringLogFile);
	print STDOUT "[removing background level genes and flat pattern genes for DESeq2]\n";
	print $logfh (Miscellaneous->getCurrentDateAndTime())."[removing background level genes and flat pattern genes for DESeq2]\n";	
	print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[removing background level genes and flat pattern genes for DESeq2]\n";
	
	if(-d $deseqOutDir && -d $backgroundFiltered_AND_flatPatternFiltered_DIR){
		my $backgroundLevel=$deseqParams->[0]->{backgroundLevel};
		my $required_PercentageOfSamplesOverBackgroundLevel=$deseqParams->[0]->{requiredPercentageOfSamplesOverBackgroundLevel};
		my $IQRthresdhold_forFlatPatternFiltering=$deseqParams->[0]->{IQRthresdhold_forFlatPatternFiltering};
	
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[background expression threshold] ".$backgroundLevel."\n";
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[required percentage of samples for background expression threshold] ".$required_PercentageOfSamplesOverBackgroundLevel."\n";
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."[IQR threshold for flat pattern filtering] ".$IQRthresdhold_forFlatPatternFiltering."\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[background expression threshold] ".$backgroundLevel."\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[required percentage of samples for background expression threshold] ".$required_PercentageOfSamplesOverBackgroundLevel."\n";
		print $logfh (Miscellaneous->getCurrentDateAndTime())."[IQR threshold for flat pattern filtering] ".$IQRthresdhold_forFlatPatternFiltering."\n";
		print STDOUT "[background expression threshold] ".$backgroundLevel."\n";
		print STDOUT "[required percentage of samples for background expression threshold] ".$required_PercentageOfSamplesOverBackgroundLevel."\n";
		print STDOUT "[IQR threshold for flat pattern filtering] ".$IQRthresdhold_forFlatPatternFiltering."\n";	
		
		my $inputFileForBackgroundFiltering=$deseqOutDir."ALLsamples.normalizedCounts.xls";
		my $newTable_background_genes_OUT=$backgroundFiltered_AND_flatPatternFiltered_DIR."ALLsamples.normalizedCounts_table.background_genes_removed.xls";
		my $newGTF_without_background_AND_flatpattern_genes=$backgroundFiltered_AND_flatPatternFiltered_DIR."GTF_without_background_AND_flatpatternGenes.gtf";
		my $list_with_discardedBackgroundGenes=$backgroundFiltered_AND_flatPatternFiltered_DIR."background_level_removed_genes.txt";
		my $nTotalGenes=0;
		my $nSelectedGenes=0;
		
		system("rm -f ".$newTable_background_genes_OUT);
		system("rm -f ".$newGTF_without_background_AND_flatpattern_genes);
		system("rm -f ".$list_with_discardedBackgroundGenes);
		
		if(open(IN,"<",$inputFileForBackgroundFiltering)){
			open(OUT,">",$newTable_background_genes_OUT);
			open(DISCARDED,">",$list_with_discardedBackgroundGenes);
			
			my $nSamples=0;
			
			############ (1) Discards genes under background expression level ##################
			
			foreach my $line(<IN>){
				chomp($line);
				
				#prints out header line				
				if($line=~ /^Name/){
					print OUT $line."\n";
					
					#counts the number of total samples
					my @tokens=split('\t',$line);	
					
					for (my $j=1;$j<@tokens;$j++){
						$nSamples++;						
					}
					
					print $logfh (Miscellaneous->getCurrentDateAndTime())." [Number of total samples] ".$nSamples."\n";
					print STDOUT "[Number of total samples] ".$nSamples."\n";
					print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [Number of total samples] ".$nSamples."\n";
						
				}else{ #lines with data values
					$nTotalGenes++;
					my $nSamplesOverBackgroundLevel=0;								
					my @tokens=split('\t',$line);
					my $gene=$tokens[0];
										
					#reads expression values in all the samples for the current gene					
					for (my $j=1;$j<@tokens;$j++){
						my $value=$tokens[$j];
						
						if($value>$backgroundLevel){
							$nSamplesOverBackgroundLevel++;
						}						
					}					
					
					my $percentageOfSamplesOverBackgroundLevel=($nSamplesOverBackgroundLevel/$nSamples)*100;					
					
					# if the current gene is over background expression level for the required number of samples
					if($percentageOfSamplesOverBackgroundLevel>=$required_PercentageOfSamplesOverBackgroundLevel){
						#gene is selected as it has passed background expression level filter
						print OUT $line."\n";
						$nSelectedGenes++;
					}
					else{
						print DISCARDED $gene."\n";
					}			
				}			
			}
			
			close(OUT);
			close(DISCARDED);
			close(IN);

			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of total genes] ".$nTotalGenes."\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to background level filtering] ".($nTotalGenes-$nSelectedGenes)."\n";
			print STDOUT "[number of total genes] ".$nTotalGenes."\n";
			print STDOUT "[number of genes discarded due to background level filtering] ".($nTotalGenes-$nSelectedGenes)."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of total genes] ".$nTotalGenes."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to background level filtering] ".($nTotalGenes-$nSelectedGenes)."\n";
			
			############ (2) Discards flat pattern genes ##################
			
			#done with an R script
			my $outputFileWithoutFlatPatterns=$backgroundFiltered_AND_flatPatternFiltered_DIR."ALLsamples.normalizedCounts_table.background_genes_AND_flat_patterns_removed.xls";
			my $listOfDiscardedGenes=$backgroundFiltered_AND_flatPatternFiltered_DIR."flat_pattern_removed_genes.txt";
			system("rm -f ".$outputFileWithoutFlatPatterns);
			system("rm -f ".$listOfDiscardedGenes);
			system("head -n 1 ".$newTable_background_genes_OUT." > ".$outputFileWithoutFlatPatterns);
			my $fiteringRScritpt="countTable=read.table(\"".$newTable_background_genes_OUT."\",header=TRUE,row.names=1)\n";
			$fiteringRScritpt.="nRows=dim(countTable)[1]\n";
			$fiteringRScritpt.="nColumns=dim(countTable)[2]\n";
			$fiteringRScritpt.="for(i in 1:nRows){\n";
			$fiteringRScritpt.="\tQ1<-quantile(as.numeric(countTable[i,]))[2] #25% quantile\n";
			$fiteringRScritpt.="\tmedian<-quantile(as.numeric(countTable[i,]))[3] #50% quantile\n";
			$fiteringRScritpt.="\tQ3<-quantile(as.numeric(countTable[i,]))[4] #75% quantile\n";
			$fiteringRScritpt.="\tIQR=Q3-Q1\n";			
			$fiteringRScritpt.="\tif(IQR>".$IQRthresdhold_forFlatPatternFiltering."){write.table(countTable[i,],file=\"".$outputFileWithoutFlatPatterns."\",append=TRUE,sep=\"\t\",row.names=TRUE,col.names=FALSE)}\n";
			$fiteringRScritpt.="\telse{write.table(row.names(countTable[i,]),file=\"".$listOfDiscardedGenes."\",append=TRUE,sep=\"\t\",row.names=FALSE,col.names=FALSE,quote=FALSE)}\n";
			$fiteringRScritpt.="}\n";
			#guarda script para filtrar outliers
			my $RfilteredScript=$backgroundFiltered_AND_flatPatternFiltered_DIR."flatPatterns_filter.R";
			open(RFILTERED,">",$RfilteredScript);
			print RFILTERED $fiteringRScritpt;
			close(RFILTERED);
			
			#ejecuto R
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [Running R for] ".$RfilteredScript."\n";
			print STDOUT "[Running R for] ".$RfilteredScript."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [Running R for] ".$RfilteredScript."\n";
			
			my $command="R --vanilla < ".$RfilteredScript." > ".$RfilteredScript.".Rlog 2>&1";
			system($command);
		
			#counts the number of discarded genes
			my $nDiscarded=0;
			if(-e $listOfDiscardedGenes){
				open(NUMBER,"wc -l $listOfDiscardedGenes | cut -d' ' -f 1 |");
				$nDiscarded=<NUMBER>;
				close(NUMBER);
				chomp($nDiscarded);
			}
			
			#counts the number of genes lost due to flat pattern filtering
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to flat pattern filtering)] ".$nDiscarded."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of genes discarded due to flat pattern filtering)] ".$nDiscarded."\n";
			print STDOUT "[number of genes discarded due to flat pattern filtering)] ".$nDiscarded."\n";
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [total number of discarded genes] ".(($nTotalGenes-$nSelectedGenes)+$nDiscarded)."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [total number of discarded genes] ".(($nTotalGenes-$nSelectedGenes)+$nDiscarded)."\n";
			print STDOUT "[total number of discarded genes] ".(($nTotalGenes-$nSelectedGenes)+$nDiscarded)."\n";			
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [number of remaining/selected genes] ".($nTotalGenes-(($nTotalGenes-$nSelectedGenes)+$nDiscarded))."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [number of remaining/selected genes] ".($nTotalGenes-(($nTotalGenes-$nSelectedGenes)+$nDiscarded))."\n";
			print STDOUT "[number of remaining/selected genes] ".($nTotalGenes-(($nTotalGenes-$nSelectedGenes)+$nDiscarded))."\n";
	
			print $logfh (Miscellaneous->getCurrentDateAndTime())." [creating new GTF] ".$newGTF_without_background_AND_flatpattern_genes."\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())." [creating new GTF] ".$newGTF_without_background_AND_flatpattern_genes."\n";
			print STDOUT "[creating new GTF] ".$newGTF_without_background_AND_flatpattern_genes."\n\n";
		
			#generates the new GTF file that excludes background level genes + flat pattern genes
			if(open(REDUCEDLIST,"<",$outputFileWithoutFlatPatterns)){			
				foreach my $line(<REDUCEDLIST>){
					chomp($line);
				
					if($line!~ /^Name/){
						my @tokens=split('\t',$line);					
						my $command="grep 'gene_id ".$tokens[0].";' ".$initialGTF." >> ".$newGTF_without_background_AND_flatpattern_genes;
						system($command);
					}
				}
			}			
			
		}#if(open(IN,"<",$inputFileForBackgroundFiltering))
		else{
			print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_deseq_branch function]: ".$inputFileForBackgroundFiltering." does not exist\n";
			print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_deseq_branch function]: ".$inputFileForBackgroundFiltering." does not exist\n";
			print STDERR "\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_deseq_branch function]: ".$inputFileForBackgroundFiltering." does not exist\n\n";
		exit(-1);
		}
			
	}#if(-d $deseqOutDir && -d $backgroundFiltered_AND_flatPatternFiltered_DIR)
	else{
		print $logfh (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_deseq_branch function]: ".$deseqOutDir." or ".$backgroundFiltered_AND_flatPatternFiltered_DIR." do not exist\n";
		print LOGFILEOUT (Miscellaneous->getCurrentDateAndTime())."\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_deseq_branch function]: ".$deseqOutDir." or ".$backgroundFiltered_AND_flatPatternFiltered_DIR." do not exist\n";
		print STDERR "\n[ERROR in removeBackgroundLevelGenesANDFlatPatternGenes_for_deseq_branch function]: ".$deseqOutDir." or ".$backgroundFiltered_AND_flatPatternFiltered_DIR." do not exist\n\n";
		exit(-1);		
	}
		
	close(LOGFILEOUT);

}


1;
