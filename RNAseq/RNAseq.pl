#!/usr/bin/perl -w

# nextpresso

# RNAseq.pl 
# Author: Osvaldo Grana
# Description : RNA-seq analysis pipeline
# v1.9.1	oct2017
#
# v1.9.2	ene2018 - removes genes with expression levels below background + removes genes with flat pattern expression.
#			  After this removal, a new GTF is created that contains only those genes that passed the filtering: cuffquant+cuffdiff+cuffnorm
#			  or htseqcount+deseq is run again using this reduced gene annotation (new GTF)
#			  (see addings in level 5 and level 6 here, plus the addings in ExecutionLevels.pm


my $version="v1.9.2, ene2018";

#new in v1.9.1: creates the temporal directory (tmp) inside the workspace, and not inside the /tmp of the machine

use strict;
use warnings;
use autodie;

#determines current working directory
#use Cwd qw(cwd); 
#my $bin=cwd;


use File::Temp qw(tempdir); #assigns a temp. subdir inside the main temp. dir

use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function

#loads own packages
use lib catdir(dirname($Bin), 'Utils'); #finds out Utils dir from parent dir
use xmlParRNAseq 'processXML'; #loads modules from Utils
use ExecutionLevels;
use Miscellaneous;

use Getopt::Long; #to get options
#use File::Basename;
use File::Spec;
use File::Path qw(make_path);
use FileHandle;

#use File::Copy;
#use Sys::Hostname;

use Carp qw( confess ); # to verbose stack traces
use Config; # to check if perl was compiled with thread support


#subroutine prototypes
sub main();
sub help();
sub checkProgramPaths($);
sub checkSampleFiles($$);

main(); #calls main function


sub main(){
	
	system("clear");
	print "*****************************************************************************\n";
	print "*                                                                           *\n";
	print "*    nextpresso: next generation sequencing expression analysis pipeline    *\n";
	print "*    ".$version."                                                        *\n";
	print "*    Author: Osvaldo Grana                                                  *\n";
	print "*                                                                           *\n";
	print "*****************************************************************************\n\n";
	
	#creates a temporal subdirectory in the temporal directory of the machine (like '/tmp')
	#with an unrecognized name like ('/tmp/mY0qHO36dP')
	# CLEANUP => 1 implies that this subdirectory is removed when the execution finishes
	# CLEANUP => 0 implies that this subdirectory is NOT removed when the execution finishes
	#my $executionCreatedTempDir = tempdir( CLEANUP => 0 );
	

	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running nextpresso.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;
	
	
	my $level = "1345"; # default value;
	
	my $configXMLSchema  = $Bin."config/config.xsd";
	my $experimentXMLSchema = $Bin."config/experiment.xsd";
	my $configXMLDocument = undef;	
	my $experimentXMLDocument = undef;

	GetOptions(
		"configDoc=s"=>\$configXMLDocument,	#string
		"expDoc=s"=>\$experimentXMLDocument,	#string				
		"step=i"=>\$level 		#numeric
	);

	#it always checks level 0
	$level="0".$level;
	
	if(defined($configXMLDocument) && defined($configXMLDocument)){
		if((!-e $configXMLSchema) || (!-e $configXMLDocument) || (!-e $experimentXMLSchema) || (!-e $experimentXMLDocument)){
			print "**** One or several of the following files do not exist:\n\n";
			print "File 1: ".$configXMLSchema."\nFile 2: ".$configXMLDocument."\nFile 3: ".$experimentXMLSchema."\nFile 4: ".$experimentXMLDocument."\n";
			print "\n[Execution finished]\n\n";
		
			help();
		}
	}else{
		help();
	}


	#<------------- GETS THE DATA FROM THE XML DOCUMENTS ------------->
	#validates XML documents and returns them as hash tables
	my $configHashRef=xmlParRNAseq::processXML($configXMLSchema, $configXMLDocument);
	my $experimentHashRef=xmlParRNAseq::processXML($experimentXMLSchema, $experimentXMLDocument);
	
	#returns a new hash ref, with "/" added at the end of path lines in the case it was not present
	$configHashRef=checkProgramPaths($configHashRef);
	
	
	#getting config info
	my $fastQCpath=${$configHashRef}{fastQCpath}[0];
	my $fastQScreenPath=${$configHashRef}{fastQScreen}[0]{path}[0];
	my $fastQScreenConf=${$configHashRef}{fastQScreen}[0]{configurationFile}[0];
	my $fastQScreenSubset=${$configHashRef}{fastQScreen}[0]{subset}[0];
	my $bedtoolsPath=${$configHashRef}{bedtoolsPath}[0];
	my $samtoolsPath=${$configHashRef}{samtoolsPath}[0];
	my $bowtiePath=${$configHashRef}{bowtiePath}[0];
	my $tophatPath=${$configHashRef}{tophatPath}[0];
	my $peakAnnotatorPath=${$configHashRef}{peakAnnotatorPath}[0];
	my $htseqCountPath=${$configHashRef}{htseqCount}[0]{path}[0];
	my $htseqCountPythonpath=${$configHashRef}{htseqCount}[0]{pythonpath}[0];
	my $tophatFusion=${$configHashRef}{tophatFusion}[0]{path}[0];
	my $cufflinksPath=${$configHashRef}{cufflinks}[0]{path}[0];
	my $bedGraphToBigWigPath=${$configHashRef}{bedGraphToBigWig}[0]{path}[0];
	my $maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep=${$configHashRef}{maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep};
	my $seqtkPath=${$configHashRef}{seqtk}[0]{path}[0];
	my $seqtk_maximunNumberOfInstancesForDownSampling=${$configHashRef}{seqtk}[0]{maximunNumberOfInstancesForDownSampling}[0];
	my $queueSystem=${$configHashRef}{queueSystem}[0];
	my $queueName=${$configHashRef}{queueName}[0];	
	my $multicore=${$configHashRef}{multicore}[0];
	my $perl5lib=${$configHashRef}{PERL5LIB}[0];
	my $gseaPath=${$configHashRef}{gsea}[0]{path}[0];
	my $gseaChip=${$configHashRef}{gsea}[0]{chip}[0];
	my $gseamaxMemory=${$configHashRef}{gsea}[0]{maxMemory}[0];
	my $extraPathsRequired=${$configHashRef}{extraPathsRequired}[0];
	if($extraPathsRequired=~ /HASH\(/){$extraPathsRequired="NO_EXTRA_PATHS"}
	
	#getting experiment info	
	my $experimentName=${$experimentHashRef}{projectName};
	my $queueProject=$experimentName;	
	
	
	my $workspace=${$experimentHashRef}{workspace};
	#creates a temporal directory
	my $executionCreatedTempDir=$workspace."/tmp/".$experimentName;
	my $initialTimeANDdate=Miscellaneous::getCurrentDateAndTime();
	$initialTimeANDdate=~ s/ //g;
	my $temporaryFilesPrefix=$initialTimeANDdate;
     	$temporaryFilesPrefix=~ s/\[//;
     	$temporaryFilesPrefix=~ s/\]//;
     	$temporaryFilesPrefix=~ s/-//g;
     	$temporaryFilesPrefix=~ s/,/_/;
     	$temporaryFilesPrefix=~ s/\:/-/g;   
	$executionCreatedTempDir.="_".$temporaryFilesPrefix;
	
	#removes it and create it again
	File::Path::make_path($executionCreatedTempDir,$executionCreatedTempDir);
	

	
	
	
	my $referenceSequence=${$experimentHashRef}{referenceSequence};
	my $indexPrefixForReferenceSequence=${$experimentHashRef}{referenceSequence};
	$indexPrefixForReferenceSequence=~ s/\.fa$//;
	$indexPrefixForReferenceSequence=~ s/\.fasta$//;
	my $GTF=${$experimentHashRef}{GTF};
	# one copy of the original GTF is preserved below in mind (initialGTF). This is valid when having spikes in the created annotation.
	# It is more appropriated, in this case, to run htseqcount and cuffdiff withou including spikes in the GTF annotation,
	# as their quantification values could affect and modify FPKM normalization for the rest of the genes.
	my $initialGTF=$GTF;
	my $samples=${$experimentHashRef}{library}; # hashRef
	my $comparisons=${$experimentHashRef}{comparison}; # hashRef
	my $tophatParams=${$experimentHashRef}{tophat}; # hashRef
	my $pairedEnd=${$experimentHashRef}{pairedEnd};
	my $fileWithChecksumCodesToValidate=${$experimentHashRef}{fileWithChecksumCodesToValidate};
	my $cufflinksParams=${$experimentHashRef}{cufflinks}; # hashRef
	my $cuffmergeParams=${$experimentHashRef}{cuffmerge}; # hashRef
	my $cuffquantParams=${$experimentHashRef}{cuffquant}; # hashRef	
	my $cuffnormParams=${$experimentHashRef}{cuffnorm}; # hashRef
	my $cuffdiffParams=${$experimentHashRef}{cuffdiff}; # hashRef
	my $htseqcountParams=${$experimentHashRef}{htseqcount}; # hashRef
	my $deseqParams=${$experimentHashRef}{deseq2};
	my $bedGraphToBigWigParams=${$experimentHashRef}{bedGraphToBigWig};
	my $gseaParams=${$experimentHashRef}{gsea};
	my $tophatfusionParams=${$experimentHashRef}{tophatfusion};
	my $spikeInControlMixesParams=${$experimentHashRef}{spikeInControlMixes};
	my $doSpikesAndGenomeRefIndexing=lc($spikeInControlMixesParams->[0]->{do});
	
	#checks that the sample files exist
	checkSampleFiles($pairedEnd,$samples);
	
	
	# <----------------- only for testing ----------------------->	 	
	#use Data::Dumper;
	#print STDERR "FILE:\n".Dumper($experimentHashRef)."\n";



	#<------------- ANALYSIS OF THE DATA STARTS ------------->

	
	
	if (! File::Spec->file_name_is_absolute($workspace)){
		print STDERR "\n[ERROR]: Incorrect absolute path for output directory: $workspace\n\n";
		exit(-1);
	}
	
	# the workspace directory is created
	if($workspace!~ /\/$/){
		$workspace.="/";
	}
	File::Path::make_path($workspace);
	
	if(!-d $workspace){
		print STDERR "\n[ERROR]: Cannot create workspace $workspace\n\n";
		exit(-1);
	}
	
	my $logFile=$workspace."run.log";	
	my $logfh=FileHandle->new(">".$logFile);
	if(!-e $logfh){
		print STDERR "\n[ERROR]: Cannot create ".$logfh."\n\n";
		exit(-1);
	}else{
		print $logfh $initialTimeANDdate." RNAseq pipeline: Starting analysis for $experimentName experiment\n";
		if(-d $executionCreatedTempDir){
			print $logfh (Miscellaneous::getCurrentDateAndTime())."[DONE]: created temporal subdirectory for this execution ".$executionCreatedTempDir."\n";}
		else{
			print $logfh (Miscellaneous::getCurrentDateAndTime())."[ERROR]: could not create temporal subdirectory for this execution ".$executionCreatedTempDir."\n";
			print STDERR "\n[ERROR]: could not create temporal subdirectory for this execution ".$executionCreatedTempDir."\n";
			exit(-1);
		}
	}

	########## level 0: converts raw read bam files to fastq files if needed or/and prepares reference and GTF index files in case of having spike-in control mixes	
	if($level=~ /0/){
		#**** IMPORTANT: since this level is always done, steps within this level are done ONLY if they were not done before.
		#**** for example: suppose that the user wants to do differential expression again (adding a new comparison), in this case it wouldn't make sense
		#**** to perform again bam files conversion to fastq files or indexing again the reference for spike-in controls, as these steps were already done
		#**** with the first execution of the analysis
		
		#**** IMPORTANT: in the case of using spike-in control mixes, the reference file and the GTF file are different and new files,
		#**** that's why they are recovered here
		($referenceSequence,$GTF,$indexPrefixForReferenceSequence)=ExecutionLevels::level_0($fileWithChecksumCodesToValidate,$workspace,$experimentName,$logfh,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$spikeInControlMixesParams,$referenceSequence,$GTF,$samples,$bedtoolsPath,$pairedEnd,$bowtiePath,$indexPrefixForReferenceSequence,
			$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);
	}

	########## level 1: sequencing quality and contamination check ##########		
	if($level=~ /1/){
		ExecutionLevels::level_1($perl5lib,$fastQCpath,$fastQScreenPath,$fastQScreenConf,$bowtiePath,$experimentName,$workspace,$referenceSequence,
		$GTF,$samples,$logfh,$executionCreatedTempDir,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$fastQScreenSubset,$pairedEnd,$queueSystem,$queueName,$multicore,$queueProject);
	}
	
	########## level 2: trimming && downsampling ##########		
	if($level=~ /2/){
		ExecutionLevels::level_2($fastQCpath,$fastQScreenPath,$fastQScreenConf,$bowtiePath,$experimentName,$workspace,$referenceSequence,
		$GTF,$samples,$logfh,$executionCreatedTempDir,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$fastQScreenSubset,$pairedEnd,$seqtkPath,
		$seqtk_maximunNumberOfInstancesForDownSampling,$queueSystem,$queueName,$multicore,$queueProject);

		#for those samples that were trimmed, or trimmed & downsampled (but not only downsampled)
		#it must perform level1 again
		my %auxHash=%$samples; # hash unref	
		my $auxiliarSamples=\%auxHash; # hash ref again

		foreach my $key (keys %$auxiliarSamples){
			my $trimming=$auxiliarSamples->{$key}{trimming}[0]{do};
			if($trimming eq "false"){
				delete $auxiliarSamples->{$key};
			}
		}
		
		#if the size of the hash %auxiliarSamples==0, i.e., none of the samples required trimming => no additional FASTQC is required
		if(keys(%$auxiliarSamples)>0){
			ExecutionLevels::level_1($perl5lib,$fastQCpath,$fastQScreenPath,$fastQScreenConf,$bowtiePath,$experimentName,$workspace,$referenceSequence,
			$GTF,$auxiliarSamples,$logfh,$executionCreatedTempDir,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$fastQScreenSubset,$pairedEnd,$queueSystem,$queueName,$multicore,$queueProject);
		}

	}else{	
		#in case that a previous analysis was done for this experiment, requiring trimming/downsampling for maybe some of the samples,
		#and if a new re-analysis is started just after this step (level 3 for example), it implies that the program has to be aware of
		#what are the proper trimmed/downsampled samples
		ExecutionLevels::changeSampleNamesInCaseTheyWereTrimmedAndOrDownsampledBefore($samples,$workspace,$pairedEnd);	
	}

	########## level 3: aligning of reads ##########		
	if($level=~ /3/){
		ExecutionLevels::level_3($tophatPath,$bowtiePath,$samtoolsPath,$bedtoolsPath,$peakAnnotatorPath,$referenceSequence,
		$indexPrefixForReferenceSequence,$samples,$GTF,$tophatParams,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,$experimentName,$logfh,
		$executionCreatedTempDir,$pairedEnd,$queueSystem,$queueName,$multicore,$queueProject);
	}

	########## level 4: transcripts assembly and quantification (cufflinks and cuffmerge) ##########
	if($level=~ /4/){
		ExecutionLevels::level_4($extraPathsRequired,$spikeInControlMixesParams,$cufflinksPath,$samtoolsPath,
		$bedtoolsPath,$referenceSequence,$indexPrefixForReferenceSequence,$samples,$GTF,$cufflinksParams,
		$cuffmergeParams,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,
		$workspace,$experimentName,$logfh,$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);		
	}			
		
	########## level 5: differential expression (cuffquant, cuffdiff and cuffnorm) ##########
	if($level=~ /5/){
		ExecutionLevels::level_5($doSpikesAndGenomeRefIndexing,$initialGTF,$extraPathsRequired,$comparisons,$cufflinksPath,$samtoolsPath,$bedtoolsPath,$referenceSequence,$indexPrefixForReferenceSequence,
		$samples,$GTF,$cuffquantParams,$cuffnormParams,$cuffdiffParams,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,$experimentName,$logfh,
		$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);
		
		#removes back ground level genes + flat pattern genes
		ExecutionLevels::removeBackgroundLevelGenesANDFlatPatternGenes_for_cuffdiff_branch($cuffnormParams,$initialGTF,$workspace,$logfh);
		
		#executes all again using the reduced GTF (without back ground level genes + without flat pattern genes)
		my $originalAlignmentsDir=$workspace."alignments/";		
		my $new_workspace=$workspace."cuffdiff_backgroundFiltered_AND_flatPatternFiltered/";

		# creates a symbolic link to the alignments dir, to emulate its presence in the new workspace directory		
		my $reducedGTF=$new_workspace."GTF_without_background_AND_flatpatternGenes.gtf";
		my $command="ln -s ".$originalAlignmentsDir." ".$new_workspace;
		system($command);

		ExecutionLevels::level_5($doSpikesAndGenomeRefIndexing,$reducedGTF,$extraPathsRequired,$comparisons,$cufflinksPath,$samtoolsPath,$bedtoolsPath,$referenceSequence,$indexPrefixForReferenceSequence,
		$samples,$reducedGTF,$cuffquantParams,$cuffnormParams,$cuffdiffParams,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$new_workspace,$experimentName,$logfh,
		$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);		
		
	}
	
	########## level 6: runs htseq-count (gets read counts for genes) + DESeq2 differential expression
	if($level=~ /6/){
		if($doSpikesAndGenomeRefIndexing eq "false"){		
			ExecutionLevels::level_6($perl5lib,$comparisons,$deseqParams,$extraPathsRequired,$htseqCountPath,$htseqCountPythonpath,$htseqcountParams,$samtoolsPath,$samples,$GTF,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,$experimentName,$logfh,
			$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);
		}else{ # when having spikes, it is more appropriate to not consider them for htseqcount as they could affect
			# normalization values for regular genes. So in this case, the original GTF is given instead of the
			#one with the combined annotation (genes+spikes) 
			ExecutionLevels::level_6($perl5lib,$comparisons,$deseqParams,$extraPathsRequired,$htseqCountPath,$htseqcountParams,$samtoolsPath,$samples,$initialGTF,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,$experimentName,$logfh,
			$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);		
		}
	}
	
	########## level 7: creates wiggle files from bam alignments
	if($level=~ /7/){
		ExecutionLevels::level_7($bedGraphToBigWigPath,$bedGraphToBigWigParams,$bedtoolsPath,$samtoolsPath,$samples,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$workspace,$experimentName,$logfh,
		$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);		
	}
	
	########## level 8: preRanked GSEA
	if($level=~ /8/){
		#****rnk files are directly taken from cuffdiff output files
		ExecutionLevels::level_8($gseaChip,$gseamaxMemory,$gseaPath,$gseaParams,$workspace,$experimentName,$logfh,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$comparisons,
			$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);		
	}
	
	########## level 9: gene fusion prediction with Tophat-fusion
	if($level=~ /9/){
		ExecutionLevels::level_9($tophatfusionParams,$tophatPath,$bowtiePath,$samtoolsPath,$workspace,$experimentName,$logfh,$maximunNumberOfInstancesAllowedToRunSimultaneouslyInOneParticularStep,$samples,$referenceSequence,
			$executionCreatedTempDir,$queueSystem,$queueName,$multicore,$queueProject);		
	}	
	
	#close the log file
	$logfh->close;
	
}

sub help(){
	my $usage = qq{
		
		perl RNAseq.pl --configDoc configDocFile --expDoc expDocFile --step step_number
		

			Example:
			
				a) complete execution of all steps in each workflow level
				perl RNAseq.pl --configDoc config/configurationParameters.xml --expDoc config/experimentParameters.xml --step 123456789 
				
				b) execution of some detailed steps
				perl RNAseq.pl --configDoc config/configurationParameters.xml --expDoc config/experimentParameters.xml --step 1345


			Steps Description:
			
				Step 1: sequencing quality && contamination check (fastQC & fastQScreen)
				Step 2: trimming && downsampling (seqtk)
				Step 3: Aligning (tophat)
				Step 4: transcripts assembly && quantification (cufflinks and cuffmerge)
				Step 5: differential expression (cuffquant, cuffdiff and cuffnorm)
				Step 6: htseq-count (gets read counts for genes) + DESeq2 differential expression
				Step 7: BedGraph and BigWig files for genome browsers
				Step 8: GSEA for specific gene sets over the different comparisons done with cuffdiff
				Step 9: gene fusion prediction
				
				
};
 
	print STDERR $usage;
	exit(1);
}

sub checkProgramPaths($){
	my ($configHashRef) = @_;
	
	if(${$configHashRef}{fastQCpath}[0]!~ /\/$/){
		${$configHashRef}{fastQCpath}[0]=${$configHashRef}{fastQCpath}[0].="/";
	}
	my $program=${$configHashRef}{fastQCpath}[0]."fastqc";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{fastQScreen}[0]{path}[0]!~ /\/$/){
		${$configHashRef}{fastQScreen}[0]{path}[0]=${$configHashRef}{fastQScreen}[0]{path}[0].="/";
	}
	$program=${$configHashRef}{fastQScreen}[0]{path}[0]."fastq_screen";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}
	
	$program=${$configHashRef}{fastQScreen}[0]{configurationFile}[0];
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}	

	if(${$configHashRef}{seqtk}[0]{path}[0]!~ /\/$/){
		${$configHashRef}{seqtk}[0]{path}[0]=${$configHashRef}{seqtk}[0]{path}[0].="/";
	}
	
	$program=${$configHashRef}{seqtk}[0]{path}[0]."seqtk";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{bedtoolsPath}[0]!~ /\/$/){
		${$configHashRef}{bedtoolsPath}[0]=${$configHashRef}{bedtoolsPath}[0].="/";
	}
	$program=${$configHashRef}{bedtoolsPath}[0]."bedtools";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{samtoolsPath}[0]!~ /\/$/){
		${$configHashRef}{samtoolsPath}[0]=${$configHashRef}{samtoolsPath}[0].="/";
	}
	$program=${$configHashRef}{samtoolsPath}[0]."samtools";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{bowtiePath}[0]!~ /\/$/){
		${$configHashRef}{bowtiePath}[0]=${$configHashRef}{bowtiePath}[0].="/";
	}
	$program=${$configHashRef}{bowtiePath}[0]."bowtie";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}
	
	if(${$configHashRef}{tophatPath}[0]!~ /\/$/){
		${$configHashRef}{tophatPath}[0]=${$configHashRef}{tophatPath}[0].="/";
	}
	$program=${$configHashRef}{tophatPath}[0]."tophat";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{htseqCount}[0]{path}[0]!~ /\/$/){
		${$configHashRef}{htseqCount}[0]{path}[0]=${$configHashRef}{htseqCount}[0]{path}[0].="/";
	}
	$program=${$configHashRef}{htseqCount}[0]{path}[0]."htseq-count";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{tophatFusion}[0]{path}[0]!~ /\/$/){
		${$configHashRef}{tophatFusion}[0]{path}[0]=${$configHashRef}{tophatFusion}[0]{path}[0].="/";
	}
	$program=${$configHashRef}{tophatFusion}[0]{path}[0]."tophat-fusion-post";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{cufflinks}[0]{path}[0]!~ /\/$/){
		${$configHashRef}{cufflinks}[0]{path}[0]=${$configHashRef}{cufflinks}[0]{path}[0].="/";
	}
	$program=${$configHashRef}{cufflinks}[0]{path}[0]."cufflinks";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}
	$program=${$configHashRef}{cufflinks}[0]{path}[0]."cuffmerge";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}
	$program=${$configHashRef}{cufflinks}[0]{path}[0]."cuffquant";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}$program=${$configHashRef}{cufflinks}[0]{path}[0]."cuffnorm";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}$program=${$configHashRef}{cufflinks}[0]{path}[0]."cuffdiff";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}
	
	if(${$configHashRef}{bedGraphToBigWig}[0]{path}[0]!~ /\/$/){
		${$configHashRef}{bedGraphToBigWig}[0]{path}[0]=${$configHashRef}{bedGraphToBigWig}[0]{path}[0].="/";
	}
	$program=${$configHashRef}{bedGraphToBigWig}[0]{path}[0]."bedGraphToBigWig";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}

	if(${$configHashRef}{peakAnnotatorPath}[0]!~ /\/$/){
		${$configHashRef}{peakAnnotatorPath}[0]=${$configHashRef}{peakAnnotatorPath}[0].="/";
	}
	$program=${$configHashRef}{peakAnnotatorPath}[0]."PeakAnnotator.jar";
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}
	
	$program=${$configHashRef}{gsea}[0]{path}[0];
 	if (!-e $program){
		print STDERR "\n[ERROR]: The program ".$program." doesn't exist\n\n";
		exit(-1);
	}
	return($configHashRef);	
}

sub checkSampleFiles($$){
	my($pairedEnd,$samples)=@_;
	
	foreach my $key (keys %$samples){
		my $type=$samples->{$key}{type}[0];
		my $leftFastqFile=$samples->{$key}{leftFile};

		if(!-e $leftFastqFile){
			print STDERR "\n[ERROR]: The left sample file ".$leftFastqFile." doesn't exist\n\n";
			exit(-1);		
		}
		
		#if paired end experiment
		if($pairedEnd eq "true"){
			my $rightFastqFile=$samples->{$key}{rightFile}[0];
			
			if(!-e $leftFastqFile){
				print STDERR "\n[ERROR]: The right sample file ".$rightFastqFile." doesn't exist\n\n";
				exit(-1);		
			}				
		}
	}
		
}


