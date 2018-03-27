#!/usr/bin/perl
# FileName :
# Author : Osvaldo Graña (a modified version of the one written by Miriam Rubio)
# v0.2 ene2017 - adds PERL5LIB as argument for fastqscreen

use strict;
use warnings;

package quality;

# Function: fastqcAnalysis
# Description: This function executes a fastqc analysis
# $fastqcPath: string. path to fastqc binary
# $files_in: fastq input files
# $dirOut: Output directory
# $nthr:thread number
sub fastqcAnalysis($$$$){
	my($fastqcPath,$files_in,$dirOut,$nthr)= @_;
        my $files= join(" ",@{$files_in});
	$dirOut .= "/FastQC";
	filesFunc::makeDir($dirOut);
	print "\n Fastqc Analysis for $files";
	system("time $fastqcPath/fastqc -t $nthr -o $dirOut --noextract $files");
	print "\n End Fastqc Analysis for $files";
}

# Function: singleFastqc
# Description: This function executes a fastqc analysis for one fastqc file
# $fastqcPath: string. path to fastqc binary
# $inputFile: fastq input file
# $outDir: output directory
sub singleFastqc($$$){
	my($fastqcPath,$inputFile,$outDir)= @_;
        
	#my $command="perl ".$fastqcPath."fastqc --nogroup -o ".$outDir." --noextract ".$inputFile;
	my $command="perl ".$fastqcPath."fastqc -o ".$outDir." --noextract ".$inputFile;
	print "\n\t[executing] ".$command."\n";
	system($command);	
}

# Function: fastqScreen
# Description: This function executes a fastqScreen analysis for one fastq file(single-end) or two fastq-files (paired end)
# $fastqScreenPath: string. path to fastqc binary
# $inputFile: one or two input files, depending if it is single end or paired end
# $outDir: output directory
# $fastqScreenConf: path to fastqScreen configuration file
# $subset: number of reads to test
sub fastqScreen($$$$$$$){
	my($perl5lib_nextpresso,$fastqScreenPath,$inputFile,$outDir,$fastQScreenConf,$subset,$fastqFileIlluminaQualityEncodingForFastqScreen)=@_;

	use Env qw(PERL5LIB);
	#my $command="export PERL5LIB=".$PERL5LIB."; ".$fastqScreenPath."fastq_screen ";
	
	my $command="";
	
	if(defined($PERL5LIB) && $PERL5LIB ne ""){
		$command.="export PERL5LIB=".$perl5lib_nextpresso.":".$PERL5LIB."; ";
	}else{
		$command.="export PERL5LIB=".$perl5lib_nextpresso."; ";
	}
	
	
	$command.=$fastqScreenPath."fastq_screen ";
	if($fastqFileIlluminaQualityEncodingForFastqScreen eq "illumina1_3"){
		$command.="--illumina1_3 ";
	}
	
	#if paired end experiment, $inputFile must have 'left file name+space+right file name'
	#The space in between is used to decide the way of launching it
	if($inputFile=~ /.+ .+/){
		$command.="--paired ";
	}
	
	$command.="--outdir ".$outDir." --conf ".$fastQScreenConf." --aligner bowtie2 --nohits --subset ".$subset." ".$inputFile;
	print "\n\t[executing] ".$command."\n";
	system($command);
}


1
