#!/usr/bin/perl

# Author : Osvaldo Grana
# Description: runs htseqCount to count reads in transcripts
# v0.2		nov2016

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";


package htseqCount;

sub runSamtools(){
	my($samtoolsPath,$bamFile,$samFile,$executionCreatedTempDir)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	my $command="LC_ALL=C; export LC_ALL; export PATH=\$PATH:".$samtoolsPath."; ";
	#$command.="samtools view -h ".$bamFile." -o ".$samFile;
	
	# the awk command inserted in the line below is to avoid lines like the following:
	# HWI-D00689_0078:7:1106:8961:27289#29222_GCCAAT	272	*	156	1	29M	*	0	0	TTAAAATGAACCTGCCGGCTGATCGTTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:-5	XM:i:1	XO:i:0	XG:i:0	MD:Z:38C11	NM:i:1	XF:Z:1 ERCC-00096-chr5 156 29M50643494F21m TTAAAATGAACCTGCCGGCTGATCGTTTTTTTTAGGATATTGTGAGTAAT FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBB	NH:i:3	CC:Z:=	CP:i:156	XS:A:+	HI:i:0
	# these lines with '*' in the 3rd column cause htseqcount to break and stop execution (they are a very small number)
	
	# old way: removed because it does not sort reads in bam files by read name
	# This sorting is required by htseqcount for paired-end experiments, being mandatory that both reads of the pair
	# appear one after the other in the generated sam file
	#$command.="samtools view -h ".$bamFile." | awk '{if(\$3!=\"*\"){print \$0}}' > ".$samFile;
	
	
	# new way: sort reads by read name
	my @auxBAM=split('\/',$bamFile);
	my $auxBAMname=$auxBAM[@auxBAM-2]."_".$auxBAM[@auxBAM-1];
	my $auxSortedBam=$auxBAMname.".finallySorted";
	
	#This line below causes it to crash
	#$command.="samtools sort -no ".$bamFile." ".$executionCreatedTempDir."/temporarySORT_".$auxBAMname." | samtools view -F 4 - > ".$samFile;
	
	# Another way to avoid crashing
	$command.="samtools sort -no ".$bamFile." ".$executionCreatedTempDir."/temporarySORT_".$auxBAMname." > ".$executionCreatedTempDir."/temporarySORT_".$auxSortedBam."; ";
	$command.="samtools view ".$executionCreatedTempDir."/temporarySORT_".$auxSortedBam." | awk '{if(\$3!=\"*\"){print \$0}}' > ".$samFile."; ";
	
	$command.="rm -rf ".$executionCreatedTempDir."/temporarySORT_".$auxSortedBam;
	
	print "\n\t[executing] ".$command."\n";	
	system($command);	
	
}

sub runHtseqCount(){
#**********if problems running htseq-count, follow Simon Anders suggestion below:
#Install it locally to the user: 'python setup.py install --user'

	my($perl5lib,$extraPathsRequired,$mode,$minaqual,$featuretype,$idattr,$htseqcountPath,$htseqCountPythonpath,$GTF,
	$library,$samFile,$htseqCountOutputFile)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	#my $command="export PATH=\$PATH:".$htseqscount."; ";
	#$command.="htseq-count ";
	#$command.="python -m HTSeq.scripts.count ";

	my $command="LC_ALL=C; export LC_ALL; export PERL5LIB=".$perl5lib.":\$PERL5LIB; ";
        
        if($extraPathsRequired ne "NO_EXTRA_PATHS"){
		$extraPathsRequired=~ s/\'//g;
		
                my @exports=split(';',$extraPathsRequired);
                foreach my $export(@exports){
                        $command.="export ".$export."; ";
                }
        }
        
        $command.="export PYTHONPATH=\$PYTHONPATH:".$htseqCountPythonpath."; python ".$htseqcountPath."htseq-count ";
	
	#posibilities: unstranded, firststrand, secondstrand    
        if(lc($library) eq "unstranded"){
                $command.="--stranded no ";
        }elsif(lc($library) eq "firststrand"){
                $command.="--stranded reverse ";
	}elsif(lc($library) eq "secondstrand"){
                $command.="--stranded yes ";
	}

	$command.="--mode ".$mode." ";
	$command.="--minaqual ".$minaqual." ";	
	$command.="--type ".$featuretype." ";
	$command.="--idattr ".$idattr." ";
	$command.=$samFile." ";
	$command.=$GTF." ";
	$command.=" > ".$htseqCountOutputFile;
	
	print "\n\t[executing] ".$command."\n";	
	system($command);
	
	$command="rm -f ".$samFile;
	print "\n\t[deleting intermediate SAM file] ".$command."\n";	
	system($command);
}

1
