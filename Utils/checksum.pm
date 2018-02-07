#!/usr/bin/perl

# Author : Osvaldo Grana
# Description: perform sample validation doing checksum
# v0.1		may2016

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";

package checksum;

sub performChecksum($$){
	my($inputFile,$outputFile)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	#if first checks if 'md5sum' is installed
	my $command="md5sum --version";
	system($command);
	if($? == 0){
		$command="md5sum ".$inputFile." >> ".$outputFile;	
		print "\n\t[executing] ".$command."\n";	
		system($command);
	}else{
		print "\n[ERROR checksum.pl]:\n";
		print "\tmd5sum is not installed, cannot perform check\n\n\n";			
	}
}

sub validate($$){
	my($fileWithChecksumCodesToValidate,$samples)=@_;
	
	my @files=split(',',$samples);
	
	#just in case that the different samples were stored in different places,
	#(it would have to perform checksum validation in each one of these directories)
	my $lastPath="";;
	foreach my $file(@files){
		chomp($file);
		
		my @tokens=split(/\//,$file);
		
		#retrieves the file path of each sample
		my $filePath="";
		for(my $i=0;$i<@tokens-1;$i++){
			$filePath.="/".$tokens[$i];
			$filePath=~ s/\/\//\//g;
		}
	
		if($lastPath!~ /$filePath/){
			print "\n\t[checksum:validate()] setting path to verify samples checksum: ".$filePath;
			my $command="cd ".$filePath."; ";
			$command.="md5sum -c --warn --strict ".$fileWithChecksumCodesToValidate;	
			print "\n\t[executing] ".$command."\n";	
			system($command);
			
			#if($? != 0){
#				print STDERR "\n[ERROR checksum.pl]:\n";
#				print STDERR $?."\n";
#				print STDERR "[Execution stopped]\n\n";
#				exit(-1);					
#			}
			
			$lastPath.="=".$filePath."=";
		}
		
	}

}

1
