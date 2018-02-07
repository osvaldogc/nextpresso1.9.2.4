#!/usr/bin/perl
# FileName :
# Author : Miriam Rubio
# Description:

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";
use File::Path qw(make_path);

package filesFunc;

sub getFiles($);
sub is_casava1_8($$);
sub makeDir($);

sub getFiles($){
	my ($directory) = @_;
	(-d $directory) or die "\n getFiles:: Directory $directory does not exist";
	my @entries = <$directory/*>;
	my @fileSet = ();
	foreach my $entry (@entries) {
		if ( -f $entry){
			push (@fileSet,$entry);
			#print $entry . "\n";
		}elsif ( -d $entry){
			getFiles($entry);
		}
	} 
 	return \@entries;
 }
 
# Function: is_casava1_8
# Description:This function returns 1 if $file_in is in casava format
# $filein: Input file name
sub is_casava1_8($$){
	my ($filein, $checkCasava) =@_;
	if ($checkCasava eq '0'){
		#Quality in ascii33, headers without specific format.
		return 2;
	}
	local *FDin;
	my $line = undef;
	open (FDin , "<".$filein) or die "\n is_casava1_8::Can't open file $filein. Reason: $!";
	$line=<FDin>;
	close(FDin);
	if ($line !~ m{^@.+:.+:.+:.+:.+:.+:.+:[YN]:*}){
	 print"\n INFO is_casava1_8:: Fastq file $filein is not in 1.8 version";
	 return 0;
	}
	print"\n INFO is_casava1_8:: Fastq file $filein in 1.8 version";
	return 1;
}



# Function: makeDir
# Auxiliary function.It creates a $dir directory if it doesn't exist.
# $dir: Directory name.
sub makeDir($){
	my ($dir) = @_;
	if(-d $dir){
		print"\n Directory $dir exists";
	}else{
		File::Path::make_path($dir);
	}
}
 
1
