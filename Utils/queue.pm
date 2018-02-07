#!/usr/bin/perl
$!=1;
# FileName : 
# Author : Miriam Rubio
# Description :

package queue;

use strict;
use FindBin qw($Bin);
use lib "$Bin";


# Function: parseOut
# Description: This function extracts pids from qsub output sentences. 
# $filename: Path to input filename

sub parseOut($$) {
	my ($filename, $queueSystem) = @_;
	local(*FDes);	
	open(FDes,$filename)or die $!;
	my $line;
	my %pids=();
	while($line = <FDes>){
		if ($queueSystem eq 'PBS'){
			my @read = split(/\./,$line);
			#print "\n".$line."pid ".$read[0];
			$pids{$read[0]} = 1;#active
		}elsif ($queueSystem eq 'SGE'){
			my @read = split(" ",$line);
                	#print "\n".$line."pid ".$read[2];
               		$pids{$read[2]} = 1;#active
		}else{
			print "\n Error in queue system manager name", $queueSystem;
			exit();
		}

	}
	close(FDes);
	return \%pids;
}

# Function: waitPidsQUEUE
# Description: This function parses the pid file and returns the SGE wait parameter.
# $file_in: Path to input filename
sub waitPidsQUEUE($$){
	my ($file_in, $queueSystem) = @_;
	my $wait_pid = "";
	if ($queueSystem ne 'none'){
  	   my $pidsSystem = parseOut($file_in, $queueSystem);
  	   if (scalar(keys(%{$pidsSystem})) > 0){
  	  	   my $pids_list = join(",",keys(%{$pidsSystem}));
  	  	   if ($queueSystem eq 'SGE'){
  	  		   my $pids_list = join(",",keys(%{$pidsSystem}));
  	  		   $wait_pid ="-hold_jid $pids_list";
  	  	   }elsif ($queueSystem eq 'PBS'){
  	  		   my $pids_list = join(":",keys(%{$pidsSystem}));
  	  		   $wait_pid ="-W depend=afterany:$pids_list";
  	  	   }else{
  	  		   print "\n Error in queue system manager name", $queueSystem;
  	  		   exit();
  	  	   }
  	   }
  	   unlink($file_in);
  	   return ($wait_pid);
	}
}

sub createScript($$$$$$$){
	my ($queueSystem, $queueName, $queueSGEProject, $projectName, $out, $error, $commPerl) = @_;
	my $fileOut = "/tmp/script.".$ENV{USER}.".sh";
	local *FDout;
	open(FDout, ">".$fileOut) or die ("\n createScript :: Error creating $fileOut");
	my $command = "";
	if ($queueSystem eq 'PBS'){
	$command = qq(
#!/bin/sh

# Job Name
#PBS -N  $projectName

#PBS -q $queueName
##PBS -l ncpus=4
# my default stdout, stderr files
#PBS -o $out
#PBS -e $error
$commPerl
);
	}elsif ($queueSystem eq 'SGE'){
		if ($queueSGEProject ne 'none'){
	$command = qq(
#!/bin/csh

# Job name
#\$ -N $projectName 
#\$ -P $queueSGEProject

#Standard Error
#\$ -e $error

# Output
#\$ -o $out

#\$ -q $queueName

# Command
$commPerl
);
		}else{
$command = qq(
#!/bin/csh

# Job name
#\$ -N $projectName 

#Standard Error
#\$ -e $error

# Output
#\$ -o $out

#\$ -q $queueName

# Command
$commPerl
);		
		}
	}else{
		print "\n Error in queue system manager name", $queueSystem;
		exit();
	}
	print FDout $command;
	close(FDout);
	
	return $fileOut;
}

sub executeScript($$$$$$$$$$){
	my ($queueSystem, $queue, $queueSGEProject, $projName, $outFile, $errFile, $redirectOut,$commandPerl, $waitPid, $multiCFlag) = @_;
	my $command = "";
	my $scriptCommand = undef;
	if ($queueSystem  ne 'none'){
		$scriptCommand = createScript($queueSystem, $queue,  $queueSGEProject,$projName, $outFile, $errFile, $commandPerl);
		print"\n Perl command $commandPerl";
		$command = "qsub $waitPid $multiCFlag ".$scriptCommand.$redirectOut;
	}else{
		$command = $commandPerl;
	}
	print "Executed ".$command."\n\n";
	(system($command)) == 0 or die "[ERROR] nextpresso execution aborted, check log files please\n\n";#"\nfailed to execute: $!\n" ;
	if ($queueSystem  ne 'none'){
		unlink($scriptCommand);
	}
#	if ($? == -1) {
#		print "executeScript:: failed to execute: $!\n";
#	}
}

sub multicoreFlag($$){
	my ($queueSystem, $multicore) = @_;
	my $multiCFlag = "";
	if ($multicore != -1){
		if ($queueSystem eq 'SGE'){
			$multiCFlag = "-pe multicore $multicore";
		}elsif ($queueSystem eq 'PBS'){
			$multiCFlag = "-l select=1:ncpus=$multicore";
		}
	}
	return $multiCFlag;
}

1
