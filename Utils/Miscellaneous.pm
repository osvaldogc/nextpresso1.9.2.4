#!/usr/bin/perl -W

# Levels.pm 
# Author: Osvaldo Graña
# Description : contains different helpful functions

# v0.1		8mar2014

package Miscellaneous;
use strict;
use Date::Format;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function


sub getCurrentDateAndTime(){
	my ($second,$minute,$hour,$day, $mon, $year) = (localtime(time))[0..5];
	
	$mon+=1;
	$year+=1900;
	if($second<10){$second="0".$second}
	if($minute<10){$minute="0".$minute}
	if($mon<10){$mon="0".$mon}
	if($day<10){$day="0".$day}
	
	return "[".$day."-".$mon."-".$year.",".$hour.":".$minute.":".$second."] ";
}


1;
