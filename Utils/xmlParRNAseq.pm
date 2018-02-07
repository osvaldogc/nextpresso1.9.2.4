#!/usr/bin/perl -w

# xmlParRNAseq.pm
# Author: Osvaldo Graña
# Description : RNA-seq analysis branch
# 
# v0.1		9feb2013
# v0.2		4ene2017 - all the code lines that are not needed (commented in previous versions) were deleted here


package xmlParRNAseq;

use strict;
use FindBin qw($Bin);
use XML::Simple;
use XML::Validator::Schema;
use lib "$Bin";


sub getDataFromConfigXMLdocument($){
	my ($elemento) = shift;		
	
	my $hashRef=printElement("Init",$elemento,"");
	
}

sub processXML($$){
	my ($XMLschema, $XMLdocument) = @_;
	
	if (!(-e $XMLschema)){
		print STDERR "\n ERROR processXML::$XMLschema does not exist";
		exit(-1);
	}
	if (!(-e $XMLdocument)){
		print STDERR "\n ERROR processXML::$XMLdocument does not exist";
		exit(-1);
	}

 
	# valida el documento XML con el schema
	my $validator = XML::Validator::Schema->new(file => $XMLschema);
	my $parser = XML::SAX::ParserFactory->parser(Handler => $validator);
	eval { $parser->parse_uri($XMLdocument) };
	die "File failed validation: $@" if $@;
	my $hashRef = XMLin($XMLdocument, forcearray=>1);
	
	return $hashRef;
}


1; # tells perl that the module was correctly loaded


