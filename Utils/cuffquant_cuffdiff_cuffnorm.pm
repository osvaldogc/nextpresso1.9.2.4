#!/usr/bin/perl

# Author : Osvaldo Grana
# Description: runs cuffquant, cuffdiff, cuffnorm and prepares preranked 'rnk' files for GSEA
# v0.3		oct2016
# v0.4		dic2016 - adds cumulative variance for the samples PCAs
#			- adds colors to sample names in PCA plots
# v0.5		oct2017 - Bug solved in PCA construction
# v1.9.2	ene2018 - added scale=TRUE in prcomp function used to build PCA
#			  new XLSX files are now created without background genes and flat pattern genes
#			  added new output file with PCA component values

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";


package cuffquant_cuffdiff_cuffnorm;

sub runCuffquant(){
	my($extraPathsRequired,$cufflinksPath,$samtoolsPath,$alignmentsDir,$referenceSequence,$samples,$GTF,$libraryType,$nThreads,
		$cuffquantNThreads,$cuffquantSeed,$cuffquantFragBiasCorrect,$cuffquantMultiReadCorrect,$cuffquantOutDir,$cuffquantMaxBundleFrags,
		$cuffquant_noEffectiveLengthCorrection,$cuffquant_noLengthCorrection,$inputFile)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	my $command="export PATH=\$PATH:".$samtoolsPath.":".$cufflinksPath."; ";
	
	if($extraPathsRequired ne "NO_EXTRA_PATHS"){
		$extraPathsRequired=~ s/\'//g;
		
		my @exports=split(';',$extraPathsRequired);
		foreach my $export(@exports){
			$command.="export ".$export."; ";
		}
	}
	
	$command.="cuffquant -p ".$cuffquantNThreads." ";
	
	#posibilities: unstranded, firststrand, secondstrand    
        if(lc($libraryType) eq "unstranded"){
                $command.="--library-type fr-unstranded ";
        }elsif(lc($libraryType) eq "firststrand"){
                $command.="--library-type fr-firststrand ";
        }elsif(lc($libraryType) eq "secondstrand"){
                $command.="--library-type fr-secondstrand ";
        }

	if($cuffquantFragBiasCorrect eq "true"){
		$command.="--frag-bias-correct ".$referenceSequence." ";
	}

	if($cuffquantMultiReadCorrect eq "true"){
		$command.="--multi-read-correct ";
	}
	
	if($cuffquant_noEffectiveLengthCorrection eq "true"){
		$command.="--no-effective-length-correction ";
	}
	
	if($cuffquant_noLengthCorrection eq "true"){
		$command.="--no-length-correction ";
	}	
	
	$command.="--max-bundle-frags ".$cuffquantMaxBundleFrags." ";
	$command.="--seed ".$cuffquantSeed." ";
	$command.="-o ".$cuffquantOutDir.$inputFile." ";
	
	$command.=$GTF." ";
	$command.=$alignmentsDir.$inputFile."/accepted_hits.bam";
	
	print "\n\t[executing] ".$command."\n";	
	system($command);
	
}

sub runCuffdiff(){
	my($extraPathsRequired,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$cuffdiffOutDir,$GTF,$cuffdiffLibraryType,$cuffdiffNThreads,
		$cuffdiffFragBiasCorrect,$cuffdiffMultiReadCorrect,$cuffdiffMinAlignmentCount,$cuffdiffSeed,
		$cuffdiffInputFiles,$referenceSequence,$cuffdiffFDR,$cuffdiffLibraryNormalizationMethod,$cuffdiffComparisonLabels,$cuffdiffMaxBundleFrags,
		$cuffdiff_noEffectiveLengthCorrection,$cuffdiff_noLengthCorrection,$cuffdiff_dispersionMethod)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	my $command="export PATH=\$PATH:".$samtoolsPath.":".$cufflinksPath."; ";
	
	if($extraPathsRequired ne "NO_EXTRA_PATHS"){
		$extraPathsRequired=~ s/\'//g;
		
		my @exports=split(';',$extraPathsRequired);
		foreach my $export(@exports){
			$command.="export ".$export."; ";
		}
	}
	
	$command.="cuffdiff -p ".$cuffdiffNThreads." ";
	$command.="-c ".$cuffdiffMinAlignmentCount." ";
	
	#posibilities: unstranded, firststrand, secondstrand    
        if(lc($cuffdiffLibraryType) eq "unstranded"){
                $command.="--library-type fr-unstranded ";
        }elsif(lc($cuffdiffLibraryType) eq "firststrand"){
                $command.="--library-type fr-firststrand ";
        }elsif(lc($cuffdiffLibraryType) eq "secondstrand"){
                $command.="--library-type fr-secondstrand ";
        }

	if($cuffdiffFragBiasCorrect eq "true"){
		$command.="--frag-bias-correct ".$referenceSequence." ";
	}

	if($cuffdiffMultiReadCorrect eq "true"){
		$command.="--multi-read-correct ";
	}

	if($cuffdiff_noEffectiveLengthCorrection eq "true"){
		$command.="--no-effective-length-correction ";
	}
	
	if($cuffdiff_noLengthCorrection eq "true"){
		$command.="--no-length-correction ";
	}	
	
	$command.="--dispersion-method ".$cuffdiff_dispersionMethod." ";
	$command.="--max-bundle-frags ".$cuffdiffMaxBundleFrags." ";
	$command.="--seed ".$cuffdiffSeed." ";
	$command.="--FDR ".$cuffdiffFDR." ";
	$command.="--library-norm-method ".$cuffdiffLibraryNormalizationMethod." ";
	$command.="-o ".$cuffdiffOutDir." "; # includes the comparison name (subdir)
	$command.="--labels ".$cuffdiffComparisonLabels." "; 
	$command.=$GTF." ";
	$command.=$cuffdiffInputFiles;
	
	print "\n\t[executing] ".$command."\n";	
	system($command);
	
}

sub runCuffnorm(){
	my($extraPathsRequired,$cuffnormNThreads,$cufflinksPath,$samtoolsPath,$cuffquantOutDir,$GTF,$cuffnormOutputFormat,
		$cuffnormLibraryNormalizationMethod,$cuffnormSeed,$cuffnormNormalization,
		$cuffnormLibraryType,$cuffnormOutDir,$cuffnormLabels,$cuffnormInputFiles)=@_;

	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	my $command="export PATH=\$PATH:".$samtoolsPath.":".$cufflinksPath."; ";
	
	if($extraPathsRequired ne "NO_EXTRA_PATHS"){
		$extraPathsRequired=~ s/\'//g;
		
		my @exports=split(';',$extraPathsRequired);
		foreach my $export(@exports){
			$command.="export ".$export."; ";
		}
	}
	
	$command.="cuffnorm -p ".$cuffnormNThreads." ";
	
	#posibilities: unstranded, firststrand, secondstrand    
        if(lc($cuffnormLibraryType) eq "unstranded"){
                $command.="--library-type fr-unstranded ";
        }elsif(lc($cuffnormLibraryType) eq "firststrand"){
                $command.="--library-type fr-firststrand ";
        }elsif(lc($cuffnormLibraryType) eq "secondstrand"){
                $command.="--library-type fr-secondstrand ";
        }

	$command.="--seed ".$cuffnormSeed." ";
	$command.="--library-norm-method ".$cuffnormLibraryNormalizationMethod." ";
	if($cuffnormNormalization eq "compatibleHits"){
		$command.="--compatible-hits-norm ";
	}else{
		$command.="--total-hits-norm ";	
	}
	
	$command.="-o ".$cuffnormOutDir." "; # includes the comparison name (subdir)
	$command.="--labels ".$cuffnormLabels." "; 
	$command.=$GTF." ";
	$command.=$cuffnormInputFiles;
	
	print "\n\t[executing] ".$command."\n";	
	system($command);

}

sub createExcel(){
	my($cuffdiffOutDir,$FDR,$cuffdiffComparisonLabels,$FPKMthreshold,$executionCreatedTempDir)=@_;
	
	use Excel::Writer::XLSX;	
	

	########## GENE LEVEL ##########	

	my($file,$outFile,$outFilteredFile);
	$file=$cuffdiffOutDir."gene_exp.diff";
	
	if($cuffdiffOutDir!~ /cuffdiff_backgroundFiltered_AND_flatPatternFiltered/){		
		$outFile=$cuffdiffOutDir;
		$outFilteredFile=$cuffdiffOutDir;
		$outFile=~ s/\/$//;
		$outFile.=".gene_exp.xlsx";
		$outFilteredFile=~ s/\/$//;
		
		#this file is finally removed in this new version
		$outFilteredFile.=".gene_exp_FPKMthreshold_FILTERED.xlsx";
	}else{
		$outFile=$cuffdiffOutDir;
		$outFilteredFile=$cuffdiffOutDir;
		$outFile=~ s/\/$//;
		$outFile.=".gene_exp_backgroundANDflatpatterns_FILTERED.diff.xlsx";
		$outFilteredFile=~ s/\/$//;
		
		#this file is finally removed in this new version
		$outFilteredFile.=".gene_exp_backgroundANDflatpatternsFiltered.diff.FPKMthreshold_FILTERED.xlsx";	
	}
	
	#the input files are sorted by q-value before doing the work,
	#so the output files are already sorted by q-value ascending
	my $fileSORTED=$file.".SORTED";
	system("head -n 1 ".$file." > ".$fileSORTED);
	#LC_ALL=C -> forces the local configuration to use '.' to represent decimal numbers
	system("LC_ALL=C; export LC_ALL; grep -v 'value_1'  ".$file." | sort -k13g >> ".$fileSORTED);
	
	#finds out what is the case/condition name (to write the legend)
	my $case=(split(',',$cuffdiffComparisonLabels))[1];	
	
	#opens gene_exp file
	if(open(IN,$fileSORTED)){

		#creates a new Excel workbook
    		my $workbook=Excel::Writer::XLSX->new($outFile);
		my $workbook2= Excel::Writer::XLSX->new($outFilteredFile);
		
		$workbook->set_tempdir($executionCreatedTempDir);
		$workbook2->set_tempdir($executionCreatedTempDir);
		
    		#adds a worksheet
    		my $worksheet = $workbook->add_worksheet();
		my $worksheet2 = $workbook2->add_worksheet();
		
    		#add and define a format
    		my $format;
		my $format2;
		my $row = 0;
		
		#legend Up/Down
		$format= $workbook->add_format();
		$format->set_color( 'red' );
		$format->set_bold(0);  # Turn bold on
		$format->set_bg_color( 'yellow' );
		$worksheet->merge_range($row,0,$row,2,"Upregulated in ".$case, $format );
		$format2= $workbook2->add_format();
		$format2->set_color( 'red' );
		$format2->set_bold(0);  # Turn bold on
		$format2->set_bg_color( 'yellow' );
		$worksheet2->merge_range($row,0,$row,2,"Upregulated in ".$case, $format2 );
		$row++;
		
		$format= $workbook->add_format();
		$format->set_color( 'green' );
		$format->set_bold(0);  # Turn bold on
		$format->set_bg_color( 'yellow' );
		$worksheet->merge_range($row,0,$row,2,"Downregulated in ".$case, $format );
		$format2= $workbook2->add_format();
		$format2->set_color( 'green' );
		$format2->set_bold(0);  # Turn bold on
		$format2->set_bg_color( 'yellow' );
		$worksheet2->merge_range($row,0,$row,2,"Downregulated in ".$case, $format2 );
		$row++;		
		
		$format= $workbook->add_format();
		$format->set_color( 'black' );
		$format->set_bold(1);  # Turn bold on
		$format->set_bg_color(-1);
		$worksheet->merge_range($row,0,$row,4,"FDR=".$FDR, $format );
		$format2= $workbook2->add_format();
		$format2->set_color( 'black' );
		$format2->set_bold(1);  # Turn bold on
		$format2->set_bg_color(-1);
		$worksheet2->merge_range($row,0,$row,4,"FDR=".$FDR, $format2 );		
		$row++;
		$row++;
		$row++;
		
		my $row2=$row;# to be used by format2
		
		my $memoryFirstRow=0;
		my $memoryLastRow=0;
		my $memoryFirstCol=0;
		my $memoryLastCol=0;	
			
		foreach my $line(<IN>){
			chomp($line);
			my @tokens=split('\t+',$line);

			#header line
			if($line=~ /sample_/){
				$format= $workbook->add_format();
				$format->set_color( 'black' );
				$format->set_bold(1);  # Turn bold on
				
				$format2= $workbook2->add_format();
				$format2->set_color( 'black' );
				$format2->set_bold(1);  # Turn bold on					
				
				$tokens[7]="FPKM_1";
				$tokens[8]="FPKM_2";
				
				$memoryFirstRow=$row;
				$memoryLastRow=$row;
				$memoryFirstCol=0;
				$memoryLastCol=@tokens-1;
				
			}else{#data line

				$format= $workbook->add_format();
				$format->set_bold(0);  # Turn bold on
				$format2= $workbook2->add_format();
				$format2->set_bold(0);  # Turn bold on

				#if log2(fold-change)>=0, then red font
				if($tokens[9]!~ /^-/){
					$format->set_color( 'red' );
					$format2->set_color( 'red' );

				}else{#green font
					$format->set_color( 'green' );
					$format2->set_color( 'green' );
				}

				#yellow background if significant difference
				if($line=~ /yes$/){
					$format->set_bg_color( 'yellow' );
					$format2->set_bg_color( 'yellow' );
				}		
				
			}

			my $done=0;
			
			#prints line with data			
			for(my $col=0;$col<@tokens;$col++){
				$worksheet->write( $row, $col,$tokens[$col], $format );
				
				#FILTERED
				if(($tokens[6] eq "OK") && ($tokens[7]>$FPKMthreshold || $tokens[8]>$FPKMthreshold)){
					$worksheet2->write( $row2, $col,$tokens[$col], $format2 );
					$done=1;					
				}elsif($line=~ /sample_1/){
					$worksheet2->write( $row2, $col,$tokens[$col], $format2 );
					$done=1;
				}
			}

			$row++;
			
			if($done){
				$row2++;
			}	
		}
		close(IN);
				
		#adds autofilter		
		$worksheet->autofilter($memoryFirstRow,$memoryFirstCol,$memoryLastRow,$memoryLastCol);
		$worksheet2->autofilter($memoryFirstRow,$memoryFirstCol,$memoryLastRow,$memoryLastCol);
		
		#set columns width
		$worksheet->set_column('A:B' , 10 );
		$worksheet->set_column('C:C' , 15 );
		$worksheet->set_column('D:D' , 7);
		$worksheet->set_column('E:F' ,12);	
		$worksheet->set_column('G:G' , 8 );	
		$worksheet->set_column('H:I' , 10 );
		$worksheet->set_column('J:J' , 19);
		$worksheet->set_column('K:M' , 10 );
		$worksheet->set_column('N:N' , 13);

		$worksheet2->set_column('A:B' , 10 );
		$worksheet2->set_column('C:C' , 15 );
		$worksheet2->set_column('D:D' , 7);
		$worksheet2->set_column('E:F' ,12);	
		$worksheet2->set_column('G:G' , 8 );	
		$worksheet2->set_column('H:I' , 10 );
		$worksheet2->set_column('J:J' , 19);
		$worksheet2->set_column('K:M' , 10 );
		$worksheet2->set_column('N:N' , 13);

			
		$workbook->close();
		$workbook2->close();
		
		#this file is not required after removing background genes and flat patterns genes
		#it is then removed
		system("rm -f ".$outFilteredFile);		
		
	}else{
		print STDERR "\n[ERROR createExcel]: $file file does not exist\n\n";
		exit(-1);
	}
	 



	########## ISOFORM LEVEL ##########
	$file=$cuffdiffOutDir."isoform_exp.diff";
	
	if($cuffdiffOutDir!~ /cuffdiff_backgroundFiltered_AND_flatPatternFiltered/){		
		$outFile=$cuffdiffOutDir;
		$outFilteredFile=$cuffdiffOutDir;
		$outFile=~ s/\/$//;
		$outFile.=".isoform_exp.xlsx";
		$outFilteredFile=~ s/\/$//;
		
		#this file is finally removed in this new version
		$outFilteredFile.=".isoform_exp_FPKMthreshold_FILTERED.xlsx";
	}else{
		$outFile=$cuffdiffOutDir;
		$outFilteredFile=$cuffdiffOutDir;
		$outFile=~ s/\/$//;
		$outFile.=".isoform_exp_backgroundANDflatpatterns_FILTERED.diff.xlsx";
		$outFilteredFile=~ s/\/$//;
		
		#this file is finally removed in this new version
		$outFilteredFile.=".isoform_exp_backgroundANDflatpatternsFiltered.diff.FPKMthreshold_FILTERED.xlsx";	
	}

	#the input files are sorted by q-value before doing the work,
	#so the output files are already sorted by q-value ascending
	$fileSORTED=$file.".SORTED";
	system("head -n 1 ".$file." > ".$fileSORTED);
	#LC_ALL=C -> forces the local configuration to use '.' to represent decimal numbers
	system("LC_ALL=C; export LC_ALL; grep -v 'value_1'  ".$file." | sort -k13g >> ".$fileSORTED);
		
	#opens gene_exp file
	if(open(IN,$fileSORTED)){

		#creates a new Excel workbook
    		my $workbook=Excel::Writer::XLSX->new($outFile);
		my $workbook2= Excel::Writer::XLSX->new($outFilteredFile);
    		#adds a worksheet
    		my $worksheet = $workbook->add_worksheet();
		my $worksheet2 = $workbook2->add_worksheet();
		
    		#add and define a format
    		my $format;
		my $format2;
		my $row = 0;
		
		#legend Up/Down
		$format= $workbook->add_format();
		$format->set_color( 'red' );
		$format->set_bold(0);  # Turn bold on
		$format->set_bg_color( 'yellow' );
		$worksheet->merge_range($row,0,$row,4,"Upregulated in ".$case, $format );
		$format2= $workbook2->add_format();
		$format2->set_color( 'red' );
		$format2->set_bold(0);  # Turn bold on
		$format2->set_bg_color( 'yellow' );
		$worksheet2->merge_range($row,0,$row,4,"Upregulated in ".$case, $format2 );		
		$row++;
		
		$format= $workbook->add_format();
		$format->set_color( 'green' );
		$format->set_bold(0);  # Turn bold on
		$format->set_bg_color( 'yellow' );
		$worksheet->merge_range($row,0,$row,4,"Downregulated in ".$case, $format );
		$format2= $workbook2->add_format();
		$format2->set_color( 'green' );
		$format2->set_bold(0);  # Turn bold on
		$format2->set_bg_color( 'yellow' );
		$worksheet2->merge_range($row,0,$row,4,"Downregulated in ".$case, $format2 );					
		$row++;		
		
		$format= $workbook->add_format();
		$format->set_color( 'black' );
		$format->set_bold(1);  # Turn bold on
		$format->set_bg_color(-1);
		$worksheet->merge_range($row,0,$row,4,"FDR=".$FDR, $format );
		$format2= $workbook2->add_format();
		$format2->set_color( 'black' );
		$format2->set_bold(1);  # Turn bold on
		$format2->set_bg_color(-1);
		$worksheet2->merge_range($row,0,$row,4,"FDR=".$FDR, $format2 );		
		$row++;
		$row++;
		$row++;
		
		my $row2=$row; # to be used by format2

		my $memoryFirstRow=0;
		my $memoryLastRow=0;
		my $memoryFirstCol=0;
		my $memoryLastCol=0;	
				
		foreach my $line(<IN>){
			chomp($line);
			my @tokens=split('\t+',$line);

			#header line
			if($line=~ /sample_/){
				$format= $workbook->add_format();
				$format->set_color( 'black' );
				$format->set_bold(1);  # Turn bold on
				
				$format2= $workbook2->add_format();
				$format2->set_color( 'black' );
				$format2->set_bold(1);  # Turn bold on					
				
				$tokens[7]="FPKM_1";
				$tokens[8]="FPKM_2";

				$memoryFirstRow=$row;
				$memoryLastRow=$row;
				$memoryFirstCol=0;
				$memoryLastCol=@tokens-1;				
				
			}else{#data line

				$format= $workbook->add_format();
				$format->set_bold(0);  # Turn bold on
				$format2= $workbook2->add_format();
				$format2->set_bold(0);  # Turn bold on

				#if log2(fold-change)>=0, then red font
				if($tokens[9]!~ /^-/){
					$format->set_color( 'red' );
					$format2->set_color( 'red' );

				}else{#green font
					$format->set_color( 'green' );
					$format2->set_color( 'green' );
				}

				#yellow background if significant difference
				if($line=~ /yes$/){
					$format->set_bg_color( 'yellow' );
					$format2->set_bg_color( 'yellow' );
				}		
				
			}

			my $done=0;
			
			#prints line with data			
			for(my $col=0;$col<@tokens;$col++){
				$worksheet->write( $row, $col,$tokens[$col], $format );
				
				#FILTERED
				if(($tokens[6] eq "OK") && ($tokens[7]>$FPKMthreshold || $tokens[8]>$FPKMthreshold)){
					$worksheet2->write( $row2, $col,$tokens[$col], $format2 );
					$done=1;					
				}elsif($line=~ /sample_1/){
					$worksheet2->write( $row2, $col,$tokens[$col], $format2 );
					$done=1;
				}
			}

			$row++;
			
			if($done){
				$row2++;
			}	
		}
		close(IN);
		
		#adds autofilter		
		$worksheet->autofilter($memoryFirstRow,$memoryFirstCol,$memoryLastRow,$memoryLastCol);
		$worksheet2->autofilter($memoryFirstRow,$memoryFirstCol,$memoryLastRow,$memoryLastCol);
		
		#set columns width
		$worksheet->set_column('A:B' , 10 );
		$worksheet->set_column('C:C' , 15 );
		$worksheet->set_column('D:D' , 7);
		$worksheet->set_column('E:F' ,12);	
		$worksheet->set_column('G:G' , 8 );	
		$worksheet->set_column('H:I' , 10 );
		$worksheet->set_column('J:J' , 19);
		$worksheet->set_column('K:M' , 10 );
		$worksheet->set_column('N:N' , 13);

		$worksheet2->set_column('A:B' , 10 );
		$worksheet2->set_column('C:C' , 15 );
		$worksheet2->set_column('D:D' , 7);
		$worksheet2->set_column('E:F' ,12);	
		$worksheet2->set_column('G:G' , 8 );	
		$worksheet2->set_column('H:I' , 10 );
		$worksheet2->set_column('J:J' , 19);
		$worksheet2->set_column('K:M' , 10 );
		$worksheet2->set_column('N:N' , 13);
		
		$workbook->close();
		$workbook2->close();
		
		#this file is not required after removing background genes and flat patterns genes
		#it is then removed
		system("rm -f ".$outFilteredFile);		
		
	}else{
		print STDERR "\n[ERROR createExcel]: $file file does not exist\n\n";
		exit(-1);
	}	 
}

sub createGSEArnkFile(){
# CONDITIONS:
#
# RNK file is now created from log2 fold-change (instead of statistic)
# 1. When creating the ranked file (rnk) for GSEA, all FPKM values are transformed as follows FPKM=original_FPKM+1. This transformation affects both FPKM1 and FPKM2.
# It tries to avoid cases like the following:
# test_id	gene_id	gene	locus	sample_1	sample_2	status	FPKM_1	FPKM_2	log2(fold_change)	test_stat	p_value	q_value	significant
# ARHGAP19-SLIT1	ARHGAP19-SLIT1	ARHGAP19-SLIT1	chr10:98757794-99052430	Control	Treatment	NOTEST	1.10925E-42	5.20103E-136	-310.032	0	1	1	no
# where two real very low FPKMs give raise to an unreal and extremely negative and log2(Fold-change), that in GSEA is placed on the edge of one of the phenotypes, representing and artifact.
# No other filtering is considered (like the one used in version 1.5 for cases with FPKM=0 in both conditions).
#
# 2. Transcript annotations where the names starts with MIR, SNORD, SNORA, and SCARNA are going to be filtered out, as they correspond to transcripts
# that are not properly detected by RNAseq, and, on the other hand, not properly quantified by cufflinks (lower than 1KB)

	my($cuffdiffOutDir)=@_;


	########## GENE LEVEL ##########	

	my $inFile=$cuffdiffOutDir."gene_exp.diff";
	my $outFile=$cuffdiffOutDir;
	$outFile=~ s/\/$//;
	$outFile.=".rnk";
	my $outFileaux=$outFile.".aux";
	
	#my $inFileRESORTED=$inFile.".RESORTED";
	
	#it first does some kind of 'random sorting', to alter inital lines sorted by gene name
	#my $command1="LC_ALL=C; export LC_ALL; sort --random-sort ".$inFile." -o ".$inFileRESORTED;	
	#print "\n\t[executing] ".$command1."\n";
	#system($command1);
	
	#opens gene_exp file RESORTED
	#if(open(IN,$inFileRESORTED)){
	if(open(IN,$inFile)){		
		my $output="";

		#runs through the expression file
		foreach my $line(<IN>){
			chomp($line);
			
			#if not in header line
			if($line!~ /test_stat/){			
				my @tokens=split('\t',$line);
				my $geneID=$tokens[2];
		
				#filters out MIR, SNORD, SNORA, and SCARNA transcripts
				if(uc($geneID)!~ /^MIR/ && uc($geneID)!~ /^SNORD/ && uc($geneID)!~ /^SNORA/ && uc($geneID)!~ /^SCARNA/){						

					my $FPKM1=$tokens[7]+1;
					my $FPKM2=$tokens[8]+1;
					
					my $log2_FoldChange=(log($FPKM2/$FPKM1))/(log(2));
					
					$output.=uc($geneID)."\t".$log2_FoldChange."\n";
					
				}
			}
		}
		close(IN);
		
		#saves the results in the auxiliar file
		open(OUT,">",$outFileaux);
		print OUT $output;
		close(OUT);
		
		#sorts the auxiliar file, saves the sorted file in the output file and deletes the auxiliar file
		#LC_ALL=C -> forces the local configuration to use '.' to represent decimal numbers
		#--stable avoids default (unrequested) sorting by gene name when log2foldchange=0
		system("LC_ALL=C; export LC_ALL; sort -k2g --stable ".$outFileaux." -o ".$outFile."; rm -f ".$outFileaux);
	}
	else{
		print STDERR "\n[ERROR createGSEArnkFiles_from_Cuffdiff]: ".$inFile." file does not exist\n\n";
		exit(-1);		
	}
	
	#deletes the auxiliar file
	#$command1="rm -f ".$inFileRESORTED;
	#print "\n\t[executing] ".$command1."\n";
	#system($command1);
}

sub calculateCorrelationsAndPCA_GeneLevel(){
	my($extraPathsRequired,$cuffnormOutDir,$cuffnormInputFiles,$cuffnormColors,$executionCreatedTempDir)=@_;

	#checks out the names of the comparisons that there are in the cuffnorm directory
	open(IN,"ls -d ".$cuffnormOutDir." | ");
	my @files=<IN>;
	close(IN);
	
	chomp(@files);	
	#performs everything for each comparison that there is in the cuffnorm directory
	foreach my $file(@files){
		chomp($file);
		
		my $prefix=$file;
		$prefix=~ s/\/$//;
		
		#as 'ls' was executed with -d above, the answer already includes the path to the file and the '/' at the end.
		#for example:
		#/mnt/bigbrother/ograna/Analysis/cuffnorm/Meki_vs_DMSO/
		#so there is no need to add it below when defining $path=$file

		my %names=();		

		#OPENS THE TABLE with the sample names
		my $numberOfSamples=0; #tracks the number of samples that there are to create the gct file
		open(IN,$file."samples.table");	
		foreach my $line(<IN>){
			chomp($line);
		
			if($line!~ /sample_id/){
				$numberOfSamples++;
				my $cuffnormName=(split('\t',$line))[0];
				my $originalName=(split('\t',$line))[1];
			
				$originalName=(split('cuffquant\/',$originalName))[1];
				$originalName=(split('\/abundances\.cxb',$originalName))[0];
			
				#recover the names
				$names{$cuffnormName}=$originalName;					
			}#if($line!~ /sample_id/)
		}#foreach my $line(<IN>)			
		close(IN);	
					
               #CREATES NEW FILES with original (appropriated) names
               open(HEADER,"head -n 1 ".$file."genes.count_table |");
               my $header=<HEADER>;
               close(HEADER);
               chomp($header);
               my @tokens=split('\t',$header);
               	      
	       #prepares new header
               my $newHeader="Name";
	       my $newHeaderForGCTfile=$newHeader."\tDescription";
               for(my $p=1; $p<@tokens;$p++){
	       		my $head=$tokens[$p];
			chomp($head);
			
			$newHeader.="\t".$names{$head};
			$newHeaderForGCTfile.="\t".$names{$head};
		}
		$newHeader.="\n";
 		$newHeaderForGCTfile.="\n";
 		
#		my $newGeneCountTable=$prefix.".genes.count_table.xls";
#		open(OUT,">",$newGeneCountTable);
#		print OUT $newHeader;
#		close(OUT);
#		
#		print "\n\t[creating] ".$newGeneCountTable."\n";
#		system("grep -v 'tracking_id' ".$file."genes.count_table >> ".$newGeneCountTable);
		
		#counts the number of genes that there are to create the gct file
		open(IN,"grep -v 'tracking_id' ".$file."genes.fpkm_table | wc -l |");
		my $nlines=<IN>;
		close(IN);
		chomp($nlines);
		my $nGenes=(split(' ',$nlines))[0];
				
		$newHeaderForGCTfile="#1.2\n".$nGenes."\t".$numberOfSamples."\n".$newHeaderForGCTfile;
		my $newFPKMTable=$prefix.".genes.fpkm_table.xls";
		open(OUT,">",$newFPKMTable);
		print OUT $newHeader;
		close(OUT);
		
		print "\n\t[creating] ".$newFPKMTable."\n";
		system("grep -v 'tracking_id' ".$file."genes.fpkm_table >> ".$newFPKMTable);
		
		#creates a GCT FPKM file for gene-e in gct file format
		my $FPKMgctFile=$newFPKMTable;
		$FPKMgctFile=~ s/xls/gct/;
		
		open(IN,$newFPKMTable);
		open(OUT,">",$FPKMgctFile);
		print OUT $newHeaderForGCTfile;
		#for each line of the FPKM file
		foreach my $line(<IN>){
			chomp($line);
			
			if($line!~ /^Name/){
				@tokens=split('\t',$line);			
			
				for(my $k=0;$k<@tokens;$k++){
					my $token=$tokens[$k];
					chomp($token);

					if($k<@tokens-1){
						print OUT $token."\t";
						
						#just after first column adds no description
						if($k==0){
							print OUT "No description added\t";
						}
					}else{
						print OUT $token."\n";
					}

				}#for(my $k=0;$k<@tokens;$k++)
			}#if($line!~ /^Name/)		
		}
		close(IN);
		close(OUT);
		

		# CORRELATION TABLE AND PCA WITH R

		#first creates R script
		my $Rscriptfile=$prefix.".correlationAndPCA.R";
		my $corrFile=$prefix.".pearsonCorrelationsAmongSamples.xls";
		my $variance=$prefix.".proportionOfVariance.pdf";
		my $PCAfile=$prefix.".samplesPCA.pdf";
		my $cumVarFile=$prefix.".cumulativeVariance.xls";
		my $PCAcompFile=$prefix.".PCAcomponent_values.xls";
		
		open(RSCRIPT,">",$Rscriptfile);	
		print RSCRIPT "FPKMs=read.table(\"".$newFPKMTable."\",header=T)\n";
		print RSCRIPT "dim(FPKMs)\n";
		print RSCRIPT "FPKMmatrix <- as.matrix(FPKMs[,2:".($numberOfSamples+1)."])\n";	
		#pearsonCorrelation
		print RSCRIPT "FPKMcorr <- cor(FPKMmatrix,method=\"pearson\")\n";
		print RSCRIPT "write.table(FPKMcorr,file=\"".$corrFile."\",sep=\"\\t\",col.names=NA)\n";
		#PCA
		#print RSCRIPT "FPKMpc <- prcomp(t(FPKMmatrix),scale=TRUE)\n";
		print RSCRIPT "possibleError<-tryCatch(FPKMpc <- prcomp(t(FPKMmatrix),scale=TRUE),error=function(e) e)\n";
		print RSCRIPT "if(inherits(possibleError, \"error\")){FPKMpc <- prcomp(t(FPKMmatrix),scale=FALSE)}\n";
		#proportion of variance
		print RSCRIPT "pdf(\"".$variance."\")\n";
		print RSCRIPT "screeplot(FPKMpc)\n";
		print RSCRIPT "dev.off()\n";
		#Cumulative variance
		print RSCRIPT "eig <- (FPKMpc\$sdev)^2\n";
		print RSCRIPT "variance <- eig*100/sum(eig)\n";
		print RSCRIPT "cumvar <- cumsum(variance)\n";
		print RSCRIPT "write.table(cumvar,file=\"".$cumVarFile."\",sep=\"\\t\",col.names=NA)\n";		
		#PCA graph
		print RSCRIPT "pdf(\"".$PCAfile."\")\n";
		print RSCRIPT "PCAcomponent <- FPKMpc\$x\n";		

		#defines colors
		#my @colors=("lightblue","lightgreen","lightpink","orange","blue","green");
		#****NOTE: as the order of the samples is not known, using colors could be confusing because
		# unrelated samples could have the same color by chance. So, better to use black for all.

		my $colorArray="col=c(".$cuffnormColors.")";

		print RSCRIPT "plot(PCAcomponent[,1],PCAcomponent[,2],ylab=\"PC2\",xlab=\"PC1\",main=c(\"PCA\"),pch=16,";
		print RSCRIPT $colorArray.")\n"; #",xlim=c(-1,1), ylim=c(-1,1))\n";
		print RSCRIPT "abline(h=0,v=mean(PCAcomponent[,1],h=0,col=\"grey\"))\n";
		print RSCRIPT "text(PCAcomponent[,1],PCAcomponent[,2],label=rownames(PCAcomponent),font=0,";
		print RSCRIPT $colorArray;
		print RSCRIPT ",cex=0.75,pos=ifelse(PCAcomponent[,1]<mean(PCAcomponent[,1]),yes=4,no=2))\n";
		print RSCRIPT "dev.off()\n";
		print RSCRIPT "write.table(PCAcomponent,file=\"".$PCAcompFile."\",sep=\"\\t\",col.names=NA)\n";	

		#3D PCA, there are several files with several different views
		my $PCAfile3D=$prefix.".samplesPCA3D_60.pdf";	
		print RSCRIPT "# 3D PCA\n";
		print RSCRIPT "library(\"scatterplot3d\")\n";
		print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
		print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=60)\n";
		print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
		print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
		print RSCRIPT "dev.off()\n";		

		$PCAfile3D=$prefix.".samplesPCA3D_minus60.pdf";	
		print RSCRIPT "# 3D PCA\n";
		print RSCRIPT "library(\"scatterplot3d\")\n";
		print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
		print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=-60)\n";
		print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
		print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
		print RSCRIPT "dev.off()\n";	

		close(RSCRIPT);

		#executes R
		my $command="export TMPDIR=".$executionCreatedTempDir."; R --vanilla < ".$Rscriptfile;
		print "\n\t[executing] ".$command."\n";		
		system($command);		
		

	}#foreach my $file(@files)	
	
}

1
