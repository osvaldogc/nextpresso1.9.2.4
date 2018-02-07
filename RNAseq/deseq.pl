#!/usr/bin/perl -w

# Author: Osvaldo Grana
# Description: runs DESeq
# v0.1		apr2015
# v0.2		dic2016 - adds cumulative variance for the samples PCAs
#			- creates rnk files for GSEA for all comparisons
#			- creates a normalized count table with all the samples
# v0.3		oct2017 - (1) Bug solved in PCA construction
#			  (2) Filters out MIR, SNORD, SNORA, and SCARNA transcripts for the RNK file
# v1.9.2	ene2018 - added scale=TRUE in prcomp function used to build PCA
#			  added new output file with PCA component values

use strict;
use FindBin qw($Bin); #finds out script path
use File::Basename qw(dirname); #calls dirname function to find out the parent dir below
use File::Spec::Functions qw(catdir); #calls catdir function
use warnings;
use File::Path qw(make_path remove_tree);

#loads own packages
use lib catdir(dirname($Bin), 'Utils');#finds out Utils dir from parent dir
use Config; # to check if perl was compiled with thread support
use Getopt::Long; #to get options

sub main();
sub runDeseq($$$$$$$$$$);
sub createFullTable($$$$);
sub executeDESeqANDcomputeCorrelationAndPCA($$$$$$$$);
sub createExcel($$$$);
sub createGSEArnkFile($$);
sub help();

main();



sub main(){
	#Checks if perl was compiled with thread support
	$Config{useithreads} or die("\n\n**** Please recompile Perl with thread support before running this program.\n\n");
	
	# Verbose stack traces
	$SIG{__DIE__} =  \&confess;
	$SIG{__WARN__} = \&confess;


	my($minimalNormalizedReadCountsPerSample,$deseqOutDir,$comparisonName,$deseqInputFiles,$deseqComparisonLabels,$extraPathsRequired,$GTF,$tmpDir,$nThreads);
	my($alpha,$pAdjustMethod);
	
	undef($alpha);
	undef($pAdjustMethod);
	undef($deseqOutDir);
	undef($comparisonName);
	undef($deseqInputFiles);
	undef($deseqComparisonLabels);
	undef($extraPathsRequired);
	undef($GTF);
	undef($tmpDir);

	GetOptions(
		"nThreads=i"=>\$nThreads,
		"alpha=f"=>\$alpha,
		"pAdjustMethod=s"=>\$pAdjustMethod,
		"deseqOutDir=s"=>\$deseqOutDir, #string
		"comparisonName=s"=>\$comparisonName,	#string
		"deseqInputFiles=s"=>\$deseqInputFiles,	#string
		"deseqComparisonLabels=s"=>\$deseqComparisonLabels,	#string
		"extraPathsRequired=s"=>\$extraPathsRequired,	#string
		"GTF=s"=>\$GTF,	#string
		"tmpDir=s"=>\$tmpDir,	#string
	);

#	if(!defined($mode) || !defined($tophatPath) || !defined($bowtiePath) || !defined($samtoolsPath) || !defined($referenceSequence) || !defined($indexPrefixForReferenceSequence) || !defined($samples) || !defined($outputFileNames) || !defined($GTF))
#	{
#		help();
#	} 

					
	runDeseq($nThreads,$alpha,$pAdjustMethod,$deseqOutDir,$comparisonName,$deseqInputFiles,$deseqComparisonLabels,$extraPathsRequired,$GTF,$tmpDir);
	
	
		
}

sub runDeseq($$$$$$$$$$){
	my($nThreads,$alpha,$pAdjustMethod,$deseqOutDir,$comparisonName,$deseqInputFiles,$deseqComparisonLabels,$extraPathsRequired,$GTF,$tmpDir)=@_;

	createFullTable($deseqInputFiles,$tmpDir,$comparisonName,$deseqOutDir);
	executeDESeqANDcomputeCorrelationAndPCA($nThreads,$deseqComparisonLabels,$deseqInputFiles,$tmpDir,$comparisonName,$deseqOutDir,$alpha,$pAdjustMethod);
	createGSEArnkFile($comparisonName,$deseqOutDir);
				
}

sub createFullTable($$$$){
	my($deseqInputFiles,$tmpDir,$comparisonName,$deseqOutDir)=@_;
	
	$deseqInputFiles=~ s/\:/,/g;
	
	my @files=split(',',$deseqInputFiles);
	my $firstSample=1;
	my $header="Name";
	my $commandToJoinAllStuffInOneTable="paste -d '\t' ";
	my @listOfSamplesToCheckIFwellSorted=();
	
	foreach my $sample(@files){
		print "[deseq] processing ".$sample."\n";
		chomp($sample);
		
		my $sortedSample=$sample.".sorted";
		#deletes undesired stuff from original files
		
		my @tokens=split('/',$sample);
		my $sampleName=$tokens[@tokens-1];
		$sampleName=~ s/\.xls//;
		system("LC_ALL=C; export LC_ALL; grep -v 'no_feature' ".$sample." | grep -v 'ambiguous' | grep -v 'too_low_aQual' | grep -v 'not_aligned' | grep -v 'alignment_not_unique' | sort -k1,1 > ".$sortedSample);
		print "[sorting] LC_ALL=C; export LC_ALL; grep -v 'no_feature' ".$sample." | grep -v 'ambiguous' | grep -v 'too_low_aQual' | grep -v 'not_aligned' | grep -v 'alignment_not_unique' | sort -k1,1 > ".$sortedSample."\n";

		if($firstSample){
			#gets row names (transcript names)
			system("cut -f 1 ".$sortedSample." > ".$tmpDir."/transcriptNames_".$comparisonName);
			$commandToJoinAllStuffInOneTable.=$tmpDir."/transcriptNames_".$comparisonName." ";
			$firstSample=0;
		}
		
		#column names
		$header.="\t".$sampleName;
		
		push(@listOfSamplesToCheckIFwellSorted,$sortedSample);
		
		system("cut -f 2 ".$sortedSample." > ".$tmpDir."/".$sampleName."_".$comparisonName);
		$commandToJoinAllStuffInOneTable.=$tmpDir."/".$sampleName."_".$comparisonName." ";
	}
	
	#joins partial files to create the full table
	$commandToJoinAllStuffInOneTable.="> ".$deseqOutDir.$comparisonName.".txt";
	print "[deseq] joining: ".$commandToJoinAllStuffInOneTable."\n";
	system($commandToJoinAllStuffInOneTable);
	
	print "[deseq] adding header: echo '".$header."' > ".$deseqOutDir.$comparisonName.".txt_2\n";
	system("echo '".$header."' > ".$deseqOutDir.$comparisonName.".txt_2");
	
	print "[deseq] adding body: cat ".$deseqOutDir.$comparisonName.".txt >> ".$deseqOutDir.$comparisonName.".txt_2\n";		
	system("cat ".$deseqOutDir.$comparisonName.".txt >> ".$deseqOutDir.$comparisonName.".txt_2");
	
	print "[deseq] moving: mv ".$deseqOutDir.$comparisonName.".txt_2 ".$deseqOutDir.$comparisonName.".txt\n";
	system("mv ".$deseqOutDir.$comparisonName.".txt_2 ".$deseqOutDir.$comparisonName.".txt");


	
	print "[deseq] Checking if all the samples are equally sorted\n";
	my $auxSuffix=".toCheckSorting";
	
	#creates files to check
	foreach my $currentSample(@listOfSamplesToCheckIFwellSorted){		
		system("cut -f 1 ".$currentSample." > ".$currentSample.$auxSuffix);
		print "[deseq] cut -f 1 ".$currentSample." > ".$currentSample.$auxSuffix."\n";
	}
	#compares the files to check if they are identical
	my $referenceSampleForChecking=$listOfSamplesToCheckIFwellSorted[0].$auxSuffix;	
	my $fileWithSortingCheckResults=$deseqOutDir."htseqCount_sorting_check.txt";
	system("rm -f ".$fileWithSortingCheckResults);
	
	foreach my $currentSample(@listOfSamplesToCheckIFwellSorted){		
		system("diff -qs ".$referenceSampleForChecking." ".$currentSample.$auxSuffix." >> ".$fileWithSortingCheckResults);
		print "[deseq] diff -qs ".$referenceSampleForChecking." ".$currentSample.$auxSuffix." >> ".$fileWithSortingCheckResults."\n";
	}
	#removes intermediate files
	foreach my $currentSample(@listOfSamplesToCheckIFwellSorted){		
		system("rm -f ".$currentSample.$auxSuffix);
		print "[deseq] rm -f ".$currentSample.$auxSuffix."\n";
	}
	#if files are not identical, the execution is stopped	
	my $res="";
	open(IN,"grep 'differ' ".$fileWithSortingCheckResults." |");
	$res=<IN>;
	close(IN);
	if(defined($res)){
		chomp($res);
		if($res=~ /differ/){
			print "[deseq] sorting of samples differ among them. Check the htseqcount sorted files please.\n";
			print "[execution finished]\n\n";
			die "[deseq] sorting of samples differ among them. Check the htseqcount sorted files please.\n[execution finished]\n\n";		
		}
	}else{		
		print "[deseq] all samples were equally sorted: ok\n";
	}
}

sub executeDESeqANDcomputeCorrelationAndPCA($$$$$$$$){
	#Executes DESeq for a particular comparison
	my($nThreads,$deseqComparisonLabels,$deseqInputFiles,$tmpDir,$comparisonName,$deseqOutDir,$alpha,$pAdjustMethod)=@_;
	
	my $controlName=(split('\,',$deseqComparisonLabels))[0];
	my $treatmentName=(split('\,',$deseqComparisonLabels))[1];
	
	my $controlSamples=(split('\:',$deseqInputFiles))[0];
	my $treatmentSamples=(split('\:',$deseqInputFiles))[1];
	
	my $deseqScript="";
	
	if($nThreads>1){
		$deseqScript.="library(\"BiocParallel\")\n";
		$deseqScript.="register(MulticoreParam(".$nThreads."))\n";
	}
	
	$deseqScript.="library(\"DESeq2\")\n";
	$deseqScript.="\ncountTable = read.table(\"".$deseqOutDir.$comparisonName.".txt\",header=TRUE, row.names=1)\n";
	
	my $condition="condition = factor(c(";
	my $libType="libType = c(";
	my $nSamples=0;
	my $nControlSamples=0;
	
	#repeats the same control id for every control sample
	my @tokens=split('\,',$controlSamples);
	foreach my $sample(@tokens){
		$condition.="\"".$controlName."\",";
		$libType.="\"oneFactor\",";
		$nSamples++;
		$nControlSamples++;
	}
	
	#repeats the same treatment id for every treatment sample
	@tokens=split('\,',$treatmentSamples);
	foreach my $sample(@tokens){
		$condition.="\"".$treatmentName."\",";
		$libType.="\"oneFactor\",";
		$nSamples++;
	}
	
	$condition=~ s/\,$//;
	$libType=~ s/\,$//;
	$condition.="))\n";
	$libType.=")\n";
	$deseqScript.=$condition;
	$deseqScript.=$libType;
	$deseqScript.="condition <- relevel(condition, \"".$controlName."\")\n";
	$deseqScript.="experiment_design=data.frame(\nrow.names = colnames(countTable),\ncondition,\nlibType)\n";	
	$deseqScript.="cds <- DESeqDataSetFromMatrix(countData = countTable, colData=experiment_design, design=~condition)\n";
	if($nThreads>1){
		$deseqScript.="cds_DESeqED <- DESeq(cds,parallel = TRUE)\n";
		$deseqScript.="res <- results(cds_DESeqED,parallel = TRUE,alpha = ".$alpha.", pAdjustMethod = \"".$pAdjustMethod."\")\n";
	}else{
		$deseqScript.="cds_DESeqED <- DESeq(cds)\n";
		$deseqScript.="res <- results(cds_DESeqED,alpha = ".$alpha.", pAdjustMethod = \"".$pAdjustMethod."\")\n";	
	}
	
	my $diffExp_file=$deseqOutDir.$comparisonName.".differentialExpression.txt";
	$deseqScript.="write.table(res,file = \"".$diffExp_file."\",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = \"\\t\",eol = \"\\n\", na = \"NA\", dec = \".\")\n";
	
	my $normalizedCounts_file=$deseqOutDir.$comparisonName.".normalizedCounts.xls";
	$deseqScript.="normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)\n";
	$deseqScript.="write.table(normalizedReadCounts,file = \"".$normalizedCounts_file."\",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = \"\\t\",eol = \"\\n\", na = \"NA\", dec = \".\")\n";
	
	my $scriptName=$deseqOutDir.$comparisonName.".R";
	open(OUT,">",$scriptName);
	print OUT $deseqScript;
	close(OUT);
	
	my $logFile=$deseqOutDir.$comparisonName.".Rlog";
	my $command="export TMPDIR=".$tmpDir."; R --vanilla < ".$scriptName." > ".$logFile." 2>&1";
	
	print "[executing] ".$command."\n";
	system($command);

	# adds header to first column in diff exp and normalized files
	print "[deseq]adding header: echo -n 'Name' > ".$diffExp_file."_2\n";
	system("echo -n 'Name' > ".$diffExp_file."_2");
	
	print "[deseq] adding body: cat ".$diffExp_file." >> ".$diffExp_file."_2\n";		
	system("cat ".$diffExp_file." >> ".$diffExp_file."_2");
	
	print "[deseq] moving: mv ".$diffExp_file."_2 ".$diffExp_file."\n";
	system("mv ".$diffExp_file."_2 ".$diffExp_file);

	# adds header for first column
	print "[deseq]adding header: echo -n 'Name' > ".$normalizedCounts_file."_2\n";
	system("echo -n 'Name' > ".$normalizedCounts_file."_2");
	
	print "[deseq] adding body: cat ".$normalizedCounts_file." >> ".$normalizedCounts_file."_2\n";		
	system("cat ".$normalizedCounts_file." >> ".$normalizedCounts_file."_2");
	
	print "[deseq] moving: mv ".$normalizedCounts_file."_2 ".$normalizedCounts_file."\n";
	system("mv ".$normalizedCounts_file."_2 ".$normalizedCounts_file);
	
	
	#creates EXCEL file
	createExcel($diffExp_file,$treatmentName,$alpha,$tmpDir);	
	
	
	#correlation among samples and PCA
	
	#first creates R script
	my $corrScriptFile=$deseqOutDir.$comparisonName.".correlationAndPCA.R";
	my $corrFile=$deseqOutDir.$comparisonName.".pearsonCorrelationAmongSamples.xls";
	my $variance=$deseqOutDir.$comparisonName.".proportionOfVariance.pdf";
	my $PCAfile=$deseqOutDir.$comparisonName.".samplesPCA.pdf";
	my $cumVarFile=$deseqOutDir.$comparisonName.".cumulativeVariance.xls";
	my $PCAcompFile=$deseqOutDir.$comparisonName.".PCAcomponent_values.xls";
	
	open(RSCRIPT,">",$corrScriptFile);	
	print RSCRIPT "normalizedCounts=read.table(\"".$normalizedCounts_file."\",header=T,sep = \"\\t\")\n";
	print RSCRIPT "dim(normalizedCounts)\n";
	print RSCRIPT "normMatrix <- as.matrix(normalizedCounts[,2:".($nSamples+1)."])\n";	
	#pearsonCorrelation
	print RSCRIPT "samplesCorr <- cor(normMatrix,method=\"pearson\")\n";
	print RSCRIPT "write.table(samplesCorr,file=\"".$corrFile."\",sep=\"\\t\",col.names=NA,eol = \"\\n\",dec = ".")\n";
	#print RSCRIPT "PC <- prcomp(t(normMatrix),scale=TRUE)\n";
	print RSCRIPT "possibleError<-tryCatch(PC <- prcomp(t(normMatrix),scale=TRUE),error=function(e) e)\n";
	print RSCRIPT "if(inherits(possibleError, \"error\")){PC <- prcomp(t(normMatrix),scale=FALSE)}\n";
	#proportion of variance
	print RSCRIPT "pdf(\"".$variance."\")\n";
	print RSCRIPT "screeplot(PC)\n";
	print RSCRIPT "dev.off()\n";
	#Cumulative variance
	print RSCRIPT "eig <- (PC\$sdev)^2\n";
	print RSCRIPT "variance <- eig*100/sum(eig)\n";
	print RSCRIPT "cumvar <- cumsum(variance)\n";
	print RSCRIPT "cumvar_matrix <- matrix (cumvar,nrow=length(cumvar),ncol=2,byrow=FALSE)\n";
	print RSCRIPT "colnames(cumvar_matrix) <- c(\"PC\",\"% cumulative variance\")\n";
	print RSCRIPT "for(i in 1:length(cumvar)){ cumvar_matrix[i][1]=i }\n";	
	print RSCRIPT "write.table(cumvar_matrix,file=\"".$cumVarFile."\",sep=\"\\t\",col.names=TRUE,row.names=FALSE)\n";	
	#PCA graph
	print RSCRIPT "pdf(\"".$PCAfile."\")\n";
	print RSCRIPT "PCAcomponent <- PC\$x\n";
	
	#defines colors
	my @colors=("lightgreen","lightpink","lightblue","lightyellow","ligthred",
		"brightgreen","brightpink","brightblue","brightyellow","ligthred",
		"black","red",",green","yellow","blue","magenta","cyan");
	
	my $colorArray="col=c(";
	my $counter=0;
	my $checked=0;
	for(my $h=0;$h<$nSamples;$h++){
		if(!$checked && $h>$nControlSamples-1){$counter++; $checked=1}
		if($comparisonName ne "ALLsamples"){$colorArray.="\"".$colors[$counter]."\","}
		else{$colorArray.="\"black\","}
	}
	$colorArray=~ s/\,$//;
	$colorArray.=")";
		
	print RSCRIPT "plot(PCAcomponent[,1],PCAcomponent[,2],ylab=\"PC2\",xlab=\"PC1\",main=c(\"PCA\"),pch=16,";
	print RSCRIPT $colorArray.")\n"; #",xlim=c(-1,1), ylim=c(-1,1))\n";
	print RSCRIPT "abline(h=0,v=mean(PCAcomponent[,1],h=0,col=\"grey\"))\n";
	print RSCRIPT "text(PCAcomponent[,1],PCAcomponent[,2],label=rownames(PCAcomponent),font=0,";
	print RSCRIPT $colorArray;
	print RSCRIPT ",cex=0.75,pos=ifelse(PCAcomponent[,1]<mean(PCAcomponent[,1]),yes=4,no=2))\n";
	print RSCRIPT "dev.off()\n";
	print RSCRIPT "write.table(PCAcomponent,file=\"".$PCAcompFile."\",sep=\"\\t\",col.names=NA)\n";	

	#3D PCA, there are several files with several different views
	my $PCAfile3D=$deseqOutDir.$comparisonName.".samplesPCA3D_60.pdf";	
	print RSCRIPT "# 3D PCA\n";
	print RSCRIPT "library(\"scatterplot3d\")\n";
	print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
	print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=60)\n";
	print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
	print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
	print RSCRIPT "dev.off()\n";		

	$PCAfile3D=$deseqOutDir.$comparisonName.".samplesPCA3D_minus60.pdf";	
	print RSCRIPT "# 3D PCA\n";
	print RSCRIPT "library(\"scatterplot3d\")\n";
	print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
	print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=-60)\n";
	print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
	print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
	print RSCRIPT "dev.off()\n";	
	
	close(RSCRIPT);
	
	#executes R
	$logFile=$deseqOutDir.$comparisonName.".correlationAndPCA.Rlog";
	$command="export TMPDIR=".$tmpDir."; R --vanilla < ".$corrScriptFile." > ".$logFile." 2>&1";
	print "[executing] ".$command."\n";		
	system($command);
	
}

sub createExcel($$$$){
	my($file,$treatmentName,$alpha,$tmpDir)=@_;
	
	use Excel::Writer::XLSX;	
	
	my $inputFile=$file;
	my $outFile=$inputFile;
	$outFile=~ s/txt$/xlsx/;


	#the input files are sorted by padj before doing the work,
	#so the output files are already sorted by q-value ascending
	my $inputFileSORTED=$inputFile.".SORTED";
	system("head -n 1 ".$inputFile." > ".$inputFileSORTED);
	#LC_ALL=C -> forces the local configuration to use '.' to represent decimal numbers
	#***** filters out lines with 'NA' values
	system("LC_ALL=C; export LC_ALL; grep -v 'baseMean'  ".$inputFile." | grep -v '	NA'  | sort -k7,7g >> ".$inputFileSORTED);
	
	#opens gene_exp file
	if(open(IN,$inputFileSORTED)){

		#creates a new Excel workbook
    		my $workbook=Excel::Writer::XLSX->new($outFile);
		$workbook->set_tempdir($tmpDir);
				
    		#adds a worksheet
    		my $worksheet = $workbook->add_worksheet();		
		
    		#add and define a format
    		my $format;
		my $row = 0;
		
		#legend Up/Down
		$format= $workbook->add_format();
		$format->set_color( 'red' );
		$format->set_bold(0);  # Turn bold on
		$format->set_bg_color( 'yellow' );
		$worksheet->merge_range($row,0,$row,2,"Upregulated in ".$treatmentName, $format );
		$row++;
		
		$format= $workbook->add_format();
		$format->set_color( 'green' );
		$format->set_bold(0);  # Turn bold on
		$format->set_bg_color( 'yellow' );
		$worksheet->merge_range($row,0,$row,2,"Downregulated in ".$treatmentName, $format );
		$row++;		
		
		$format= $workbook->add_format();
		$format->set_color( 'black' );
		$format->set_bold(1);  # Turn bold on
		$format->set_bg_color(-1);
		$worksheet->merge_range($row,0,$row,4,"FDR=".$alpha, $format );
		$row++;
		$row++;
		$row++;
		
		my $memoryFirstRow=0;
		my $memoryLastRow=0;
		my $memoryFirstCol=0;
		my $memoryLastCol=0;	
			
		foreach my $line(<IN>){
			chomp($line);
			my @tokens=split('\t+',$line);

			#header line
			if($line=~ /baseMean/){
				$format= $workbook->add_format();
				$format->set_color( 'black' );
				$format->set_bold(1);  # Turn bold on
				$memoryFirstRow=$row;
				$memoryLastRow=$row;
				$memoryFirstCol=0;
				$memoryLastCol=@tokens-1;
				
			}else{#data line

				$format= $workbook->add_format();
				$format->set_bold(0);  # Turn bold on

				#if log2FoldChange>=0, then red font
				if($tokens[2]!~ /^\-/){
					$format->set_color( 'red' );
				}else{#green font
					$format->set_color( 'green' );
				}

				#yellow background if significant difference
				if($tokens[6]<=$alpha){
					$format->set_bg_color( 'yellow' );
				}		
				
			}

			
			#prints line with data			
			for(my $col=0;$col<@tokens;$col++){
				$worksheet->write( $row, $col,$tokens[$col], $format );
			}

			$row++;

		}
		close(IN);
				
		#adds autofilter		
		$worksheet->autofilter($memoryFirstRow,$memoryFirstCol,$memoryLastRow,$memoryLastCol);
		
		#set columns width
		$worksheet->set_column('A:A' , 10 );
		$worksheet->set_column('B:B' , 15 );
		$worksheet->set_column('C:C' , 20 );
		$worksheet->set_column('D:D' , 10);
		$worksheet->set_column('E:G' ,17);

		$workbook->close();
		
		
		system("rm -rf ".$inputFileSORTED);
	}else{
		print STDERR "\n[ERROR createExcel]: $inputFile file does not exist\n\n";
		exit(-1);
	}

}

sub createGSEArnkFile($$){
# A RNK file is created from log2 fold-change (instead of statistic)
	
	my($comparisonName,$deseqOutDir)=@_;
	
	my $diffExp_file=$deseqOutDir.$comparisonName.".differentialExpression.txt";
	
	my $unsortedRnkFile=$deseqOutDir.$comparisonName.".rnk.unsorted";
	my $rnkFile=$deseqOutDir.$comparisonName.".rnk";

	########## GENE LEVEL ##########	

	if(open(DIFFEXP,$diffExp_file)){
	open(UNSORTED,">",$unsortedRnkFile);
	foreach my $line(<DIFFEXP>){
		chomp($line);
		
		if($line !~ /log2FoldChange/){
			my @tokens=split('\t',$line);
			
			#log2FoldChange
			my $log2FC=$tokens[2];
			
			if($log2FC ne "NA"){
				my $geneID=$tokens[0];
				
				#filters out MIR, SNORD, SNORA, and SCARNA transcripts
				if(uc($geneID)!~ /^MIR/ && uc($geneID)!~ /^SNORD/ && uc($geneID)!~ /^SNORA/ && uc($geneID)!~ /^SCARNA/){
					print UNSORTED uc($geneID)."\t".$log2FC."\n";
				}
			}		
		}
	}	
	close(DIFFEXP);
	close(UNSORTED);
	

	#sorts the auxiliar file, saves the sorted file in the output file and deletes the auxiliar file
	#LC_ALL=C -> forces the local configuration to use '.' to represent decimal numbers
	#--stable avoids default (unrequested) sorting by gene name when log2foldchange=0
	system("LC_ALL=C; export LC_ALL; sort -k2g --stable ".$unsortedRnkFile." -o ".$rnkFile."; rm -f ".$unsortedRnkFile);
	print "\t[executing] LC_ALL=C; export LC_ALL; sort -k2g --stable ".$unsortedRnkFile." -o ".$rnkFile."; rm -f ".$unsortedRnkFile."\n";
	
	
	}
	else{
		print STDERR "\n[ERROR createGSEArnkFiles_from_DESeq2]: ".$diffExp_file." file does not exist\n\n";
		exit(-1);		
	}
}


sub help(){
	my $usage = qq{
--quality.pl--
			
	Execution error. Check out the execution modes in the following examples:
			
	Example for trimming a single end experiment
		perl preprocess.pl --nThreads 6 --mode 1 --seqtkPath /path --positions 0:3,0:3 --samples /path/file1.fq:/path/file1.trimmed.fq
				
	Example for trimming a paired end experiment
		perl preprocess.pl --nThreads 6 --mode 1 --seqtkPath /path --positions 0:3,0:3 --samples /path/file1.fq:/path/file1.trimmed.fq,/path/file2.fq:/path/file2.trimmed.fq
				
	Example for downsampling a single end experiment
		perl preprocess.pl --nThreads 6 --mode 2 --seqtkPath /path --seed 123L --nReads 10000 --samples /path/file1.fq,/path/file2.fq
				
	Example for downsampling a paired end experiment
		perl preprocess.pl --nThreads 6 --mode 2 --seqtkPath /path --seed 123L --nReads 10000 --samples /path/file1_R1.fq:/path/file1_R2.fq,/path/file2_R1.fq:/path/file2_R2.fq
				


};
 
	print STDERR $usage;
	exit(1);
}

