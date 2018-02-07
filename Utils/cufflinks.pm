#!/usr/bin/perl

# Author : Osvaldo Grana
# Description: performs transcripts quantification and assembly
# v0.3		oct2016
# v0.4		dic2016 - adds cumulative variance for the samples
# v1.9.2	ene2018 - added scale=TRUE in prcomp function used to build PCA
#			  added new output file with PCA component values

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";


package cufflinks;

sub runCufflinks(){
	my($extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$alignmentsDir,$referenceSequence,$samples,$GTF,
		$libraryType,$nThreads,$nCufflinksThreads,$fragBiasCorrect,$multiReadCorrect,$libraryNormalizationMethod,$useGTFwithCufflinks,$maxBundleFrags,
		$noEffectiveLengthCorrection,$noLengthCorrection,$normalization,$inputFile)=@_;
	
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
	
	$command.="cufflinks -p ".$nCufflinksThreads." ";
	
	 #posibilities: unstranded, firststrand, secondstrand    
        if(lc($libraryType) eq "unstranded"){
                $command.="--library-type fr-unstranded ";
        }elsif(lc($libraryType) eq "firststrand"){
                $command.="--library-type fr-firststrand ";
        }elsif(lc($libraryType) eq "secondstrand"){
                $command.="--library-type fr-secondstrand ";
        }

#if specified cufflinks does not run, although it does not return any error	
#	$command.="--library-norm-method ".$libraryNormalizationMethod." ";

	if($useGTFwithCufflinks eq "true"){$command.="--GTF ".$GTF." "}
	
	if($multiReadCorrect eq "true"){
		$command.="--multi-read-correct ";
	}
	
	if($noEffectiveLengthCorrection eq "true"){
		$command.="--no-effective-length-correction ";
	}
	
	if($noLengthCorrection eq "true"){
		$command.="--no-length-correction ";
	}
	
	if($normalization eq "compatibleHits"){
		$command.="--compatible-hits-norm ";
	}else{
		$command.="--total-hits-norm ";	
	}
	
	if($fragBiasCorrect eq "true"){
		$command.="--frag-bias-correct ".$referenceSequence." ";
	}
	
	$command.="--max-bundle-frags ".$maxBundleFrags." ";
	
	$command.="-o ".$outDir.$inputFile." ";
	$command.=$alignmentsDir.$inputFile."/accepted_hits.bam";
	
	print "\n\t[executing] ".$command."\n";	
	system($command);
	
}

sub calculateCorrelationsAndPCA_GeneLevel(){
	my($experimentName,$outDir,$samples,$executionCreatedTempDir)=@_;
	
	my @files=split(',',$samples);	
	
	my $XLSfile=$outDir."samplesFPKMs_geneLevel.xls";
	my $GCTfile=$outDir."samplesFPKMs_geneLevel.gct";
	system("rm -f ".$XLSfile." ".$GCTfile);
	
	# FIRST generates a matrix to use with R for correlation and PCA tests
	# (this matrix does not contain gene names)
	foreach my $file(@files){
		chomp($file);
	
		# gene level
		my $inputFile=$outDir.$file."/genes.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";		
			
		#sorts the genes/transcripts
		my $command="LC_ALL=C; export LC_ALL; sort -k1n ".$inputFile." -o ".$sortedFile;
		print "\n\t[executing] ".$command."\n";		
		system($command);
		
		#creates header files and body files, and joins both in '.final' file
		#**** grep -E 'chr[12]?[0-9XYxyMm]:' to avoid lines like the following with chr6_apd_hap1:174179-195170
		#$command="grep -v 'FPKM' ".$sortedFile." | grep -E 'chr[12]?[0-9XYxyMm]:' > ".$sortedFile.".fpkm; grep 'FPKM' ".$sortedFile." > ".$sortedFile.".header; cat ".$sortedFile.".header ".$sortedFile.".fpkm > ".$sortedFile.".auxiliar";
		
		#grep -E deleted: it causes problems when the reference genomes are not 'chr' based, like when having scaffolds for example
		$command="grep -v 'FPKM' ".$sortedFile." > ".$sortedFile.".fpkm; grep 'FPKM' ".$sortedFile." > ".$sortedFile.".header; cat ".$sortedFile.".header ".$sortedFile.".fpkm > ".$sortedFile.".auxiliar";
		print "\n\t[executing] ".$command."\n";		
		system($command);		
		
		#creates the final file
		$command="awk '{if(\$10==\"FPKM\"){print \"".$file."\"} else { print \$10}}' ".$sortedFile.".auxiliar > ".$sortedFile.".final";
		print "\n\t[executing] ".$command."\n";		
		system($command);

		# captures gene names
		$command="echo 'Name' > ".$sortedFile.".geneNames";
		print "\n\t[executing] ".$command."\n";		
		system($command);
		
		#**** grep -E 'chr[12]?[0-9XYxyMm]:' to avoid lines like the following with chr6_apd_hap1:174179-195170
		#$command="grep -v 'tracking_id' ".$sortedFile." | grep -E 'chr[12]?[0-9XYxyMm]:' | cut -f 1 >> ".$sortedFile.".geneNames";
		
		#grep -E deleted: it causes problems when the reference genomes are not 'chr' based, like when having scaffolds for example
		$command="grep -v 'tracking_id' ".$sortedFile." | cut -f 1 >> ".$sortedFile.".geneNames";		
		print "\n\t[executing] ".$command."\n";		
		system($command);			
	}

	#creates a FPKM matrix for gene-e in gct file format (with gene names)
	#(the gene names are taken from the first sample 'gene names' derived file '$partialFiles')
	my $partialFiles=$outDir.$files[0]."/genes.fpkm_tracking_sorted.txt.geneNames ";
	foreach my $file(@files){
		chomp($file);
		my $inputFile=$outDir.$file."/genes.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";
		
		$partialFiles.=$sortedFile.".final ";
	}
	
	#pastes everything in a common matrix	
	my $command="paste -d \"\t\" ".$partialFiles." > ".$XLSfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);


	# SECOND generates a matrix to use with gene-e
	# (this matrix contains gene names)
	foreach my $file(@files){
		chomp($file);
	
		# gene level
		my $inputFile=$outDir.$file."/genes.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";
		
		# captures gene names
		my $command="echo 'Name	Description' > ".$sortedFile.".geneNames";
		print "\n\t[executing] ".$command."\n";		
		system($command);
		
		#**** grep -E 'chr[12]?[0-9XYxyMm]:' to avoid lines like the following with chr6_apd_hap1:174179-195170
		#$command="grep -v 'tracking_id' ".$sortedFile." | grep -E 'chr[12]?[0-9XYxyMm]:' | cut -f 1 | awk '{OFS=\"\t\"}{print \$1,\"No description\"}' >> ".$sortedFile.".geneNames";		
		
		#grep -E deleted: it causes problems when the reference genomes are not 'chr' based, like when having scaffolds for example
		$command="grep -v 'tracking_id' ".$sortedFile." | cut -f 1 | awk '{OFS=\"\t\"}{print \$1,\"No description\"}' >> ".$sortedFile.".geneNames";				
		print "\n\t[executing] ".$command."\n";		
		system($command);
	}
	
	#creates a FPKM matrix for gene-e in gct file format (with gene names)
	#(the gene names are taken from the first sample 'gene names' derived file '$partialFiles')
	$partialFiles=$outDir.$files[0]."/genes.fpkm_tracking_sorted.txt.geneNames ";
	foreach my $file(@files){
		chomp($file);
		my $inputFile=$outDir.$file."/genes.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";
		
		$partialFiles.=$sortedFile.".final ";
	}
	
	#pastes everything in a common matrix	
	$command="paste -d \"\t\" ".$partialFiles." > ".$GCTfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);
	
	#finally, gct format is added
	#(using auxiliar files called "a" and "b")
	my $a=$executionCreatedTempDir."/cufflinks_genelevel_a_".$experimentName;
	my $b=$executionCreatedTempDir."/cufflinks_genelevel_b_".$experimentName;

	$command="rm -f ".$b." ".$a;
	system($command);
		
	open(AUXILIAR,">",$a);
	print AUXILIAR "#1.2\n";
	#counts the number of genes in the gct file
	open(NLINES,"wc -l ".$XLSfile." |");
	my $nlines=<NLINES>;
	close(NLINES);	
	$nlines=(split(' ',$nlines))[0];
	# prints gct format
	print AUXILIAR ($nlines-1)."\t".@files."\n";
	close(AUXILIAR);
	
	$command="mv -f ".$GCTfile." ".$b."; cat ".$a." ".$b." > ".$GCTfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);	
	
	#finally, deletes all intermediate files
	foreach my $file(@files){
		chomp($file);

		# gene level
		my $inputFile=$outDir.$file."/genes.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";

		$command="rm -f ".$sortedFile." ".$sortedFile.".fpkm ".$sortedFile.".header ".$sortedFile.".auxiliar ".$sortedFile.".final ".$sortedFile.".geneNames";
		print "\n\t[executing] ".$command."\n";
		system($command);			
	}
	
	# CORRELATION TABLE AND PCA WITH R
	
	#first creates R script
	my $Rscriptfile=$outDir."correlationAndPCA_geneLevel.R";
	my $corrFile=$outDir."pearsonCorrelationsAmongSamples_geneLevel.xls";
	my $variance=$outDir."proportionOfVariance_geneLevel.pdf";
	my $PCAfile=$outDir."samplesPCA_geneLevel.pdf";
	my $cumVarFile=$outDir."cumulativeVariance_geneLevel.xls";
	my $PCAcompFile=$outDir.".PCAcomponent_values_geneLevel.xls";
	
	open(RSCRIPT,">",$Rscriptfile);	
	print RSCRIPT "FPKMs=read.table(\"".$XLSfile."\",header=T)\n";
	print RSCRIPT "dim(FPKMs)\n";
	print RSCRIPT "FPKMmatrix <- as.matrix(FPKMs[,2:".(@files+1)."])\n";	
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
	my @colors=("black","black","black","black","black","black");
	
	my $colorArray="col=c(";
	my $counter=0;
	for(my $h=0;$h<@files;$h++){
		$colorArray.="\"".$colors[$counter]."\",";
		$counter++;
		if($counter>=(@colors-1)){$counter=0}		
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
	my $PCAfile3D=$outDir."samplesPCA3D_60_geneLevel.pdf";	
	print RSCRIPT "# 3D PCA\n";
	print RSCRIPT "library(\"scatterplot3d\")\n";
	print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
	print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=60)\n";
	print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
	print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
	print RSCRIPT "dev.off()\n";		

	$PCAfile3D=$outDir."samplesPCA3D_minus60_geneLevel.pdf";	
	print RSCRIPT "# 3D PCA\n";
	print RSCRIPT "library(\"scatterplot3d\")\n";
	print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
	print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=-60)\n";
	print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
	print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
	print RSCRIPT "dev.off()\n";	
	
	close(RSCRIPT);
	
	#executes R
	$command="export TMPDIR=".$executionCreatedTempDir."; R --vanilla < ".$Rscriptfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);	
}

sub calculateCorrelationsAndPCA_IsoformLevel(){
	my($experimentName,$outDir,$samples,$executionCreatedTempDir)=@_;
	
	my @files=split(',',$samples);	
	
	my $XLSfile=$outDir."samplesFPKMs_isoformLevel.xls";
	my $GCTfile=$outDir."samplesFPKMs_isoformLevel.gct";
	system("rm -f ".$XLSfile." ".$GCTfile);
	
	# FIRST generates a matrix to use with R for correlation and PCA tests
	# (this matrix does not contain gene names)
	foreach my $file(@files){
		chomp($file);
	
		# gene level
		my $inputFile=$outDir.$file."/isoforms.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";		
			
		#sorts the genes/transcripts
		my $command="LC_ALL=C; export LC_ALL; sort -k1n ".$inputFile." -o ".$sortedFile;
		print "\n\t[executing] ".$command."\n";		
		system($command);
		
		#creates header files and body files, and joins both in '.final' file
		#**** grep -E 'chr[12]?[0-9XYxyMm]:' to avoid lines like the following with chr6_apd_hap1:174179-195170
		#$command="grep -v 'FPKM' ".$sortedFile." | grep -E 'chr[12]?[0-9XYxyMm]:' > ".$sortedFile.".fpkm; grep 'FPKM' ".$sortedFile." > ".$sortedFile.".header; cat ".$sortedFile.".header ".$sortedFile.".fpkm > ".$sortedFile.".auxiliar";
		
		#grep -E deleted: it causes problems when the reference genomes are not 'chr' based, like when having scaffolds for example
		$command="grep -v 'FPKM' ".$sortedFile." > ".$sortedFile.".fpkm; grep 'FPKM' ".$sortedFile." > ".$sortedFile.".header; cat ".$sortedFile.".header ".$sortedFile.".fpkm > ".$sortedFile.".auxiliar";
		print "\n\t[executing] ".$command."\n";		
		system($command);		
		
		#creates the final file
		$command="awk '{if(\$10==\"FPKM\"){print \"".$file."\"} else { print \$10}}' ".$sortedFile.".auxiliar > ".$sortedFile.".final";
		print "\n\t[executing] ".$command."\n";		
		system($command);

		# captures gene names
		$command="echo 'Name' > ".$sortedFile.".geneNames";
		print "\n\t[executing] ".$command."\n";		
		system($command);
		
		#**** grep -E 'chr[12]?[0-9XYxyMm]:' to avoid lines like the following with chr6_apd_hap1:174179-195170
		#$command="grep -v 'tracking_id' ".$sortedFile." | grep -E 'chr[12]?[0-9XYxyMm]:' | cut -f 1 >> ".$sortedFile.".geneNames";
		
		#grep -E deleted: it causes problems when the reference genomes are not 'chr' based, like when having scaffolds for example
		$command="grep -v 'tracking_id' ".$sortedFile." | cut -f 1 >> ".$sortedFile.".geneNames";		
		print "\n\t[executing] ".$command."\n";		
		system($command);			
	}

	#creates a FPKM matrix for gene-e in gct file format (with gene names)
	#(the gene names are taken from the first sample 'gene names' derived file '$partialFiles')
	my $partialFiles=$outDir.$files[0]."/isoforms.fpkm_tracking_sorted.txt.geneNames ";
	foreach my $file(@files){
		chomp($file);
		my $inputFile=$outDir.$file."/isoforms.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";
		
		$partialFiles.=$sortedFile.".final ";
	}
	
	#pastes everything in a common matrix	
	my $command="paste -d \"\t\" ".$partialFiles." > ".$XLSfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);


	# SECOND generates a matrix to use with gene-e
	# (this matrix contains gene names)
	foreach my $file(@files){
		chomp($file);
	
		# gene level
		my $inputFile=$outDir.$file."/isoforms.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";
		
		# captures gene names
		my $command="echo 'Name	Description' > ".$sortedFile.".geneNames";
		print "\n\t[executing] ".$command."\n";		
		system($command);
		
		#**** grep -E 'chr[12]?[0-9XYxyMm]:' to avoid lines like the following with chr6_apd_hap1:174179-195170
		#$command="grep -v 'tracking_id' ".$sortedFile." | grep -E 'chr[12]?[0-9XYxyMm]:' | cut -f 1 | awk '{OFS=\"\t\"}{print \$1,\$4}' >> ".$sortedFile.".geneNames";		
		
		#grep -E deleted: it causes problems when the reference genomes are not 'chr' based, like when having scaffolds for example
		$command="grep -v 'tracking_id' ".$sortedFile." | cut -f 1,4 | awk '{OFS=\"\t\"}{print \$1,\$2}' >> ".$sortedFile.".geneNames";				
		print "\n\t[executing] ".$command."\n";		
		system($command);
	}
	
	#creates a FPKM matrix for gene-e in gct file format (with gene names)
	#(the gene names are taken from the first sample 'gene names' derived file '$partialFiles')
	$partialFiles=$outDir.$files[0]."/isoforms.fpkm_tracking_sorted.txt.geneNames ";
	foreach my $file(@files){
		chomp($file);
		my $inputFile=$outDir.$file."/isoforms.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";
		
		$partialFiles.=$sortedFile.".final ";
	}
	
	#pastes everything in a common matrix	
	$command="paste -d \"\t\" ".$partialFiles." > ".$GCTfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);
	
	#finally, gct format is added
	#(using auxiliar files called "a" and "b")
	my $a=$executionCreatedTempDir."/cufflinks_isoformlevel_a_".$experimentName;
	my $b=$executionCreatedTempDir."/cufflinks_isoformlevel_b_".$experimentName;

	$command="rm -f ".$b." ".$a;
	system($command);
	
	open(AUXILIAR,">",$a);
	print AUXILIAR "#1.2\n";
	#counts the number of genes in the gct file
	open(NLINES,"wc -l ".$XLSfile." |");
	my $nlines=<NLINES>;
	close(NLINES);	
	$nlines=(split(' ',$nlines))[0];
	# prints gct format
	print AUXILIAR ($nlines-1)."\t".@files."\n";
	close(AUXILIAR);
	
	$command="mv -f ".$GCTfile." ".$b."; cat ".$a." ".$b." > ".$GCTfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);	
	
	#finally, deletes all intermediate files
	foreach my $file(@files){
		chomp($file);

		# gene level
		my $inputFile=$outDir.$file."/isoforms.fpkm_tracking";
		my $sortedFile=$inputFile."_sorted.txt";

		$command="rm -f ".$sortedFile." ".$sortedFile.".fpkm ".$sortedFile.".header ".$sortedFile.".auxiliar ".$sortedFile.".final ".$sortedFile.".geneNames";
		print "\n\t[executing] ".$command."\n";
		system($command);			
	}

	# CORRELATION TABLE AND PCA WITH R
	
	#first creates R script
	my $Rscriptfile=$outDir."correlationAndPCA_isoformLevel.R";
	my $corrFile=$outDir."pearsonCorrelationsAmongSamples_isoformLevel.xls";
	my $variance=$outDir."proportionOfVariance_isoformLevel.pdf";
	my $PCAfile=$outDir."samplesPCA_isoformLevel.pdf";
	my $cumVarFile=$outDir."cumulativeVariance_isoformLevel.xls";
	my $PCAcompFile=$outDir.".PCAcomponent_values_isoformLevel.xls";
	
	open(RSCRIPT,">",$Rscriptfile);	
	print RSCRIPT "FPKMs=read.table(\"".$XLSfile."\",header=T)\n";
	print RSCRIPT "dim(FPKMs)\n";
	print RSCRIPT "FPKMmatrix <- as.matrix(FPKMs[,2:".(@files+1)."])\n";	
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
	my @colors=("black","black","black","black","black","black");
	
	my $colorArray="col=c(";
	my $counter=0;
	for(my $h=0;$h<@files;$h++){
		$colorArray.="\"".$colors[$counter]."\",";
		$counter++;
		if($counter>=(@colors-1)){$counter=0}		
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
	my $PCAfile3D=$outDir."samplesPCA3D_60_isoformLevel.pdf";	
	print RSCRIPT "# 3D PCA\n";
	print RSCRIPT "library(\"scatterplot3d\")\n";
	print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
	print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=60)\n";
	print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
	print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
	print RSCRIPT "dev.off()\n";		

	$PCAfile3D=$outDir."samplesPCA3D_minus60_isoformLevel.pdf";	
	print RSCRIPT "# 3D PCA\n";
	print RSCRIPT "library(\"scatterplot3d\")\n";
	print RSCRIPT "pdf(\"".$PCAfile3D."\")\n";
	print RSCRIPT "plot3d<-scatterplot3d(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3],bg=\"black\", ylab=\"PC2\", xlab=\"PC1\",zlab=\"PC3\",main=c(\"PCA\"),pch=20,type=\"h\",angle=-60)\n";
	print RSCRIPT "text(plot3d\$xyz.convert(PCAcomponent[,1],PCAcomponent[,2],PCAcomponent[,3]),label=rownames(PCAcomponent),font=0,";
	print RSCRIPT $colorArray.",cex=0.75,pos=ifelse(PCAcomponent[,1] < mean(PCAcomponent[,1]), yes=4, no=2))\n";
	print RSCRIPT "dev.off()\n";	
	
	close(RSCRIPT);
	
	#executes R
	$command="export TMPDIR=".$executionCreatedTempDir."; R --vanilla < ".$Rscriptfile;
	print "\n\t[executing] ".$command."\n";		
	system($command);		
}


sub runCuffmerge(){
	my($extraPathsRequired,$cufflinksPath,$samtoolsPath,$outDir,$GTF,$referenceSequence,$nCufflinksThreads,$assembliesFile,$cuffmergeOutDir)=@_;
	
	use Env qw(PATH);
	#$PATH.=":".$peakAnnotatorPath;
	#$ENV{'PATH'}.=":".$peakAnnotatorPath;
	
	# when executing cuffmerge, it first run cufflinks, searching through the path environment variable.
	# In order to identify it properly the last slash bar must be deleted, otherwise
	# cuffmerge adds one more slash bar to the path, and having to gives problems
	# something similar happens with cuffmergeOutDir
	if($cufflinksPath=~ /\/$/){
		$cufflinksPath=~ s/\/$//;
	}

	if($cuffmergeOutDir=~ /\/$/){
		$cuffmergeOutDir=~ s/\/$//;
	}
		
	my $command="export PATH=\$PATH:".$samtoolsPath.":".$cufflinksPath."; ";
	
	if($extraPathsRequired ne "NO_EXTRA_PATHS"){
		$extraPathsRequired=~ s/\'//g;
	
		my @exports=split(';',$extraPathsRequired);
		foreach my $export(@exports){
			$command.="export ".$export."; ";
		}
	}
	
	$command.="cuffmerge -p ".$nCufflinksThreads;
	$command.=" -g ".$GTF;
	$command.=" -s ".$referenceSequence;
	

	$command.=" -o ".$cuffmergeOutDir;	
	$command.=" ".$assembliesFile;
	
	print "\n\t[executing] ".$command."\n";	
	system($command);	
}

sub doSpikeInCorrection(){
## the method used here is described in:
## Lovén J, Orlando DA, Sigova AA, Lin CY, Rahl PB, Burge CB, Levens DL, Lee TI, 
##  Young RA. Revisiting global gene expression analysis. Cell. 2012 Oct
##  26;151(3):476-82. doi: 10.1016/j.cell.2012.10.012.

	my($inputFPKMmatrix,$outputFPKMmatrix,$firstPosition,$lastPosition,$dimMatrix)=@_;

	my $Rscript="data<-read.delim(\"".$inputFPKMmatrix."\")\n";
	$Rscript.="rownames(data)<-make.unique(as.character(data[,1]),sep=\"_rep\")\n";
	$Rscript.="data.m<-as.matrix(data[,2:$dimMatrix])\n";
	$Rscript.="library(\"affy\")\n";
	$Rscript.="adjusted.matrix<-loess.normalize(data.m,subset = $firstPosition:$lastPosition,epsilon = 10^-2, maxit = 1,log.it = FALSE,\nverbose = TRUE,span = 2/3,family.loess =\"symmetric\")\n";
	$Rscript.="outputMatrix <- ifelse(adjusted.matrix<0,0,adjusted.matrix)\n";
	$Rscript.="write.table(outputMatrix,file=\"".$outputFPKMmatrix."\",sep=\"\\t\",col.names=NA)\n";
	
	my $scriptName=$inputFPKMmatrix;
	$scriptName=~ s/xls/Rscript\.for\.calculating\.SpikeInCorrection\.R/;

	open(OUT,">",$scriptName);
	print OUT $Rscript;
	close(OUT);
	
	#executes R
	my $command="R --vanilla < ".$scriptName;
	print "\n\t[executing] ".$command."\n";		
	system($command);
}



1
