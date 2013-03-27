#!/usr/bin/perl

#wrapper for Sam's RNA-seq pipeline.
#Only parts of hte pipeline can be run at a time. Options are
# P: converts fastq to fasta and computes qualities stats
# C: maps, filters, and extracts un mapped reads against the genome
# M: genome mapping on
# F: genome filtering and un-mapped extraction only
# T: maps, filters, and extracts un mapped reads against the transcriptome
# V: creates hits per base and visualisation files
# R: matches to annotation and creates RPKM tables
# Rm: matches to annotation only
# Rr: creates RPKM tables only
# S: creates splicing analysis tables
# D: runs CTVRS
# U: reverses the exonerate output. Use with dUTP protocol
# B: chip-mode

##NOTES##
#Each fastq file need to be analysed in its own folder or error in depth calculation#



#inputs
use strict;
use warnings;
use Benchmark;
use Data::Dumper;
use POSIX qw(strftime);

my $scriptp="/jurg/homedir/samuel/POMBE_SEQ/analysis/SCRIPTS/";
my %status = ();
my $n = 1;
my @command=@ARGV;


my %transTime=('Jan'=>'01','Feb'=>'02','Mar'=>'03','Apr'=>'04','May'=>'05','Jun'=>'06','Jul'=>'07','Aug'=>'08','Sep'=>'09','Oct'=>'10','Nov'=>'11','Dec'=>'12');
my @today=split /\s/,(localtime);

if ($today[2] !~ /\d{2}/)
{
 $today[2]="0".$today[2];
}

my $use = "\nUse:\n-o: options (default PCTVRS, see readme.txt)\n-gp: gff path (/jurg/homedir/samuel/POMBE_SEQ/analysis/)\n-qp: fastq path (default current directory)\n-q: fastq name (NO DEFAULT)\n-g: gff name (default gff_090511.txt)\n-r: read length (default 50)\n-t: tag (default today's date (".$today[2].$transTime{$today[1]}.substr($today[4],2,2)."))\n\n";

##-o
my $options="PCTVRS";
##-gp
my $gffp = "/jurg/homedir/samuel/POMBE_SEQ/analysis/";
##-qp
my $fastqp = "./";
##-q
my $fastq;
##-g
my $gff = "gff_090511.txt";
##-r
my $read = "50";
##-t
my $tag = $today[2].$transTime{$today[1]}.substr($today[4],2,2);

my $readinput;

die $use unless (defined(@ARGV));

while (defined($readinput=shift))
{

 die $use unless ($readinput =~ /-/);


 if($readinput eq '-o')
 {
  $options=shift;
 }

 if($readinput eq '-gp')
 {
  $gffp=shift;
 }

 if($readinput eq '-qp')
 {
  $fastqp=shift;
 }

 if($readinput eq '-q')
 {
  $fastq=shift;
 }

 if($readinput eq '-g')
 {
  $gff=shift;
 }

 if($readinput eq '-r')
 {
  $read=shift;
 }

 if($readinput eq '-t')
 {
  $tag=shift;
 }
}

die "\ngff file with unexpected name structure please use ".'"'."gff_DDMMYY.txt".'"'."\n$use\n" unless ($gff =~ /gff_\d{6}\.txt/);
die "\n-q option missing with no default\n$use\n"  unless (defined($fastq));
#print "$options\n$gffp\n$fastqp\n$fastq\n$gff\n$read\n$tag\n";
#print "Let's go!\n";



if ($options eq "D")
{
 $options="CTVRS";
}

#time it
my $time0 = new Benchmark;

my $now = localtime;

print "\nCommand: $0 @command\n\ndate: $now\n\n"; 

####

die "\nTest run went OK.\n" unless ($options ne "Z");

####


#variables
my $chr_depth;
my $trs_depth;
my $depth;
my $sizeBits;
my $sizeRec;
my $depth_test;
#################QC and FASTQ processing##############################################################
if ($options =~ /P/)
{
 print "Running QC_fastqTOfastsa.pl\n"; 
 $status{$n}=system "perl ".$scriptp."QC_fastqTOfastsa.pl ".$fastqp.' '.$fastq;
 die "Problem stage $n" if $status{$n}; $n++;
}

#################ALLchr mapping and filtering#########################################################
if ($options =~ /[CMF]/)
{
 if ($options =~ /[CM]/)
 {
 print "Preparing for ALLchr mapping...\n";
 $status{$n}=system "split -d -l 2000000 ".$fastq.'.fasta';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "sh ".$scriptp."mapping.090511.ERCC.sh";
 die "Problem stage $n" if $status{$n}; $n++;

 #$status{$n}=system "wait";
 #die "Problem stage $n" if $status{$n}; $n++;

 print "Cleaning up...\n";
 $status{$n}=system "cat x*.ALLchr > ".$fastq.'.ALLchr';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mkdir bits".'.'.$fastq.'.ALLchr';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mv x* bits".'.'.$fastq.'.ALLchr';
 die "Problem stage $n" if $status{$n}; $n++;

 $sizeBits=`wc -l bits.$fastq.ALLchr/*ALLchr`;
 
 if ($sizeBits =~ /total/)
 {
 $sizeBits=~ /(\d{1,15})\stotal/;
 $sizeBits=$1;
 }
 else
 {
 $sizeBits=~ /(\d{1,15})/;
 $sizeBits=$1;
 }

 #print "$sizeBits\n";
 my $sizeRec=`wc -l $fastq.ALLchr`;
 $sizeRec=~ /(\d{1,15})/;
 $sizeRec=$1;
 #print "$sizeRec\n";
 #print $sizeBits - $sizeRec."\n";
 if ($sizeBits != $sizeRec)
 {
  die "ALLchr reconstruction failed\n";
 }
 else
 {
  print "ALLchr reconstruction successfull\n";
 }
 }
 if ($options =~ /[CF]/)
 {
####En travaux####

 if ($options =~ /[U]/)
 {
  print "Reversing coordinates for dUTP protocol...\n";
  $status{$n}=system "cp ".$fastq.'.ALLchr '.$fastq.'.U.ALLchr';
  $status{$n}=system "perl ".$scriptp."dUTPconvert.pl ".$fastq.'.U.ALLchr '.$fastq.'.ALLchr';
  die "Problem stage $n" if $status{$n}; $n++;
 }
##################
##################

 #perl ~/POMBE_SEQ/analysis/SCRIPTS/MAPfilter.pl 0510_1.fasta N0510_1.ALLchr N0510_1.5mis.ALLchr 51 10 &
 print "Filtering...\n";
 $status{$n}=system "perl ".$scriptp."MAPfilter.pl ".$fastq.'.fasta '.$fastq.'.ALLchr '.$fastq.'.5mis.ALLchr '.$read." 10";
 die "Problem stage $n" if $status{$n}; $n++;

 print "collapsing reads...\n";
 $status{$n}=system "perl ".$scriptp."MAPcollapse.pl ".$fastq.'.5mis.ALLchr';
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ~/POMBE_SEQ/analysis/SCRIPTS/diagnostic.pl N0510_1.5mis.ALLchr 51 &
 print "Running diagnostic...\n";
 $status{$n}=system "perl ".$scriptp."diagnostic.pl ".$fastq.'.5mis.ALLchr '.$read;
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ~/POMBE_SEQ/analysis/SCRIPTS/ReadFilter.pl N0510_1.5mis.ALLchr 0510_1.fasta &
 print "Extracting un-mapped reads...\n";
 $status{$n}=system "perl ".$scriptp."ReadFilter.pl ".$fastq.'.5mis.ALLchr '.$fastq.'.fasta';
 die "Problem stage $n" if $status{$n}; $n++;
 }
}
#################ALLtrs mapping and filtering#########################################################
if ($options =~ /T/)
{
 print "Preparing for ALLtrs mapping...\n";
 $status{$n}=system "split -d -l 2000000 ".$fastq.'.5mis.ALLchr.left';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "sh ".$scriptp."mapping_transcripts.090511.sh";
 die "Problem stage $n" if $status{$n}; $n++;

 print "Cleaning up...\n";
 $status{$n}=system "cat x*.ALLtrs > ".$fastq.'.ALLtrs';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mkdir bits".'.'.$fastq.'.ALLtrs';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mv x* bits".'.'.$fastq.'.ALLtrs';
 die "Problem stage $n" if $status{$n}; $n++;

 $sizeBits=`wc -l bits.$fastq.ALLtrs/*ALLtrs`;
 if ($sizeBits =~ /total/)
 {
 $sizeBits=~ /(\d{1,15})\stotal/;
 $sizeBits=$1;
 }
 else
 {
 $sizeBits=~ /(\d{1,15})/;
 $sizeBits=$1;
 }
 #print "$sizeBits\n";
 my $sizeRec=`wc -l $fastq.ALLtrs`;
 $sizeRec=~ /(\d{1,15})/;
 $sizeRec=$1;
 #print "$sizeRec\n";
 #print $sizeBits - $sizeRec."\n";
 if ($sizeBits != $sizeRec)
 {
  die "ALLtrs reconstruction failed\n";
 }
 else
 {
  print "ALLtrs reconstruction successfull\n";
 }

 if ($options =~ /[U]/)
 {
  print "Reversing coordinates for dUTP protocol...\n";
  $status{$n}=system "cp ".$fastq.'.ALLtrs '.$fastq.'.U.ALLtrs';
  $status{$n}=system "perl ".$scriptp."dUTPconvert.pl ".$fastq.'.U.ALLtrs '.$fastq.'.ALLtrs';
  die "Problem stage $n" if $status{$n}; $n++;
 }

 #perl ~/POMBE_SEQ/analysis/SCRIPTS/MAPfilter.pl 0510_1.fasta N0510_1.ALLtrs N0510_1.5mis.ALLtrs 51 10
 print "Filtering...\n";
 $status{$n}=system "perl ".$scriptp."MAPfilter.pl ".$fastq.'.fasta '.$fastq.'.ALLtrs '.$fastq.'.5mis.ALLtrs '.$read." 10";
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ~/POMBE_SEQ/analysis/SCRIPTS/diagnostic.pl N0510_1.5mis.ALLtrs 51
 print "Running diagnostic...\n";
 $status{$n}=system "perl ".$scriptp."diagnostic.pl ".$fastq.'.5mis.ALLtrs '.$read;
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ~/POMBE_SEQ/analysis/SCRIPTS/ReadFilter.pl N0510_1.5mis.ALLtrs N0510_1.5mis.ALLchr.left
 print "Extracting un-mapped reads...\n";
 $status{$n}=system "perl ".$scriptp."ReadFilter.pl ".$fastq.'.5mis.ALLtrs '.$fastq.'.5mis.ALLchr.left';
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ~/POMBE_SEQ/analysis/SCRIPTS/trsTOgen.pl ~/POMBE_SEQ/analysis/gff_090511.txt N0510_1.5mis.ALLtrs
 print "converting trs coordinates...\n";
 $status{$n}=system "perl ".$scriptp."trsTOgen.pl ".$gffp.$gff." ".$fastq.'.5mis.ALLtrs';
 die "Problem stage $n" if $status{$n}; $n++;
 
 print "Collapsing reads...\n";
 $status{$n}=system "perl ".$scriptp."MAPcollapse.pl ".$fastq.'.5mis.ALLtrs.GEN';
 die "Problem stage $n" if $status{$n}; $n++;


 ##perl ~/POMBE_SEQ/analysis/SCRIPTS/mapREADS.pl N0510_1.5mis.ALLtrs.GEN ~/POMBE_SEQ/analysis/gff_090511.txt
 ##$status{$n}=system "perl ".$scriptp."mapREADS.pl ".$fastq.'.5mis.ALLtrs.GEN '.$gffp.$gff;
 ##die "Problem stage $n" if $status{$n}; $n++;
}

##################Per bases pair hits. Needs further formating in R###########################################
if ($options =~ /V/)
{
 #perl ~/POMBE_SEQ/analysis/SCRIPTS/tocount.ALLchr.pl N0510_1.5mis.ALLchr
 print "Creating count files...\n";
 $status{$n}=system "perl ".$scriptp."tocount.ALLchr.pl ".$fastq;
 die "Problem stage $n" if $status{$n}; $n++;
 #perl ~/POMBE_SEQ/analysis/SCRIPTS/tocount.ALLtrs.pl N0510_1.5mis.ALLtrs.GEN ~/POMBE_SEQ/analysis/gff_090511.txt
 $status{$n}=system "perl ".$scriptp."tocount.ALLtrs.pl ".$fastq.' '.$gffp.$gff;
 die "Problem stage $n" if $status{$n}; $n++;
 
 print "Formating counts files and creating .rda files for visualisation...\n";
 ##use unique count from one file##
 #$chr_depth=`grep 'unique' *chr.diag*`;
 #$chr_depth =~ /\d+$/;
 #$chr_depth = $&;
 #$trs_depth=`grep 'unique' *trs.diag*`;
 #$trs_depth =~ /\d+$/;
 #$trs_depth = $&;
 #$depth=($chr_depth + $trs_depth) / 10000000;
 #print "$chr_depth\n$trs_depth\n$depth\n";
 
 ##use total number of reads after collapse for NORM file creation##
if ($options =~ /B/)
 {
 $chr_depth=`wc -l $fastq.5mis.ALLchr.col`;
 $chr_depth=~ /(\d{1,15})/;
 $chr_depth=$1;
 $depth=($chr_depth) / 10000000;
 } 

else
 {
  $depth_test=`ls -l *5mis.ALLchr.col | wc -l`;
  if ($depth_test == 1)
  {
  $chr_depth=`wc -l $fastq.5mis.ALLchr.col`;
  $chr_depth=~ /(\d{1,15})/;
  $chr_depth=$1;
  $trs_depth=`wc -l $fastq.5mis.ALLtrs.GEN.col`;
  $trs_depth=~ /(\d{1,15})/;
  $trs_depth=$1;
  }
  else
  {
  $chr_depth=`wc -l $fastq*.5mis.ALLchr.col`;
  $chr_depth=~ /(\d{1,15})\stotal/;
  $chr_depth=$1;
  $trs_depth=`wc -l $fastq*.5mis.ALLtrs.GEN.col`;
  $trs_depth=~ /(\d{1,15})\stotal/;
  $trs_depth=$1;
  }
  $depth=($chr_depth + $trs_depth) / 10000000;
 }

 open (OUT, ">", "BASES_".$fastq.'.r') or die 'could not open BASES R file';
 print OUT 'source("'.$scriptp.'CreateChromosomes.r")'."\n";
 print OUT 'source("'.$scriptp.'toLINEAR.r")'."\n";
 print OUT 'source("'.$scriptp.'format.bases.r")'."\n";
 if($options =~ /B/)
 {
 print OUT 'CreateChromosomes.CHIP(f1.dir="",f1="'.$fastq.'",depth='.$depth.')'."\n";
 }
 else
 {
 print OUT 'CreateChromosomes(f1.dir="",f1="'.$fastq.'",depth='.$depth.')'."\n";
 }
 print OUT 'q()'."\n";
 close OUT;

 $status{$n}=system 'R --vanilla --slave < BASES_'.$fastq.'.r > BASES_'.$fastq.'.out';

}

#################RPKM analysis#################################################################################
if ($options =~ /R/)
{
 unless ($options =~ /r/)
 {
 #perl ~/POMBE_SEQ/analysis/SCRIPTS/mapREADS.IV.pl N0510_1.5mis.ALLchr gff_090511.txt
 print "Mapping reads to annotation...\n";
 $status{$n}=system "perl ".$scriptp."mapREADS.IV.pl ".$fastq.'.5mis.ALLchr.col '.$gffp.$gff;
 die "Problem stage $n" if $status{$n}; $n++;
 #perl ~/POMBE_SEQ/analysis/SCRIPTS/mapREADS.V.pl N0510_1.5mis.ALLtrs.GEN gff_090511.txt
 $status{$n}=system "perl ".$scriptp."mapREADS.V.pl ".$fastq.'.5mis.ALLtrs.GEN.col '.$gffp.$gff;
 die "Problem stage $n" if $status{$n}; $n++;
 }
 unless ($options =~ /m/)
 {
 #perl SCRIPTS/READScount.pl 0510_1/0510_1.5mis.ALLchr.mapIV.gff_090511.txt  gff_090511.txt &
 print "Creating read counts tables...\n";
 $status{$n}=system "perl ".$scriptp."READScount.pl ".$fastq.' '.$gffp.$gff;
 die "Problem stage $n" if $status{$n}; $n++;
 #perl SCRIPTS/READScount.pl 0510_1/0510_1.5mis.ALLtrs.GEN.mapV.gff_090511.txt  gff_090511.txt
 #$status{$n}=system "perl ".$scriptp."READScount.pl ".$fastq.' '.$gffp.$gff;
 #die "Problem stage $n" if $status{$n}; $n++;

 print "Creating RPKM tables...\n";
 #$chr_depth=`grep 'unique' *chr.diag*`;
 #$chr_depth =~ /\d+$/;
 #$chr_depth = $&;
 #$trs_depth=`grep 'unique' *trs.diag*`;
 #$trs_depth =~ /\d+$/;
 #$trs_depth = $&;
 ##use total depth from collapsed files for RPKM calculation
 
 ##use total number of reads after collapse for NORM file creation##
 $depth_test=`ls -l *5mis.ALLchr.col | wc -l`;
 if ($depth_test == 1)
 {
 $chr_depth=`wc -l $fastq.5mis.ALLchr.col`;
 $chr_depth=~ /(\d{1,15})/;
 $chr_depth=$1;
 $trs_depth=`wc -l $fastq.5mis.ALLtrs.GEN.col`;
 $trs_depth=~ /(\d{1,15})/;
 $trs_depth=$1;
 }
 else
 {
 $chr_depth=`wc -l $fastq*.5mis.ALLchr.col`;
 $chr_depth=~ /(\d{1,15})\stotal/;
 $chr_depth=$1;
 $trs_depth=`wc -l $fastq*.5mis.ALLtrs.GEN.col`;
 $trs_depth=~ /(\d{1,15})\stotal/;
 $trs_depth=$1;
 }
 $depth=($chr_depth + $trs_depth) / 1000000;
 
 #print "$chr_depth\n$trs_depth\n$depth\n";

 open (OUT, ">", "RPKM_".$fastq.'.r') or die 'could not open RPKM R file';
 print OUT 'source("'.$scriptp.'computeRPKM.r")'."\n";
 print OUT 'gff=read.delim("'.$gffp.$gff.'",stringsAsFactors=FALSE)'."\n";
 #depth1=26.282949
 print OUT 'R'.$fastq.'=computeRPKM("'.$fastq.'.5mis",gff.name="'.$gff.'",depth='.$depth.',gff=gff)'."\n";
 print OUT 'save(list=c("R'.$fastq.'"),file="RPKM_'.$fastq.'_'.$tag.'.rda")'."\n";
 print OUT 'q()'."\n";
 close OUT;

 $status{$n}=system 'R --vanilla --slave < RPKM_'.$fastq.'.r > RPKM_'.$fastq.'.out';
 }
}

#######SPLICING####################################################################################################
if ($options =~ /S/)
{
 print "Performing Splicing analysis...\n";
 #perl ../SCRIPTS/junctionMAP.pl 0510_7.5mis ../gff_260310.txt 0
 $status{$n}=system 'perl '.$scriptp.'junctionMAP.pl '.$fastq.'.5mis '.$gffp.$gff.' 0';
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ../SCRIPTS/getREADS.pl 0510_7.5mis.junctions.intronreads 0510_7.fasta
 $status{$n}=system 'perl '.$scriptp.'getREADS.pl '.$fastq.'.5mis.junctions.intronreads '.$fastq.'.fasta';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system 'cp '.$fastq.'.5mis.junctions.intronreads.selected x00.fasta';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system 'sh '.$scriptp.'mapping_intron.090511.sh';
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ../SCRIPTS/MAPfilter.pl x00.fasta x00.ALLtrs x00.5mis.ALLtrs 50 10 &
 
 if ($options =~ /[U]/)
 {
  print "Reversing coordinates for dUTP protocol...\n";
  $status{$n}=system 'cp x00.ALLtrs x00.U.ALLtrs';
  $status{$n}=system "perl ".$scriptp.'dUTPconvert.pl x00.U.ALLtrs x00.ALLtrs';
  die "Problem stage $n" if $status{$n}; $n++;
 }

 print "filtering re-mapped intron reads...\n";
 $status{$n}=system 'perl '.$scriptp.'MAPfilter.pl x00.fasta x00.ALLtrs x00.5mis.ALLtrs '.$read.' 10';
 die "Problem stage $n" if $status{$n}; $n++;

 #perl ../SCRIPTS/junctionMAP.pl 0510_7.5mis ../gff_260310.txt x00.5mis.ALLtrs
 $status{$n}=system 'perl '.$scriptp.'junctionMAP.pl '.$fastq.'.5mis '.$gffp.$gff.' x00.5mis.ALLtrs';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mkdir bits".'.'.$fastq.'.Splicing';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mv x* bits".'.'.$fastq.'.Splicing';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mv *junction* bits".'.'.$fastq.'.Splicing';
 die "Problem stage $n" if $status{$n}; $n++;

 $status{$n}=system "mv bits".'.'.$fastq.'.Splicing/*junctions.filtered .';
 die "Problem stage $n" if $status{$n}; $n++;

 print "Creating splicing table...\n";
 open (OUT, ">", "SPLICE_".$fastq.'.r') or die 'could not open SPLICE R file';
 print OUT 'source("'.$scriptp.'CreateGeneSplicingTable.r")'."\n";
 print OUT 'trans=read.delim("'.$fastq.'.5mis.junctions.filtered",stringsAsFactors=FALSE)'."\n";
 print OUT 'load("RPKM_'.$fastq.'_'.$tag.'.rda")'."\n";
 print OUT 'S'.$fastq.'=CreateGeneSplicingTable(tr=trans,ei=R'.$fastq.')'."\n";
 print OUT 'save(list=c("S'.$fastq.'"),file="SPLICE_'.$fastq.'_'.$tag.'.rda")'."\n";
 print OUT 'q()'."\n";
 close OUT;

 $status{$n}=system 'R --vanilla --slave < SPLICE_'.$fastq.'.r > SPLICE_'.$fastq.'.out';
 die "Problem stage $n" if $status{$n}; $n++;

}
if ($options =~ /X/)
{
 #my @sub=(20000000,15000000,10000000,5000000,2500000,1000000,500000,100000);
 my @sub=(500000,100000);
 
 $status{$n}=system "mkdir subsampling".'.'.$fastq;
 die "Problem stage $n" if $status{$n}; $n++;
 
 foreach(@sub)
 {
#  print "###subsampling $_ reads\n";
#  $status{$n}=system 'perl '.$scriptp.'RandomLineSelect.pl '.$fastq.'.fasta '.$fastq.'.5mis.ALLchr '.$_.' MAP';
#  die "Problem stage $n" if $status{$n}; $n++;
#
#  print "Running diagnostic...\n";
#  $status{$n}=system "perl ".$scriptp."diagnostic.pl ".$fastq.'.5mis.ALLchr.RANDOM.'.$_.'.MAP '.$read;
#  die "Problem stage $n" if $status{$n}; $n++;
#   
#  print "Mapping reads to annotation...\n";
#  $status{$n}=system "perl ".$scriptp."mapREADS.IV.pl ".$fastq.'.5mis.ALLchr.RANDOM.'.$_.'.MAP '.$gffp.$gff;
#  die "Problem stage $n" if $status{$n}; $n++;
# 
#  print "Creating read counts tables...\n";
#  $status{$n}=system "perl ".$scriptp."READScount.pl ".$fastq.'.5mis.ALLchr.RANDOM.'.$_.'.MAP.mapIV.'.$gff." ".$gffp.$gff;
#  die "Problem stage $n" if $status{$n}; $n++;
# 
#  print "Creating RPKM tables...\n";
#  $chr_depth=$fastq.'.5mis.ALLchr.RANDOM.'.$_.'.MAP.diagnostic';
#  $chr_depth=`grep 'unique' $chr_depth`;
#  $chr_depth =~ /\d+$/;
#  $chr_depth = $&;
#  $depth=($chr_depth) / 1000000;
#
#  open (OUT, ">", "RPKM_".$fastq.'.'.$_.'.r') or die 'could not open RPKM R file';
#  print OUT 'source("'.$scriptp.'computeRPKM.r")'."\n";
#  print OUT 'gff=read.delim("'.$gffp.$gff.'",stringsAsFactors=FALSE)'."\n";
#  #depth1=26.282949
#  print OUT 'R'.$fastq.'.'.$_.'=computeRPKM1("'.$fastq.'.5mis.ALLchr.RANDOM.'.$_.'.MAP",gff.name="'.$gff.'",depth='.$depth.',gff=gff)'."\n";
#  print OUT 'save(list=c("R'.$fastq.'.'.$_.'"),file="RPKM_'.$fastq.'.'.$_.'_'.$tag.'.rda")'."\n";
#  print OUT 'q()'."\n";
#  close OUT;
#
#  $status{$n}=system 'R --vanilla --slave < RPKM_'.$fastq.'.'.$_.'.r > RPKM_'.$fastq.'.'.$_.'.out';
#  die "Problem stage $n" if $status{$n}; $n++;
  
  $status{$n}=system 'mv *'.$_.'* subsampling.'.$fastq;
  die "Problem stage $n" if $status{$n}; $n++;
 }
}
 
###clean up###
$status{$n}=system "rm -f *bases";
die "Problem stage $n" if $status{$n}; $n++;
$status{$n}=system "rm -rf bits*";
die "Problem stage $n" if $status{$n}; $n++;


#calculate time taken
my $time1 = new Benchmark;
my $timdiff = timediff($time1, $time0);
print "Running pipeline took ", timestr($timdiff), "\n";
#print "status outputs ", Dumper \%status, "\n\n";


