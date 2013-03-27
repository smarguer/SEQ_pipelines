#!/usr/bin/perl
################################
##works with .ALLtrs and .ALLchr files as sources of exon-exon and exon-intron
##reads respectively
################################

use strict;
use warnings;
use List::Util qw[min max];


if (@ARGV != 3) {die "need an experiment name and a gff file, and an option";}
(my $in, my $gff, my $option)=@ARGV;

my $in1=$gff;

open (IN, $gff) or die 'could not find the gff file';
open (IN1, "$in.ALLtrs") or die 'could not find the transcripts input';
open (IN2, "$in.ALLchr") or die 'could not find the cds input';
open (IN3, $in1) or die 'could not find the gff file';

if (-f $option)
{
 open (IN4, $option) or die 'could not find the unspecific reads file';
 open (OUT, ">", "$in.junctions.filtered") or die 'could not open output file 4';
 open (OUT1, ">", "$in.junctions.not.filtered") or die 'could not open output file 2';
 open (OUT2, ">", "$in.junctions.intronreads.filtered") or die 'could not open output file 3';
}

else
{
open (OUT, ">", "$in.junctions") or die 'could not open output file';
open (OUT1, ">", "$in.junctions.not") or die 'could not open output file 2';
open (OUT2, ">", "$in.junctions.intronreads") or die 'could not open output file 3';
}

my $line=<IN>;
my $line1;
my $line2;
my $line3=<IN3>;

my $count=0;
my $count1=0;
my $count3=0;
my $count4=0;
my $not;
my %CDS;
my %INT;
my %startG;
my %endG;
my %chrG;
my %strandG;
my %count;
my %countint;
my %readmatch;
my %intronFilter;
my @holder;
my @holder1;
my @not;
my @formatedALLchr;
my $i=0;
my $k;
my $chr;
my $chrom;
my $end1;
my $mpos;
my $lchr;
my $nchr;
my $total;
my $totalint;
my $featureNumber=0;
my $featureNumber1=0;
my @features=('CDS','misc_RNA','tRNA','snRNA','snoRNA','LTR','rRNA','mRNA');
my @excluded=('utr_intron','CDS_motif','BLASTN_HIT','gap','repeat_unit','mRNA','rep_origin','polyA_signal','LTR','tRNA','5\'UTR','3\'UTR');
my %INTcount;
my @holder2;
my $maxINT;
my $start;
my $end;

my $test='SPAC8C9.05';   
my @holder4;
my %duplicate;



###############################################################
#produce junctions data structure from gff
###############################################################
while ($line3=<IN3>)
{
 chomp ($line3);
 @holder2 = split (/\t/, $line3);

 if ($holder2[8] eq 'No comment'){next;}

 elsif ($holder2[2] eq "intron")
 {
 $INTcount{$holder2[9]}++;
 }

}
close IN3;

#print "preparing junctions data structure.\n";
while ($line=<IN>)
{
 chomp ($line);
 @holder= split (/\t/,$line);
    
 if (grep /$holder[9]/, @excluded)
 {
  next;
 }

 if ($holder[8] eq 'No comment'){next;}
	
 if ($holder[2] eq 'gene')
 { 
  $startG{$holder[9]}=$holder[3];
  $endG{$holder[9]}=$holder[4];
  $strandG{$holder[9]}=$holder[6];
  if ($strandG{$holder[9]} eq '.')
  {
   $strandG{$holder[9]} = '+';
   print "$holder[6]\n";
  }
  $chrG{$holder[9]}=$holder[12];
 }
 
 elsif (grep /$holder[2]/,@features)
 {
  unless ($CDS{$holder[9]})
  {
   $count=-1;
  }
 $count++;
 $CDS{$holder[9]}[$count]=$holder[4]-$holder[3]+1;	
 }

 elsif ($holder[2] eq 'intron')
 {
  unless ($INT{$holder[9]})
  {
   $count1=-1;
  }
  
  $count1++;
  $INT{$holder[9]}[$count1]=$holder[4]-$holder[3]-1;
  
  if ($holder[6] eq '-')
  {
   $featureNumber1=abs($INTcount{$holder[9]}-($count1+1))+1;
  }
  else
  {
   $featureNumber1=$count1+1;
  }
  
  $countint{$holder[9]}{$featureNumber1}{START}=$holder[3]+1;
  $countint{$holder[9]}{$featureNumber1}{END}=$holder[4]-1;

#print "$featureNumber1\t$count1\n";

#if ($holder[9] eq 'SPAC8C9.05')
#if ($holder[9] eq 'SPAC13G6.05c')
#{
#print "$holder[9]\t$featureNumber1\t$countint{$holder[9]}{$featureNumber1}{START}\t$countint{$holder[9]}{$featureNumber1}{END}\n";
#}
  
  $countint{$holder[9]}{$featureNumber1}{COUNTp}=0;
  $countint{$holder[9]}{$featureNumber1}{COUNTm}=0;
  
  for ($i=($holder[3]+1);$i<$holder[4];$i++)
  {
   $intronFilter{$holder[12]}{$i}=$holder[9];	
   $count4++;
  }	
 }

 $i=0;
}



for my $out (sort keys %CDS)
{
 $totalint=0;

#leave out intron-less genes
 if ($#{$CDS{$out}} == 0)
 {
 #print "####################$out#############\n";
 next;
 }


 if ($strandG{$out} eq '-')
 {
  @{$CDS{$out}} = reverse(@{$CDS{$out}});
	
  if ($INT{$out})
  {
   @{$INT{$out}} = reverse(@{$INT{$out}});
   #$maxINT=max(keys %{$countint{$out}});
   #for ($l=1; $l<=$maxINT; $l++)
   #{
   # 
   #}
  }
 }

 $end1=$#{$CDS{$out}}-1;

 for my $i (0..$end1)
 {
  $total += $_ foreach @{$CDS{$out}}[0..$i] ;

###
#unless (defined $CDS{$out}[$i])
#{
#print "$out\n";  
#}
###


  my $j=$i+1;
  
  #if ($INT{$out}[$i])
  #{
  # $totalint += $INT{$out}[$i];
  #}

  $count{$out}{$j}{POS}=$total;
  $count{$out}{$j}{COUNTp}=0;
  $count{$out}{$j}{COUNTm}=0;

#if ($out eq "SPAC8C9.05")
#if ($out eq 'SPAC13G6.05c')
#{
#print "$count{$out}{$j}{POS}\n";
#} 
 
#if ($INT{$out}[$i])
#  {
#   $countint{$out}{$j}{START}=$startG{$out}+$total+$totalint-$INT{$out}[$i]+1;
#   $countint{$out}{$j}{END}=$startG{$out}+$total+$totalint;
#   $countint{$out}{$j}{COUNTp}=0;
#   $countint{$out}{$j}{COUNTm}=0;
#  }

  $total=0;

 } 
}



close IN;


#die;
##############################################################
#creates hash if ambiguous file provied
##############################################################

if (-f $option)
{
 my $line4;
 while ($line4=<IN4>)
 {
 chomp ($line4);
 @holder4 = split (/ /, $line4);
 unless ($holder4[0]){next;}
 if ($holder4[0] ne "vulgar:"){next;}
 $duplicate{$holder4[1]}=1;
 }
}







##########################################################
#Analyse .ALLtrs exonerate output file. Counts transreads
##########################################################

#print "Analysing transreads.\n";

while ($line1=<IN1>)
{
 chomp ($line1);
 @holder = split (/ /, $line1);
 unless ($holder[0]){next;}
 if ($holder[0] ne "vulgar:"){next;}
 unless (defined %{$count{$holder[5]}})
 {
  next;
 }
	
 if ($holder[8] eq '+')
 {
  $not=0;
  $featureNumber=max(keys %{$count{$holder[5]}});
  for ($i=1; $i<=$featureNumber; $i++)
  {
###
#if ($holder[5] eq "SPAC8C9.05")
#{        
#print "$featureNumber\n";
#}
###
  if (($holder[6] < $count{$holder[5]}{$i}{POS}) && ($holder[7] > $count{$holder[5]}{$i}{POS}))
  {
   $count{$holder[5]}{$i}{COUNTp}++;
   $not++;
  }

###
#if ($holder[5] eq "SPAC8C9.05")
#{
#print "$count{$holder[5]}{$i}{COUNTp}\n";
#}
###
  }
 }	

 elsif ($holder[8] eq '-')
 {
  $featureNumber=max(keys %{$count{$holder[5]}});
  $not=0;
  for ($i=1; $i<=$featureNumber; $i++)
  {
   if (($holder[6] > $count{$holder[5]}{$i}{POS}) && ($holder[7] < $count{$holder[5]}{$i}{POS}))
   {
    $count{$holder[5]}{$i}{COUNTm}++;
    $not++;
   }
  }
 }

 if ($not==0)
 {
  push(@not, $holder[1]);
 }

##
#debug
#if ($holder[5] eq "SPAC8C9.05")
#{
#print "$line1\n";
#}

}
#####################################################################
#Analyse exonerate .ALLchr file counts exon-intron junction reads
#####################################################################

#print "Analysing exon-intron junctions\n";

$i=0;

while ($line2=<IN2>)
{
 $count3++;
 chomp ($line2);
 @holder= split (/ /,$line2);
 unless ($holder[0]){next;}
 if ($holder[0] ne "vulgar:"){next;}
 $chrom=substr($holder[5],2,1);
 
##translate hit names       
 if ($intronFilter{$chrom}{$holder[6]})
 {
  $holder[5]=$intronFilter{$chrom}{$holder[6]};
 }
 elsif ($intronFilter{$chrom}{$holder[7]})
 {
  $holder[5]=$intronFilter{$chrom}{$holder[7]};
 }
 else
 {
  next;
 }
 

#if ($holder[5] eq 'SPAC8C9.05')
#{
#print "$holder[5]\n$line2\n";
#}

#die;        
############
#records junction
###########
 unless ($countint{$holder[5]})
 {
  print "$holder[5]\n";
  next;
 }
 
 if ($duplicate{$holder[1]})
 {
  #if ($holder[5] eq $test)
  #{
  #print "REMOVED_$line2\n";
  #}  
  next;
 }

 unless ($countint{$holder[5]}){print "$holder[5]\n";}
 
 if ($holder[6] < $holder[7])
 {
  $start=$holder[6]+1;
  $end=$holder[7];
 }
 elsif ($holder[6] > $holder[7])
 {
  $start=$holder[7]+1;
  $end=$holder[6];
 }


#if ($holder[5] eq 'SPAC8C9.05')
#{
# print "$start\t$end\n";
#}

 #if ($holder[8] eq $strandG{$holder[5]})
 #{
 $featureNumber=max(keys %{$count{$holder[5]}});
 $not=0;
 for ($i=1; $i<=$featureNumber; $i++)
 {
  if ((($start < $countint{$holder[5]}{$i}{START}) && ($end > $countint{$holder[5]}{$i}{START})) ||
		  (($start < $countint{$holder[5]}{$i}{END}) && ($end > $countint{$holder[5]}{$i}{END})))
  {

   print OUT2 "$holder[1]\n";

   if ($holder[8] eq $strandG{$holder[5]})
   {
    $countint{$holder[5]}{$i}{COUNTp}++;
    $not++;
#if ($holder[5] eq $test)
#{
# print "E$i\t$countint{$holder[5]}{$i}{START}\t$countint{$holder[5]}{$i}{END}\n";
# print "$line2\n";
#}
   }
   elsif ($holder[8] ne $strandG{$holder[5]})
   {
    $countint{$holder[5]}{$i}{COUNTm}++;
    $not++;

#if ($holder[5] eq $test)
#{
# print "E$i\t$countint{$holder[5]}{$i}{START}\t$countint{$holder[5]}{$i}{END}\n";
# print "AS_$line2\n";
#}
   } 
  }
 }
}





#die;



####################################################################
#Wrinting output
####################################################################

#print "Writing output file.\n";



print OUT1 join("\n",@not)."\n"."$#not"."\n";

for my $out (sort keys %count)
{
	for my $i (sort keys %{$count{$out}})
	{

##debug
unless (defined %{$count{$out}{$i}})
{
print "CDS\t$out\t$i\n";
next;
}
unless (defined %{$countint{$out}{$i}})
{
print "INT\t$out\t$i\n";
next;
}
##
	
print OUT "$out.".$i."\t$out\t$i\t$count{$out}{$i}{POS}\t$countint{$out}{$i}{START}\t$countint{$out}{$i}{END}\t$count{$out}{$i}{COUNTp}\t".($countint{$out}{$i}{COUNTp}/2)."\t$count{$out}{$i}{COUNTm}\t".($countint{$out}{$i}{COUNTm}/2)."\n";
	} 
}
close IN1;
close IN2;
close OUT;
close OUT1;
close OUT2;
if (-f $option)
{
close IN4;
}

