#!/usr/bin/perl

use strict;
use warnings;

my $line;
my $count=0;
my $count1=0; 
my %base;
my @holder;
my @holder1;
my @toDo;
my @toDo_in;
my @chr;
my $i;
my $j;
my $pos;
my $mpos;
my $chr;
my $lchr;
my $nchr;
my %lchr;

###all in one
###processes all .ALLchr files with appropriate names in the folder

if (@ARGV != 1) {die "need 1 name";}
(my $in)=@ARGV;

@toDo=`ls $in*.5mis.ALLchr.col`;
#print "@toDo\n";

for (@toDo){$_=~s/^\s+|\s+$//g};

#print "To process: \n";
#print join("\n",@toDo)."\n";

#die;


open (OUT1, ">", "$in.5mis.ALLchr.col.chr1.bases") or die 'could not open output file';
open (OUT2, ">", "$in.5mis.ALLchr.col.chr2.bases") or die 'could not open output file';
open (OUT3, ">", "$in.5mis.ALLchr.col.chr3.bases") or die 'could not open output file';

##Create data structure
@chr=("CH1_bases","CH2_bases","CH3_bases");

$lchr{$chr[0]}=5579133;
$lchr{$chr[1]}=4539804;
$lchr{$chr[2]}=2452883;

foreach(@chr)
{
 $nchr=$_;
 for ($i =1; $i <= (2 * $lchr{$nchr}); $i++)
 {
  $base{$nchr}{$i}=0;
  #print "$i";
 }
}

foreach (@toDo)
{
 open (IN, $_) or die "could not find the input file named $_";
 print "processing $_\n";

 while ($line=<IN>)
 {
  chomp ($line);
  @holder = split (/ /, $line);
  if ($holder[0] !~ /vulgar/){next;}
  unless (grep /$holder[5]/, @chr){next;}
  $count++;
  @holder1= ($holder[6],$holder[7]);
  @holder1= sort (@holder1);
	
  if ($holder[8] eq "-")
  {
   $holder1[0]=$holder1[0]+$lchr{$holder[5]};
   $holder1[1]=$holder1[1]+$lchr{$holder[5]};
  } 		
#to adjust the start of the match coordinate to the real genome coordinate.the
 
  $holder1[0]=$holder1[0]+1;
#if ($holder1[0]==0){$holder1[0]=1;}
#if ($holder1[1]==0){$holder1[1]=1;}

  foreach $pos ($holder1[0] .. $holder1[1])
  {
   $base{$holder[5]}{$pos}++;
  }
 }
 close IN;
}	
#print "printing...\n";

##print output##
for my $out1 (sort keys %{$base{$chr[0]}})
{
 print OUT1 "$out1\t$base{$chr[0]}{$out1}\n";
}
for my $out2 (sort keys %{$base{$chr[1]}})
{
 print OUT2 "$out2\t$base{$chr[1]}{$out2}\n";
}
for my $out3 (sort keys %{$base{$chr[2]}})
{
 print OUT3 "$out3\t$base{$chr[2]}{$out3}\n";
}

close OUT1;
close OUT2;
close OUT3;





 
