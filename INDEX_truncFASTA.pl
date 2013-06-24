#!/usr/bin/perl


use strict;
use warnings;

my $line;
my $count=0;
my $count1=0;
my %counter;
my $start=35;
my $length=50;
my $header;
my $lineOut;
my $RNT=8;
#

if (@ARGV != 1) {die "need 2 file names";}
(my $in, my $out)=@ARGV;

open (IN, $in.'.fasta') or die 'could not find the input file';

while ($line=<IN>)
{
 $count++;
 if ($count == 1)
 {
  $header=$line;
  chomp $header;
 }
 elsif ($count == 2)
 {
  if($line=~/^([ATCGN]{$RNT})([ATCGN]{6})TTTTTTTTTT/)
 {
  $header="$header;".$1.';'.$2;
##
  $lineOut = reverse($');
  $lineOut =~ tr/ACGTacgt/TGCAtgca/;
  $lineOut =~ /[A]+$/;
  $lineOut =$`;
####TO ALLOW ONE MIS-MATCH###
#chop($lineOut);
#$lineOut =~ /[A]+$/;
#############################
  $lineOut = substr($lineOut,-45);
  $lineOut = $lineOut."AAAAA";
##
##
#  $line=substr($line,$start,$length);
##
 }
 else
 {
  $header="$header;NA;NA";
##
  $lineOut = substr ($line,15);
  $lineOut = reverse($lineOut);
  $lineOut =~ tr/ACGTacgt/TGCAtgca/;
  if($lineOut =~ /[A]+$/)
  {
   $lineOut =$`;
####TO ALLOW ONE MIS-MATCH###
#chop($lineOut);
#$lineOut =~ /[A]+$/;
#############################
   $lineOut = substr($lineOut,-45);
   $lineOut = $lineOut."AAAAA";
  }
  else
  {
   $lineOut = substr($lineOut,-50);
  }
##
 }
 print "$header\n$lineOut\n";
 $count=0;
 }
}
close IN;
