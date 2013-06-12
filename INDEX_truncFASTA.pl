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
  if($line=~/^([ATCGN]{6})([ATCGN]{6})TTTTTTTTTT/)
 {
  $header="$header;".$1.';'.$2;
  $line=substr($line,$start,$length);
 }
 else
 {
  $header="$header;NA;NA";
  $line=substr($line,$start,$length);
 }
 print "$header\n$line\n";
 $count=0;
 }
}
close IN;
