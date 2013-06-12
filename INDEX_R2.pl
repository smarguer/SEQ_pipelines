#!/usr/bin/perl

use strict;
use warnings;

my $line;
my $line1;

my $count=0;
my @holder;
my @holder1;
my @holder2;
my %index;

if (@ARGV != 1) {die "wrong number of files";}
(my $in)=@ARGV;

open (IN, $in.'_R1_50.fasta') or die 'could not find the input file';
open (IN1, $in.'_R2.fasta') or die 'could not find the input file';

while ($line=<IN>)
 {
  $count++;
  chomp ($line);
  if($line=~/^>/)
  {
  @holder = split (/\t/, $line);
  @holder2 = split (/;/,$holder[0]);
  $index{$holder2[0]}="$holder2[3];$holder2[4]";  
  }
 }
 while ($line1=<IN1>)
 {
  chomp $line1;
  @holder1=split(/;/,$line1);

  if($line1=~/^>/)
  {
   print "$line1;$index{$holder1[0]}\n";
  }
  else
  {
   print substr ($line1,11,50)."\n";
  }
 }

close IN;
close IN1;

