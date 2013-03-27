#!/usr/bin/perl


use strict;
use warnings;
use List::Util qw[min max];

if (@ARGV != 4) {die "need two files and a number of reads and an read option [FASTA/MAP]";}
(my $in, my $in1, my $number, my $map)=@ARGV;

unless (($map eq "FASTA") || ($map eq "MAP")) {die "Unknown read option [FASTA/MAP]"};

open (IN, $in) or die 'could not find the fasta file';
open (IN1, $in1) or die 'could not find the mapping output file';

my $line;
my $line1;
my $randomNumber;
my %read;
my %map;
my %seq;
my $out;
my $name;
my $count=0;
my $count1=0;
my $count2=0;
my $count3=0;
my @holder;

open (OUT, ">", "./$in1.RANDOM.$number.$map") or die 'could not open output file';

while ($line=<IN>)
{
 $count++;
 unless ($line=~/\>/){next;}
 $randomNumber=int(rand(1000000000));
#print"$randomNumber\n"; 
 chomp $line;
 $read{$randomNumber}=$line;
 #$read{$randomNumber}=$line;
#print "$read{$randomNumber}\n"
}
#print "fasta read\n";

while ($line1=<IN1>)
{
 $count1++;
 chomp $line1;
 @holder = split (/ /, $line1);
 unless ($holder[0]){next;}
 if ($holder[0] ne "vulgar:"){next;}
 $name=">$holder[1]";
#print "$name\n";

 if($map{$name})
 {
  $map{$name} .= "\n$line1";
 }
 else
 {
  $map{$name}=$line1; 
 }
}

#print "map read\n";

for $out (sort {$a<=>$b} keys %read)
{
 $count2++;
 if ($map eq "FASTA")
 {
  last if ($count2 > $number);
 }
 if ($map eq "MAP")
 {
  last if ($count3 >= $number);
 }
 if ($map{$read{$out}})
 {
#  print $map{$read{$out}}."\n";
  $count3++;
  print OUT $map{$read{$out}}."\n";
 }
 else
 {
  next;
 }
}


close IN;
close IN1;
close OUT;
