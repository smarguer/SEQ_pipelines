#!/usr/bin/perl


use strict;
use warnings;

my $line;
my $count=0;
my $count1=0;
my %counter;
my $start=0;
#my $length=40;
#my $start=0;
#my $length=58;
#my $start=35;

#my $start=0;
my $length=50;
#

if (@ARGV != 2) {die "need 2 file names";}
(my $in, my $out)=@ARGV;

open (IN, $in) or die 'could not find the input file';
open (OUT, ">", $out) or die 'could not open output file';

while ($line=<IN>)
{
	$count++;
	if ($count == 1)
	{
	print OUT "$line";
	#if (!$counter{$line}){$counter{$line}=0;}
	#$counter{$line}++;
	}
	elsif ($count == 2)
	{
	$line=substr($line,$start,$length);
	print OUT "$line\n";
	$count=0;
	}
	
}


#for my $out (sort keys %counter)
#{
#if ($counter{$out} > 1){$count1++;}
#}
#print "$count1\n";


close IN;
close OUT;


