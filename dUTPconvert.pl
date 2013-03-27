#!/usr/bin/perl


use strict;
use warnings;
use POSIX;
my $line;
my @holder;
my $start;
my $end;

#######################
#Takes input
#######################

if (@ARGV != 2) {die "I need 2 file names";}
(my $in, my $out)=@ARGV;

open (IN, $in) or die 'could not find the mapping output';
open (OUT, ">", $out) or die 'could not open output file';

######################################
#parses and reverses Exonerate output
######################################


while ($line=<IN>)
{
	chomp ($line);
	@holder = split (/ /, $line);
	unless ($holder[0]){next;}
	if ($holder[0] ne "vulgar:"){next;}
	$start=$holder[6];
    $end=$holder[7];
    $holder[6]=$end;
    $holder[7]=$start;
    if ($holder[8] eq '-')
    {
     $holder[8]='+';
    }
    elsif ($holder[8] eq '+')
    {
     $holder[8]='-';
    }
    
    print OUT join(' ', @holder)."\n";
}

close IN;
close OUT;
