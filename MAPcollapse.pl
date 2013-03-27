#!/usr/bin/perl


use strict;
use warnings;
use POSIX;
my $line;
my $line1;
my $count=0;
my $count1=0;
my $holder;
my @holder;
my @holder1;
my $pos;
my $mpos;
my $score=0;
my $chrom='start';
my $sense;
my $readlength;
my %readlength;
my $mis=0;
my $len;
my $steps;
my $out1;
my $up="start";
my @hold=();
my $index;
my $M;
my $input=0;
my %recordIndex;
my %inputIndex;
#######################
#Takes input
#######################

if (@ARGV != 1) {die "Need one mapping file to collapse and an output name";}
(my $in)=@ARGV;

open (IN, $in) or die 'could not find the mapping output';
open (OUT, ">", $in.'.col') or die 'could not open output file';

if (-e $in.'.recordIndex')
{
 $input=1;
 open(IN1,$in.'.recordIndex') or die 'could not open  input index file';
 while($line1=<IN1>)
 {
  @holder1=split('\t',$line1);
  $inputIndex{$holder1[0]}=$holder1[1];
 }
}
else
{
 open (OUT1, ">", $in.'.recordIndex') or die 'could not open recordIndex file';
}
#open (OUT1,">",$out.".log") or die 'could not open log file';
######################################
#parses fasta file
######################################
#while ($line1=<IN1>)
#{
##	chomp ($line1);
##	if ($line1 =~ /[\#]/)
##        {
##         next;
##        }
#    if ($line1 =~ /[\>]/)
#	{
#	 $out1=substr($line1,1);
#	}
#	else
#	{
##	 $readlength{$out1}=length $line1;
#	 $readlength{$out1}=$readlength;
#	}
#}
#
######################################
#parses Exonerate output
######################################


while ($line=<IN>)
{
	
	$count++;
	#print "$count\n";
        $count1++;
	chomp ($line);

	@holder = split (/ /, $line);
        unless ($holder[0]){next;}
	if ($holder[0] ne "vulgar:"){next;}
	
	#$len=$holder[3] - $holder[2];
################################################################################
	if ($holder[1] eq $up)
        {
         $up=$holder[1];
         push (@hold, $line);
         $M=1;
         next;
	}
        elsif (@hold ==1)
        {
         $up=$holder[1];
         print OUT "U$hold[0]\n";
         @hold=();
         push (@hold, $line);
         $M=0;
         next;         
        } 
        elsif (@hold >1)
        {
         #print "ok\n";
         $up=$holder[1];
	 if($input==1)
         {
          $index=$inputIndex{$count};
         }
         else
         {
          $index = int(rand(@hold));         
          $recordIndex{$count}=$index;
         }
         #print "$index\n";
         print OUT "M$hold[$index]\n";
         @hold=();
         push (@hold, $line);
         $M=0;
         next;         
        }
        elsif($up eq "start")
        {
         $up = $holder[1];
         push (@hold, $line);
        }
 
}

if ($M==1)
{
 $index = int(rand(@hold));         
 print OUT "M$hold[$index]\n";
}

if($input==0)
{
 for my $out (keys %recordIndex)
 {
 print OUT1 "$out\t$recordIndex{$out}\n";
 }
 close OUT1;
}
close IN;
close OUT;
