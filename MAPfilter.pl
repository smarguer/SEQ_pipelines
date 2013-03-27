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

#######################
#Takes input
#######################

if (@ARGV != 5) {die "I need 3 file names, the readlength and the number of mismatches allowed";}
(my $in1, my $in, my $out, $readlength, my $perc)=@ARGV;

open (IN, $in) or die 'could not find the mapping output';
open (IN1, $in1) or die 'could not find fasta file';
open (OUT, ">", $out) or die 'could not open output file';
#open (OUT1,">",$out.".log") or die 'could not open log file';
######################################
#parses fasta file
######################################
#while ($line1=<IN1>)
#{
#	chomp ($line1);
#	if ($line1 =~ /[\#]/)
#        {
#         next;
#        }
#    if ($line1 =~ /[\>]/)
#	{
#	 $out1=substr($line1,1);
#	}
#	else
#	{
#	 $readlength{$out1}=length $line1;
#	 $readlength{$out1}=$readlength;
#	}
#}

######################################
#parses Exonerate output
######################################


while ($line=<IN>)
{
	
	$count++;
	$count1++;
	$score=0;
	chomp ($line);

	@holder = split (/ /, $line);
	unless ($holder[0]){next;}
	if ($holder[0] ne "vulgar:"){next;}
	$len=$holder[3] - $holder[2];
################################################################################
	#$mis=floor($readlength{$holder[1]}/$perc);
	$mis=floor($readlength/$perc);
	if ($count1==100)
        {
         print 'Maximum number of mismatches allowed is: ';
         print "$mis";
         print ' with read length = ';
         print "$readlength\n";
         #print "$readlength{$holder[1]}\n";
	    }
        #$steps=$mis-($readlength{$holder[1]}-$holder[3]+$holder[2]);
        $steps=$mis-($readlength-$holder[3]+$holder[2]);
	if ($steps < 0){$steps=0; }		
#################################################################################		
		
	#if (($holder[9] >= (($readlength{$holder[1]}*5)-($mis*9))) && 
	if (($holder[9] >= (($readlength*5)-($mis*9))) && 
	    ($holder[9] <= ($len * 5)) && 
        ($holder[9] >= (($len * 5) - ($steps * 9))) && 
        #(($readlength{$holder[1]}-$len) <= $mis) &&
        (($readlength-$len) <= $mis) &&
	    ($holder[2] <= 2) &&	
        #($holder[3] >= ($readlength{$holder[1]}-2))
        ($holder[3] >= ($readlength-2))
	   )
	
	#if (($holder[9] >= (($readlength{$holder[1]}*5)-($mis*9))) && ($holder[9] <= ($len * 5)))
	{
		print OUT "$line\n";
	}
			
####################################################################################

		
        $count=0;
}
#my @time=times;
#print OUT1 $time[0]."\n".$time[1]."\n".$time[2]."\n".$time[3]."\n";

close IN;
close IN1;
close OUT;
#close OUT1;
