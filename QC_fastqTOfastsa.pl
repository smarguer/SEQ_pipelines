#!/usr/bin/perl


use strict;
use warnings;

my $line;
my $count=0;
my $count1=0; 
my $count2=0;
my %counter;
my $seq=0;
my $header=0;
my $Average;
my $i;
my $out1=0;
my $out2=0;
my @qual_table;
my @out;
my %average;
my $readaverage=0;
###ILL or Sanger quals###
my $type=33;
if (@ARGV != 2) {die "need 1 file name and a path";}
(my $path, my $in)=@ARGV;

open (IN, $path."$in.fastq") or die 'could not find the input file';
open (OUT, ">", "$in.fasta") or die 'could not open output file';
open (OUT1, ">", "$in.qualtable") or die 'could not open output file'; 
open (OUT2, ">", "$in.qualaverage") or die 'could not open output file';

while ($line=<IN>)
{
	$count++;
	if ($count == 1)
	{
	 chomp ($line);
	 #if (!$counter{$line}){$counter{$line}=0;}
	 #$counter{$line}++;
	 #if ($counter{$line} > 1)
	 #{
	 # $header = ">$line"."_"."$counter{$line}";
	 #}
	 #else
	 #{
	  $header = ">$line";
	  $header =~ s/\s/\;/g;
	 #}
	}
	
	elsif ($count == 2)
	{
	 $seq = $line;
	}
	
	elsif ($count == 4)
	{
	 $count2++;
	 chomp($line);
	 my @char=split(//,$line);
	 my @qual=map(ord,@char);
	 @qual=map($_-$type,@qual);
	
	 for ($i=0; $i<=$#qual; $i++)
	 {
	  $qual_table[$i]{$qual[$i]}++;
	  #print "$qual[$i]\n";
	  #print "$qual_table[$i]{$qual[$i]}\n";
	  #print "@{$qual_table[$i]}\n";
	  $Average+=$qual[$i];
	 }
	
	 print OUT "$header\;";
	 printf OUT ('%.1f',($Average/($#qual+1)));
	 print OUT "\n$seq";
 	 $readaverage=int((($Average/($#qual+1))+0.5));
	 if ($average{$readaverage})
	 {
	 $average{$readaverage}++;
	 }
         
	 else
	 {
	 $average{$readaverage}=1;
	 }
	 
	 #print "$average{$readaverage}\n";
	 


#print OUT1 "@{$qual_table[$i]}\n"; 
	  
	 
	 #print OUT join("\t",@qual);
	 #print OUT "\n";
	 $count=0;
	 $seq=0;
	 $header=0;
	 $Average=0;
        }
}

####print out qual table######################################

for ($i=0; $i<=$#{@qual_table}; $i++)
{
           for ($out1=1;$out1<=40;$out1++)
 	   {
	   $out[$i][($out1-1)]=0;
   	   }
           

	   for ($out1=1;$out1<=40;$out1++)
 	   {
	    if ($qual_table[$i]{$out1})
	    {
             $out[$i][($out1-1)]=$qual_table[$i]{$out1};
            }
           }
           print OUT1 join ("\t",@{$out[$i]})."\n";
}

####print average file###########################################

for ($out2=1;$out2<=40;$out2++)
{
if ($average{$out2})
{
print OUT2 "$out2\t$average{$out2}\n";
}
else
{
print OUT2 "$out2\t0\n";
}
}

#################################################################

#for my $out (sort keys %counter)
#{
#if ($counter{$out} > 1)
#{
#$count1++;
#}
#if ($count1 != 0)
#{
#print "$count1\n";
#}
#}

close IN;
close OUT;
close OUT1;
close OUT2;
