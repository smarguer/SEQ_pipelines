#!/usr/bin/perl
#!/usr/bin/perl


use strict;
use warnings;


########################
if (@ARGV != 2) {die "need 1 file name and a read length";}
(my $in, my $readlength)=@ARGV;

open (IN, $in) or die 'could not find the input file';
#open (IN1, $in) or die 'could not find the chromosome file';
#open (IN2, $reads) or die 'could not find the chromosome file';
open (OUT, ">", "$in.diagnostic") or die 'could not open output file';
#open (OUT2, ">", $out2) or die 'could not open second output file';

########################


my $line;
my $line1;
my $line2;
my $count=0;
my $count1=0;
my $count2=0;
my $count3=0;
my $count4=0;
my $count5=0;
my $count6=0;
my $count7=0;
my $count8=0;
my $count9=0;
my $count10=0;
my $count11=0;
my $count12=0;
my $count13=0;
my $count14=0;
my $count15=0;
my $count16=0;
my $count17=0;
my $count18=0;
my $count19=0;
my $chrom;
my $length=0;
my $holder;
my $match;

my $max=0;
my $min=37;


my @holder;
my @read;
my %reads;
my %partialmatch;
my %endmatch;

#####################
while ($line=<IN>)
{
        
        
        chomp ($line);

        @holder = split (/ /, $line);
        unless ($holder[0]){next;}
        if ($holder[0] !~ /vulgar/){next;}
        
        $count++;


#########counts reads and multiple matches#################
	if ($reads{$holder[1]})
        {
        $reads{$holder[1]}++;
        #print "$holder[1]\n";
        }
        else 
        {
        $count2++;
        $reads{$holder[1]}=1;
        }
########################################################### 



	if ($holder[9])
        {
        
	$match = ($holder[3]-$holder[2]);
	if($match == $readlength)
	{
	 $count3++;
	}

	if($match < $readlength)
	 {
	  $count4++;
          if($partialmatch{$match})
	  {
	   $partialmatch{$match}++;
          }
          else
          {
  	   $partialmatch{$match}=1; 
	  }
         
	  if($endmatch{$holder[2]}{$holder[3]})
	  {
	   $endmatch{$holder[2]}{$holder[3]}++;
  	  }
	  else
	  {
	   $endmatch{$holder[2]}{$holder[3]}=1;
          }
	 }
	#$count6++;
        #if ($match == 50){$count10++;}
	#if ($holder[9] == 250){$count11++;}
	#if ($holder[9] == 241){$count12++;}
	#if ($holder[9] == 232){$count13++;}
	#if ($holder[9] == 223){$count14++;}
	#if ($holder[9] == 214){$count15++;}
	#if ($holder[9] == 205){$count16++;}
	#if (($match >= 40) && ($match < 50)){$count17++;}
	#if (($match >= 30) && ($match < 40)){$count18++;}
	#if (($match >= 20) && ($match < 30)){$count19++;}
	#if ($holder[9] == (5*(abs($holder[3]-$holder[2])))) {$count3++;}
        #if ($holder[9] == ((5*(abs($holder[3]-$holder[2])))-9)) {$count4++;}
        #if ($holder[9] == ((5*(abs($holder[3]-$holder[2])))-18)) {$count5++;}
        #if ($holder[9] == ((5*(abs($holder[3]-$holder[2])))-27)) {$count8++;}
        #if ($holder[9] == ((5*(abs($holder[3]-$holder[2])))-36)) {$count9++;}
        
	#if (abs($holder[3]-$holder[2])>$max){$max=abs($holder[3]-$holder[2]);}
        #if (abs($holder[3]-$holder[2])<$min){$min=abs($holder[3]-$holder[2]);}
        #if (abs($holder[3]-$holder[2])>=20){$count7++;}
        #$length=$length+abs($holder[3]-$holder[2]);
        }
}


#########################Counts unique reads#################################################
for my $out (keys %reads)
{
 if ($reads{$out}==1)
 {
  $count1++;
 }
 else
 {
  $count5++;
 }
}
#############################################################################################



print OUT "total: ".$count."\n";
print OUT "reads: ".$count2."\n";
print OUT "Reads with unique match:     ".$count1."\n";
print OUT "Reads with multiple matches: ".$count5."\n";
print OUT "***********\n";
print OUT "full length $readlength match: ".$count3."\n";
print OUT "partial match:           ".$count4."\n";
#################
for my $out (sort keys %partialmatch)
{
 print OUT "$out nt partial match: $partialmatch{$out}\n";
}
################
for my $out5 (sort keys %endmatch)
{
 for my $out3 (sort keys %{$endmatch{$out5}})
 {
  print OUT "5': $out5 nt, 3': $out3 -> $endmatch{$out5}{$out3}\n";
 }
}









#print "\t0mis:".$count11."\n";
#print "\t1mis:".$count12."\n";
#print "\t2mis:".$count13."\n";
#print "\t3mis:".$count14."\n";
#print "\t4mis:".$count15."\n";
#print "\t5mis:".$count16."\n";
#print "***********\n";
#print "40nt to 49nt long match: ".$count17."\n";
#print "30nt to 39nt long match: ".$count18."\n";
#print "20nt to 29nt long match: ".$count19."\n";
#print "***********\n";
#print "perfect partial match: ".$count3."\n";
#print "\t1mis partial match: ".$count4."\n";
#print "\t2mis partial match: ".$count5."\n";
#print "\t3mis partial match: ".$count8."\n";
#print "\t4mis partial match: ".$count9."\n";
#print "total mapped: ".$count6."\n";
#print "matches >= 20nt: ".$count7."\n";
#print "max match: ".$max."\n";
#print "min match: ".$min."\n";
#print "average match length: ".$length/$count6."\n";

close IN;
close OUT;
