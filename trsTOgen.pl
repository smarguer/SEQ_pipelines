


#!/usr/bin/perl


use strict;
use warnings;
use List::Util qw[min max];

if (@ARGV != 2) {die "need a gff file and a .ALLtrs file";}
(my $gff, my $in)=@ARGV;
my $in1=$gff;

open (IN, $gff) or die 'could not find the gff file';
open (IN1, $in) or die 'could not find the input file';
open (IN2, $in1) or die 'could not find the input file';
open (OUT, ">", "$in.GEN") or die 'could not open output file'; 
my $line=<IN>;
my $line1;
my $line2=<IN2>;
my $count=0;
my $count1=0;
#my $count3=0;
my $out=0;
my $out1=0;
my $out2=0;
my $out3=0;
my $out4=0;
my $out5=0;
my $i=0;
my $j=0;
my $imax=0;
my $jmax=0;
my $out1max=0;
my @holder;
my @holder1;
my @holder2;

my $translength=0;
my $transHashLength=0;
my $length=0;
my $matchS;
my $matchE;
my $featureNumber;
my %GENE;
my %CDS;
my %INTRON;
my %trans;
my %count;
my %count1;
my %INTcount;
my %CDScount;

my $REV;
my %REVtrans;

my @features=('CDS','misc_RNA','tRNA','snRNA','snoRNA','LTR','rRNA','mRNA');
my @excluded=('dgI repeat','dgII repeat','dgIII repeat','dhI repeat', 'dhII repeat','dhIII repeat','CDS_motif','BLASTN_HIT','gap','repeat_unit','mRNA','rep_origin','polyA_signal','LTR','tRNA', 'repeat_region');
####################################################################
##Creates ref unsed to corrects for strand in exon and intron number
#####################################################################

while ($line2=<IN2>)
{
 chomp ($line2);
 @holder2 = split (/\t/, $line2);
 
 if ($holder2[8] eq 'No comment'){next;}

 unless ($CDScount{$holder2[9]})
 {
  $CDScount{$holder2[9]}=0;
 }
 unless ($INTcount{$holder2[9]})
 {
  $INTcount{$holder2[9]}=0;
 }
 if (grep /$holder2[2]/,@features)
 {
 $CDScount{$holder2[9]}++;
 }
 elsif ($holder2[2] eq "intron")
 {
 $INTcount{$holder2[9]}++;
 }

}
close IN2;





#reads gff information.
#print "preparing junctions data structure.\n";
while ($line=<IN>)
{

    chomp ($line);
    @holder= split (/\t/,$line);

    if ($holder[8] eq 'No comment'){next;}

    unless ($count{$holder[9]})
    {
    $count{$holder[9]}=0;
    }
    unless ($count1{$holder[9]})
    {
    $count1{$holder[9]}=0;
    }

    if (grep /$holder[9]/, @excluded)
    {
    next; 
    }

#creates GENE, CDS and INTRON hashes, %{Name}{featureNumber}=featureLength.
    if ($holder[2] eq 'gene')
    {
     $GENE{$holder[9]}{START}=$holder[3];
     $GENE{$holder[9]}{END}=$holder[4];
	 $GENE{$holder[9]}{STRAND}=$holder[6];
	 $GENE{$holder[9]}{CHR}=$holder[12];
    }

    elsif (grep /$holder[2]/,@features)
    {


        $count{$holder[9]}++;
        if ($holder[6] eq "+")
        {
         $featureNumber=$count{$holder[9]};
        }
        elsif ($holder[6] eq "-")
        {
         $featureNumber=abs($count{$holder[9]}-$CDScount{$holder[9]})+1;
        }
 		else
        {
         $featureNumber=$count{$holder[9]};
        }   
        $CDS{$holder[9]}{$featureNumber}=$holder[4]-$holder[3]+1;    

##debugging
#print "CDS\t$holder[9]\t$featureNumber\t$CDS{$holder[9]}{$featureNumber}\t$count{$holder[9]}\t$CDScount{$holder[9]}\n";
##
    }    

 	elsif ($holder[2] eq 'intron')
    {
        $count1{$holder[9]}++;
        if ($holder[6] eq "+")
        {
         $featureNumber=$count1{$holder[9]};
        }
        elsif ($holder[6] eq "-")
        {
         $featureNumber=abs($count1{$holder[9]}-$INTcount{$holder[9]})+1;
        }
    	else
        {
         $featureNumber=$count1{$holder[9]};
        }

        $INTRON{$holder[9]}{$featureNumber}=$holder[4]-$holder[3]-1;

##debugging
#print "INTRON\t$holder[9]\t$featureNumber\t$INTRON{$holder[9]}{$featureNumber}\t$count1{$holder[9]}\t$INTcount{$holder[9]}\n";
##    
    }
}

#create transformation hash.
for $out (keys %GENE)
{
#print "$out\n"; 
 $out1max=max(keys %{$CDS{$out}});

unless($out1max)
{
print "$CDS{$out}{1}\n";
print "$out\n";
}

 $translength=0;
 for ($out1=1;$out1<=$out1max;$out1++ )
 {

#print "$out1\n";
#print "$CDS{$out}{$out1}\n";
  unless(defined $CDS{$out}{$out1})
  {
   print "$out\t$out1\n";
  }
  $translength+=$CDS{$out}{$out1};
#print "$translength\n";
 }
 for ($out2=1;$out2<=$translength;$out2++)
 {
  $trans{$out}{$out2}=$out2;
 }
}

#transforms coordinates from transcript to unspliced transcript.
for $out3 (keys %trans)
{

 $imax=max(keys %{$CDS{$out3}});
 $jmax=max(keys %{$trans{$out3}});

 $length=0;
 for ($i=1;$i<=$imax;$i++)
 {
  $length += $CDS{$out3}{$i}; 
  for ($j=1;$j<=$jmax;$j++)
  {
   #if ($trans{$out3}{$j} >= $length)
    if ($j > $length)
    {
    if ($INTRON{$out3}{$i})
    {
     $trans{$out3}{$j}+=$INTRON{$out3}{$i};
    }
   }
  }
 }
###keep for debugging
#for ($out5=1;$out5<=$jmax;$out5++)
#{
#print "$out5\t$trans{$out3}{$out5}\n";  
#}


}

####################################################
#Processes .ALLtrs file
#####################################################

while ($line1=<IN1>)
{
 chomp ($line1);
 @holder1 = split (/ /, $line1);
 unless ($holder1[0]){next;}
 if ($holder1[0] ne "vulgar:"){next;}

unless(defined %{$trans{$holder1[5]}})
{
print "$line1\n";
next;
} 


if ($holder1[8] eq '-')
{
 $matchS=$holder1[7]+1;
 $matchE=$holder1[6];
}
else
{
 $matchS=$holder1[6]+1;
 $matchE=$holder1[7];
}

#if (($holder1[5] eq 'SPBC651.06')||($holder1[5] eq 'SPAC25H1.10c'))
#{
#next;
#}

######################
#for set up
#unless(defined %{$trans{$holder1[5]}}){next;}
######################

unless($trans{$holder1[5]}{1})
{ 
print "$holder1[5]\tno trans\n";
}
unless($GENE{$holder1[5]}{START} && $GENE{$holder1[5]}{END})
{
print "$holder1[5]\tno GENE\n";
}
unless($trans{$holder1[5]}{$matchS} && $trans{$holder1[5]}{$matchE})
{
if ($holder1[8] eq '-')
{
print "MINUS$holder1[5]\tmissing coord\t$holder1[6]\t$holder1[7]\t".min(keys %{$trans{$holder1[5]}})."\t".max(keys %{$trans{$holder1[5]}})."\n";
}
else
{
print "$holder1[5]\tmissing coord\t$holder1[6]\t$holder1[7]\t".min(keys %{$trans{$holder1[5]}})."\t".max(keys %{$trans{$holder1[5]}})."\n";
}
next;
}

#if ($holder1[8] eq '+')
# {
  if ($GENE{$holder1[5]}{STRAND} eq '+')
  {
  $holder1[6]=$trans{$holder1[5]}{$matchS}+$GENE{$holder1[5]}{START}-1;
  $holder1[7]=$trans{$holder1[5]}{$matchE}+$GENE{$holder1[5]}{START}-1;
  if ($holder1[8] eq '-')
  {
  $holder1[8]='-';
  }
  }
  elsif ($GENE{$holder1[5]}{STRAND} eq '-')
  {

  $holder1[6]=$GENE{$holder1[5]}{END}-$trans{$holder1[5]}{$matchS}+1;
  $holder1[7]=$GENE{$holder1[5]}{END}-$trans{$holder1[5]}{$matchE}+1;

  if ($holder1[8] eq '-')
  {
   $holder1[8]='+';
  }
  else
  {
   $holder1[8]='-';
  }
  }
# }
 #elsif ($holder1[8] eq '-')
 #{
 #antisense transreads
 #}

if ($GENE{$holder1[5]}{CHR}==1){$holder1[5]="CH1_bases";}
elsif ($GENE{$holder1[5]}{CHR}==2){$holder1[5]="CH2_bases";}
elsif ($GENE{$holder1[5]}{CHR}==3){$holder1[5]="CH3_bases";}
elsif ($GENE{$holder1[5]}{CHR}==4){$holder1[5]="CH4_bases";}
elsif ($GENE{$holder1[5]}{CHR}==5){$holder1[5]="CH5_bases";}
elsif ($GENE{$holder1[5]}{CHR}==6){$holder1[5]="CH6_bases";}


print OUT join (" ",@holder1)."\n";
}




close IN;
close IN1;
close OUT;

