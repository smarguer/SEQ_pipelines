use strict;
use warnings;
use POSIX;
use Storable qw(dclone);

########################

if (@ARGV != 2) {die "need a file and a gff file";}
(my $in, my $gff)=@ARGV;
my $in1=$gff;

open (IN, $gff) or die 'could not find the input file';
open (IN2, $in1) or die 'could not find the input file';

$gff =~ /(gff_\d{6}\.txt)/;
my $gff1=$1;
#print "$gff1\n";
#die;
########################

my $line=<IN>;
my $line1;
my $line2=<IN2>;
my $mid;
my $gene;
my $strand;
my $genestrand1;
my $genestrand2;
my $count=0;
my $count1=0;
my $check=0;
my @holder;
my @holder1;
my @holder2;
my $name;
my $length;
my $seg;
my $countmatch;
my $countstrand=1;
my $matchStart;
my $matchEnd;
my $start;
my $end;
my $i;
my $what;

my %gffStart;
my %gffEnd;
my %gffChr;
my %gffStrand;
my %match;
my %out;
my $chr;
my $featureNumber;
my $first;
my $last;
my $out;
my $out2;
my %HIT;
my %HITstart;
my %HITend;
my %HITmid;
my %HITfirst;
my %HITlast;
my %FEATURE;
my %GENE;
my %FEATUREgff;
my @done;
my @done1;

my %sense;
my %antisense;
my %count;
my %count1;
my %CDScount;
my %INTcount;
my @features=('CDS','misc_RNA','tRNA','snRNA','snoRNA','LTR','rRNA','mRNA');
#my @excluded=('CDS_motif','BLASTN_HIT','gap','repeat_unit','mRNA','rep_origin','polyA_signal','LTR','tRNA');
my @excluded=('CDS_motif','BLASTN_HIT','gap','repeat_unit','mRNA','promoter','rep_origin','polyA_signal','polyA_site','3\'UTR','5\'UTR','tRNA');
#my @excluded=('CDS_motif','BLASTN_HIT','gap','repeat_unit','mRNA','rep_origin','polyA_signal','3\'UTR','5\'UTR','tRNA');
my %multiple;
my @out;
my @toDo;
my %CHRF;
my %CHRG;
my %TRSF;
my %TRSG;




######################################################################
#Creates ref unsed to corrects for strand in exon and intron numbering
######################################################################

while ($line2=<IN2>)
{
 chomp ($line2);
 @holder2 = split (/\t/, $line2);
 if (grep /$holder2[9]/, @excluded)
 {
  next; 
 } 	
 
 if ($holder2[9] eq "LTR")
 {
  $holder2[9] = "LTR.$holder2[3].$holder2[12]";
 }

 unless ($CDScount{$holder2[9]})
 {
  $CDScount{$holder2[9]}=0;
 }
 unless ($INTcount{$holder2[9]})
 {
  $INTcount{$holder2[9]}=0;
 }
 if (grep /$holder2[2]/, @features)
 {
 $CDScount{$holder2[9]}++;
 }
 elsif ($holder2[2] eq "intron")
 {
 $INTcount{$holder2[9]}++;
 }

}
close IN2;

#############################################################
##Reads gff and creates data structures for segment mapping
#############################################################

while ($line=<IN>)
{
 chomp ($line);
 @holder = split (/\t/, $line);
 if (grep /$holder[9]/, @excluded)
 {
  next; 
 } 
 
 if ($holder[9] eq "LTR")
 {
  $holder[9] = "LTR.$holder[3].$holder[12]";
 }

 unless ($count{$holder[9]})
 {
  $count{$holder[9]}=0;
 }
 unless ($count1{$holder[9]})
 {
  $count1{$holder[9]}=0;
 }      
   
 if ($holder[2] eq "gene")
 {
  $gffStart{$holder[9]}=$holder[3];
  $gffEnd{$holder[9]}=$holder[4];
  $gffStrand{$holder[9]}=$holder[6];
  $gffChr{$holder[9]}=$holder[12];
  $count1=0;
  $GENE{$holder[9]}{START}=$holder[4];
  $GENE{$holder[9]}{END}=$holder[3];
  $GENE{$holder[9]}{SEG}=0;  
  $GENE{"AS.$holder[9]"}{START}=$holder[4];
  $GENE{"AS.$holder[9]"}{END}=$holder[3];
  $GENE{"AS.$holder[9]"}{SEG}=0;
 }
#####
#####
 elsif (grep /$holder[2]/, @features)
 {
  $count{$holder[9]}++;
  if ($holder[6] eq "+")
  {
   $sense{$holder[9]} = "+";
   $antisense{$holder[9]} = "-";
   $featureNumber=$count{$holder[9]};
  }
  elsif ($holder[6] eq "-")
  {
   $sense{$holder[9]} = "-";
   $antisense{$holder[9]} = "+";
   $featureNumber=abs($count{$holder[9]}-$CDScount{$holder[9]})+1;
  }
   $FEATURE{"$holder[9].E$featureNumber"}{SEG}=0;
   $FEATURE{"$holder[9].E$featureNumber"}{START}=$holder[4];
   $FEATURE{"$holder[9].E$featureNumber"}{END}=$holder[3];
   
   $FEATURE{"AS.$holder[9].E$featureNumber"}{SEG}=0;      
   $FEATURE{"AS.$holder[9].E$featureNumber"}{START}=$holder[4];      
   $FEATURE{"AS.$holder[9].E$featureNumber"}{END}=$holder[3];      
   
   $FEATUREgff{"$holder[9].E$featureNumber"}{SEG}=0;
   $FEATUREgff{"$holder[9].E$featureNumber"}{START}=$holder[3];
   $FEATUREgff{"$holder[9].E$featureNumber"}{END}=$holder[4];
   
   $FEATUREgff{"AS.$holder[9].E$featureNumber"}{SEG}=0;      
   $FEATUREgff{"AS.$holder[9].E$featureNumber"}{START}=$holder[3];      
   $FEATUREgff{"AS.$holder[9].E$featureNumber"}{END}=$holder[4];      
}
 
 elsif ($holder[2] eq "intron")
 {
  $count1++;
  $name="INTRON_".$holder[9]."_".$count1;
  $gffStart{$name}=$holder[3];
  $gffEnd{$name}=$holder[4];
  $gffStrand{$name}=$holder[6];
  $gffChr{$name}=$holder[12];       
  $count1{$holder[9]}++;
 
  if ($holder[6] eq "+")
  {
   $sense{$holder[9]} = "+";
   $antisense{$holder[9]} = "-";
   $featureNumber=$count1{$holder[9]};
  }
  elsif ($holder[6] eq "-")
  {
   $sense{$holder[9]} = "-";
   $antisense{$holder[9]} = "+";
   $featureNumber=abs($count1{$holder[9]}-$INTcount{$holder[9]})+1;
  }

  $FEATURE{"INTRON_$holder[9].I$featureNumber"}{SEG}=0;
  $FEATURE{"INTRON_$holder[9].I$featureNumber"}{START}=$holder[4];
  $FEATURE{"INTRON_$holder[9].I$featureNumber"}{END}=$holder[3];
  
  $FEATURE{"AS.INTRON_$holder[9].I$featureNumber"}{SEG}=0;
  $FEATURE{"AS.INTRON_$holder[9].I$featureNumber"}{START}=$holder[4];
  $FEATURE{"AS.INTRON_$holder[9].I$featureNumber"}{END}=$holder[3];
  
  $FEATUREgff{"INTRON_$holder[9].I$featureNumber"}{SEG}=0;
  $FEATUREgff{"INTRON_$holder[9].I$featureNumber"}{START}=$holder[3];
  $FEATUREgff{"INTRON_$holder[9].I$featureNumber"}{END}=$holder[4];
  
  $FEATUREgff{"AS.INTRON_$holder[9].I$featureNumber"}{SEG}=0;
  $FEATUREgff{"AS.INTRON_$holder[9].I$featureNumber"}{START}=$holder[3];
  $FEATUREgff{"AS.INTRON_$holder[9].I$featureNumber"}{END}=$holder[4];
  
  $featureNumber=0;
 }
}

#%FEATUREgff = dclone (\%FEATURE);
#%FEATUREgff = %{dclone (\%FEATURE) };
#%{$FEATUREgff}=%{$FEATURE};


%CHRF=%{dclone(\%FEATURE)};
%TRSF=%{dclone(\%FEATURE)};
%CHRG=%{dclone(\%GENE)};
%TRSG=%{dclone(\%GENE)};
%FEATURE=();
%GENE=();

print "gff parsed\n";
#print "#####start#####\n";
#print "$FEATURE{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$CHRF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$TRSF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$GENE{'SPAC5H10.06c'}{SEG}\n";
#print "$CHRG{'SPAC5H10.06c'}{SEG}\n";
#print "$TRSG{'SPAC5H10.06c'}{SEG}\n";

#die;
#####################################
#Reads .map files and maps reads 
#####################################

@toDo=`ls $in*.5mis.ALL*.col.map*$gff1`;
#print "@toDo\n";

for (@toDo){$_=~s/^\s+|\s+$//g};

#print "To process: \n";
#print join("\n",@toDo)."\n";


foreach (@toDo)
{
 $what=$_;
 open (IN1, $what) or die "could not find the input file named $_";
 $line1=<IN1>;
 print "processing $what\n";
 if ($what =~ /chr/)
 {
  %FEATURE=%{dclone(\%CHRF)};
  %GENE=%{dclone(\%CHRG)};
  #print "processing ALLchr...\n";
 }
 elsif ($what =~ /trs/)
 {
  %FEATURE=%{dclone(\%TRSF)};
  %GENE=%{dclone(\%TRSG)};
  #print "processing ALLtrs...\n";
 }
 else
 {
  die 'Not quite sure about this file...\n';
 }


#print "#####in loop start#####\n";
#print "$FEATURE{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$CHRF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$TRSF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$GENE{'SPAC5H10.06c'}{SEG}\n";
#print "$CHRG{'SPAC5H10.06c'}{SEG}\n";
#print "$TRSG{'SPAC5H10.06c'}{SEG}\n";


 while ($line1=<IN1>)
 {
  @done=();
  @done1=();
  @holder1=();
  $count++;
  chomp ($line1);

  @holder1 = split (/\t/, $line1);
 
  my $size=$#holder1;

  for ($i=5;$i<=$size;$i++)
  {
   if ($holder1[$i] eq "NEW")
   {
    next;
   } 
  
   unless (defined $FEATURE{$holder1[$i]}{SEG})
   {
    print "$line1\n$holder1[$i] wasn't included in the data structure\n";
    next;
   }   

#######
#for gene table
#######
   $gene=$holder1[$i];
  
   if ($gene=~/INTRON_/)
   {
    $gene=$`.$';
   }
  
   if ($gene=~/\.I/)
   {
    $gene=$`;
   }    
   elsif ($gene=~/\.E/)
   {
    $gene=$`;
   }
####

   unless (grep /$holder1[$i]/, @done)
   { 
    $FEATURE{$holder1[$i]}{SEG}++;
    push (@done, $holder1[$i]);  
   }

   unless (grep /$gene/, @done1)
   {
    $GENE{$gene}{SEG}++;
    push (@done1, $gene);
   }
  #if (($FEATURE{$holder1[$i]}{START} eq "NON_UNIQUE")||($GENE{$gene}{START} eq "NON_UNIQUE")||($holder1[1] eq "MULTIPLE"))
  #if (defined $multiple{$holder1[0]})
   if ($holder1[0] =~ /M@/)
   {
    $FEATURE{$holder1[$i]}{START} = "NON_UNIQUE";
    $FEATURE{$holder1[$i]}{END} = "NON_UNIQUE";
    $GENE{$gene}{START} = "NON_UNIQUE";
    $GENE{$gene}{END} = "NON_UNIQUE";
   }
  
   if ($FEATURE{$holder1[$i]}{START} ne "NON_UNIQUE")
   {
    if ($FEATURE{$holder1[$i]}{START} > $holder1[1])
    {
     $FEATURE{$holder1[$i]}{START} = $holder1[1];
    }
  
    if ($FEATURE{$holder1[$i]}{END} < $holder1[2])
    { 
     $FEATURE{$holder1[$i]}{END} = $holder1[2];
    }
   }  
   if($GENE{$gene}{START} ne "NON_UNIQUE")
   {
    if ($GENE{$gene}{START} > $holder1[1])
    { 
     $GENE{$gene}{START} = $holder1[1];
    }
  
    if ($GENE{$gene}{END} < $holder1[2])
    {
     $GENE{$gene}{END} = $holder1[2];
    }
   }
  }
#$multiple{$holder1[0]}=1;
 }
#print "#####in loop end I#####\n";
#print "$FEATURE{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$CHRF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$TRSF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$GENE{'SPAC5H10.06c'}{SEG}\n";
#print "$CHRG{'SPAC5H10.06c'}{SEG}\n";
#print "$TRSG{'SPAC5H10.06c'}{SEG}\n";


 if ($what =~ /chr/)
 {
  %CHRF=%{dclone(\%FEATURE)};
  %CHRG=%{dclone(\%GENE)};
  #print "done with ALLchr...\n";
 }
 elsif ($what =~ /trs/)
 {
  %TRSF=%{dclone(\%FEATURE)};
  %TRSG=%{dclone(\%GENE)};
  #print "done with ALLtrs...\n";
 }
 else
 {
  die 'Not quite sure about this file (at the end)...\n';
 }

%FEATURE=();
%GENE=();

close IN1;
#print "#####in loop end II#####\n";
#print "$FEATURE{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$CHRF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$TRSF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$GENE{'SPAC5H10.06c'}{SEG}\n";
#print "$CHRG{'SPAC5H10.06c'}{SEG}\n";
#print "$TRSG{'SPAC5H10.06c'}{SEG}\n";


}

#print "#####after loop#####\n";
#print "$FEATURE{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$CHRF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$TRSF{'SPAC5H10.06c.E1'}{SEG}\n";
#print "$GENE{'SPAC5H10.06c'}{SEG}\n";
#print "$CHRG{'SPAC5H10.06c'}{SEG}\n";
#print "$TRSG{'SPAC5H10.06c'}{SEG}\n";


#die;

####################
#print out
####################

@out=("chr","trs");

foreach (@out)
{

 %FEATURE=();
 %GENE=();

 if($_ eq "chr")
 {
  open (OUT, ">", "$in.5mis.ALLchr.col.mapIV.$gff1.READSFEATURETABLE") or die 'could not open output file';
  open (OUT1, ">", "$in.5mis.ALLchr.col.mapIV.$gff1.READSGENETABLE") or die 'could not open output file';
  %FEATURE=%{dclone(\%CHRF)};
  %GENE=%{dclone(\%CHRG)};
 }
 elsif($_ eq "trs")
 {
  open (OUT, ">", "$in.5mis.ALLtrs.GEN.col.mapV.$gff1.READSFEATURETABLE") or die 'could not open output file';
  open (OUT1, ">", "$in.5mis.ALLtrs.GEN.col.mapV.$gff1.READSGENETABLE") or die 'could not open output file';
  %FEATURE=%{dclone(\%TRSF)};
  %GENE=%{dclone(\%TRSG)};
 } 
 else
 {
  die 'Not sure what to print out...';
 }


 print OUT "Name\tFEAT.start\tFEAT.end\treads\tstart\tend\tAS.reads\tAS.start\tAS.end\n";
 print OUT1 "Name\tchr\tstrand\tORF.start\tORF.end\treads\tstart\tend\tAS.reads\tAS.start\tAS.end\n";

 for $out (keys %FEATURE)
 {
  if ($out=~/^AS/)
  {
   next;
  }
  my $out1="AS.$out";
  if ($out =~ /INTRON_/)
  {
   $out2=$';
  }
  else
  {
   $out2=$out;
  }
  if ($FEATURE{$out1}{SEG}==0)
  {
   $FEATURE{$out1}{START}=0;
   $FEATURE{$out1}{END}=0;
  }
  if ($FEATURE{$out}{SEG}==0)
  {
   $FEATURE{$out}{START}=0;
   $FEATURE{$out}{END}=0;
  }
  print OUT "$out2\t$FEATUREgff{$out}{START}\t$FEATUREgff{$out}{END}\t$FEATURE{$out}{SEG}\t$FEATURE{$out}{START}\t$FEATURE{$out}{END}\t$FEATURE{$out1}{SEG}\t$FEATURE{$out1}{START}\t$FEATURE{$out1}{END}\n";
 }

 for $out (keys %GENE)
 {
  if ($out=~/AS/)
  {
   next;
  }
  my $out1="AS.$out";
  if ($GENE{$out1}{SEG}==0)
  {
   $GENE{$out1}{START}=0;
   $GENE{$out1}{END}=0;
  }
 
  if ($GENE{$out}{SEG}==0)
  {
   $GENE{$out}{START}=0;
   $GENE{$out}{END}=0;
  }



  print OUT1 "$out\t$gffChr{$out}\t$gffStrand{$out}\t$gffStart{$out}\t$gffEnd{$out}\t$GENE{$out}{SEG}\t$GENE{$out}{START}\t$GENE{$out}{END}\t$GENE{$out1}{SEG}\t$GENE{$out1}{START}\t$GENE{$out1}{END}\n";
 }

 close OUT;
 close OUT1;
}
