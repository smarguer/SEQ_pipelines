#!/usr/bin/perl

use strict;
use warnings;

#echo "index 1 - CGTGAT"
#echo "index 2 - ACATCG"
#echo "index 3 - GCCTAA"
#echo "index 4 - TGGTCA"
#echo "index 5 - CACTGT"
#echo "index 6 - ATTGGC"
#echo "index 7 - GATCTG"
#echo "index 8 - TCAAGT"
#echo "index 9 - CTGATC"
#echo "index 10 - AAGCTA"
#echo "index 11 - GTAGCC"
#echo "index 12 - TACAAG"

if (@ARGV != 1) {die "wrong number of files";}
(my $in)=@ARGV;


my $line;
my $count=0;
my @holder;
my @holder1;
my %count;
my %rand;
my $out;
my $out1;
my $index;
my $feature;
my @index=("CGTGAT","ACATCG","GCCTAA","TGGTCA","CACTGT","ATTGGC","GATCTG","TCAAGT","CTGATC","AAGCTA","GTAGCC","TACAAG","Unmapped");
my @features=("ERCC-00130",
"ERCC-00096",
"ERCC-00074",
"ERCC-00002",
"ERCC-00004",
"ERCC-00171",
"ERCC-00113",
"ERCC-00046",
"ERCC-00136",
"ERCC-00108",
"ERCC-00009",
"ERCC-00145",
"ERCC-00003",
"ERCC-00116",
"ERCC-00042",
"ERCC-00111",
"ERCC-00043",
"ERCC-00092",
"ERCC-00060",
"ERCC-00076",
"ERCC-00022",
"ERCC-00095",
"ERCC-00131",
"ERCC-00035",
"ERCC-00044",
"ERCC-00112",
"ERCC-00062",
"ERCC-00025",
"ERCC-00051",
"ERCC-00162",
"ERCC-00071",
"ERCC-00165",
"ERCC-00079",
"ERCC-00019",
"ERCC-00144",
"ERCC-00053",
"ERCC-00084",
"ERCC-00078",
"ERCC-00170",
"ERCC-00148",
"ERCC-00126",
"ERCC-00099",
"ERCC-00054",
"ERCC-00163",
"ERCC-00059",
"ERCC-00154",
"ERCC-00085",
"ERCC-00034",
"ERCC-00157",
"ERCC-00160",
"ERCC-00028",
"ERCC-00150",
"ERCC-00067",
"ERCC-00143",
"ERCC-00039",
"ERCC-00014",
"ERCC-00077",
"ERCC-00033",
"ERCC-00134",
"ERCC-00031",
"ERCC-00058",
"ERCC-00069",
"ERCC-00147",
"ERCC-00109",
"ERCC-00073",
"ERCC-00120",
"ERCC-00040",
"ERCC-00137",
"ERCC-00013",
"ERCC-00097",
"ERCC-00156",
"ERCC-00158",
"ERCC-00164",
"ERCC-00168",
"ERCC-00123",
"ERCC-00104",
"ERCC-00142",
"ERCC-00024",
"ERCC-00016",
"ERCC-00041",
"ERCC-00081",
"ERCC-00017",
"ERCC-00138",
"ERCC-00012",
"ERCC-00086",
"ERCC-00117",
"ERCC-00098",
"ERCC-00061",
"ERCC-00083",
"ERCC-00075",
"ERCC-00057",
"ERCC-00048",
"CH1_bases",
"CH2_bases",
"CH3_bases",
"CH4_bases",
"CH5_bases",
"CH6_bases",
"Unknown");

##
for $feature (@features)
{
 for $index (@index)
  {
   $count{$feature}{$index}{ALL}=0;
   $count{$feature}{$index}{UNI}=0;
  }
}

$feature=0;
$index=0;

##
open (IN, $in) or die 'could not find the input file';
open (OUT1,'>', $in.".ERCCall") or die 'could not find the input file';
open (OUT2,'>', $in.".ERCCuni") or die 'could not find the input file';
##

while ($line=<IN>)
 {
  $count++;
  chomp ($line);
##
  @holder = split (/ /, $line);
  @holder1=split (';', $holder[1]);
  if(grep /$holder1[4]/, @index)
  {
   $index=$holder1[4];
  }
  else
  {
   $index="Unmapped";
  }
##
#print "$index\n";
##
  if(grep /$holder[5]/, @features)
  {
   $feature=$holder[5];
  }
  else
  {
   $feature="Unknown";
  }
##
#print "$feature\n";
##

  $count{$feature}{$index}{ALL}++;
  unless($rand{$index}{$holder1[3]})
  {
   $count{$feature}{$index}{UNI}++;
   $rand{$index}{$holder1[3]}=1;
  }    
 }
print OUT1 "index1-CGTGAT\tindex2-ACATCG\tindex3-GCCTAA\tindex4-TGGTCA\tindex5-CACTGT\tindex6-ATTGGC\tindex7-GATCTG\tindex8-TCAAGT\tindex9-CTGATC\tindex10-AAGCTA\tindex11-GTAGCC\tindex12-TACAAG\tUnmapped\n";

print OUT2 "index1-CGTGAT\tindex2-ACATCG\tindex3-GCCTAA\tindex4-TGGTCA\tindex5-CACTGT\tindex6-ATTGGC\tindex7-GATCTG\tindex8-TCAAGT\tindex9-CTGATC\tindex10-AAGCTA\tindex11-GTAGCC\tindex12-TACAAG\tUnmapped\n";

for $out (@features)
 {
  print OUT1 "$out";
  print OUT2 "$out";
  for $out1 (@index)
  {
   print OUT1 "\t$count{$out}{$out1}{ALL}";
   print OUT2 "\t$count{$out}{$out1}{UNI}";
  }
  print OUT1 "\n";
  print OUT2 "\n";
 }
close IN;

