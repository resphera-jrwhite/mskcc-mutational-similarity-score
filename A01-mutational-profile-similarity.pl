#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw/shuffle/;
use File::Spec;
#*********************************************************************
# Author: James Robert White PhD
# Email: jwhite@respherabio.com
# Code developed for the Diaz Lab at MSKCC
#*********************************************************************
use Getopt::Std;
use warnings;

use vars qw/$opt_i $opt_j $opt_p $opt_t $opt_e/;

getopts("i:j:p:t:e:");

my $usage = "Usage:  $0 \
        -i file with Sample1 mutations (mutations formatted as: BRAF_chr7_140453136-140453136_A_T)
        -j file with Sample2 mutations (mutations formatted as: BRAF_chr7_140453136-140453136_A_T)
        -t cancer type (see below for options)
        -p panel ID (ASX, PS64, CS125, CS447, MSKIMPACT341, MSKIMPACT410, MSKIMPACT468, CAPPSeq, WEX, Sample1)
        -e list of TCGA ids to exclude from empirical analysis (for internal validation)

.OPTIONS.
  Available cancer types: ALLTCGA,     ACC,    BLCA,  BRCA,  CESC,  CHOL,  COAD,  DLBC,  
                             ESCA,     GBM,    HNSC,  KICH,  KIRC,  KIRP,   LGG,  LIHC, 
                             LUAD,    LUSC,      OV,  PAAD,  PCPG,  PRAD,  READ,  SARC, 
                             SKCM,    STAD,    TGCT,  THCA,  THYM,  UCEC,   UCS,   UVM, COADREAD, 
                         
                          COADMSS,   UCECMSS, COADREADMSS,
                          DIAZMSI, MOSAICMSI, DIAZMSIPOLE
                                  
  Available panel IDs:    ASX, PS64, CS125, CS447, MSKIMPACT341, MSKIMPACT410, 
                          MSKIMPACT468, CAPPSeq, WEX, Sample1

.DESCRIPTION.
  This code is designed to use recurrent mutations mined from TCGA (MC3)
  to predict whether two tumor samples are related in nature. The code
  computes a probability that two random tumor samples overlap by chance.
  The code adjusts for mutations that are highly recurrent as well, e.g.
  the BRAF V600E mutation in Thyroid cancer is present in approximately 60%
  of patients, so it's much more likely to occur by chance than a different
  more unique mutation.
  
  Mutation coordinates are in HG19 space.

  If the user chooses \"Sample1\" as the input panel ID, then the ROIs are limited to those specific mutations
  listed in the Sample1 input file (simulates ddPCR results).

  References for calculations:
  1. Wilson single proportional confidence interval: Newcombe RG (1998). Two-sided confidence 
  intervals for the single proportion: comparison of seven methods. Statistics in Medicine.

\n";

my $abs_path = File::Spec->rel2abs( $0 ) ;
my @abs_path = split /\//, $abs_path;
$abs_path    = join("/", @abs_path[0..($#abs_path-1)]);

die $usage unless defined $opt_i
              and defined $opt_j
              and defined $opt_t
              and defined $opt_p;

my $ifile = $opt_i;
my $jfile = $opt_j;
  # using the variance of the proportion estimate above,
  # compute the variance of product of refval*refval
  # this is E[p^4] - [E[p^2]]^2 

my $type  = $opt_t;
my $panel = $opt_p;

my %exclude = ();
if (defined($opt_e)){
  my $exclude = $opt_e;
  my @exclude = split /\,/, $exclude;
  for my $e (@exclude){
    $exclude{$e} = 1;	
  } 
}
my %panels = ('ASX'=>1, 'PS64'=>1, 'CS125'=>1, 'CS447'=>1, 'MSKIMPACT341'=>1, 'MSKIMPACT410'=>1, 'MSKIMPACT468'=>1, 'CAPPSeq'=>1, 'WEX'=>1, 'Sample1'=>1);

my %types = (
'ALLTCGA'     => 'All_cancers', 
'ACC'         => 'Adrenal',   
'BLCA'        => 'Bladder',
'BRCA'        => 'Breast',    
'CESC'        => 'Cervical',  
'CHOL'        => 'Cholangiocarcinoma',
'COAD'        => 'Colon',    
'DLBC'        => 'Diffuse_Large_B-cell_Lymphoma',  
'ESCA'        => 'Esophageal',
'GBM'         => 'Glioblastoma',     
'HNSC'        => 'Head_and_neck',  
'KICH'        => 'Kidney_Chromophobe',
'KIRC'        => 'Kidney_Renal_Clear_Cell',    
'KIRP'        => 'Kidney_Renal_Papillary_Cell',  
'LGG'         => 'Lower_grade_glioma',
'LIHC'        => 'Liver',    
'LUAD'        => 'Lung_adenocarcinoma',  
'LUSC'        => 'Lung_squamous_cell',
'OV'          => 'Ovarian',      
'PAAD'        => 'Pancreatic',  
'PCPG'        => 'Pheochromocytoma_Paraganglioma',
'PRAD'        => 'Prostate',    
'READ'        => 'Rectal',  
'SARC'        => 'Sarcoma',
'SKCM'        => 'Melanoma',   
'STAD'        => 'Stomach', 
'TGCT'        => 'Testicular',
'THCA'        => 'Thyroid', 
'THYM'        => 'Thymoma',  
'UCEC'        => 'Uterine_Corpus_Endometrial',
'UCS'         => 'Uterine_Carcinosarcoma',     
'UVM'         => 'Uveal_Melanoma',
'DIAZMSI'     => 'MSI Cases | Luis Diaz Jr. MD',
'MOSAICMSI'   => 'MSI Cases | MOSAIC algorithm',
'COADREAD'    => 'Colorectal cancers',
'COADMSS'     => 'Colon cancer MSS',    
'UCECMSS'     => 'Uterine_Corpus_Endometrial MSS',
'COADREADMSS' => 'Colorectal MSS', 
'DIAZMSIPOLE' => 'All_cancer_types MSI+POLE | Luis Diaz Jr. MD');

if (!defined($panels{$panel})){
  die "Error: The panel $panel is not available.\n";
}

if (!defined($types{$type})){
  die "Error: The cancer type $type is not available.\n";
}

if (! -e  "$ifile"){
  die "Error: Cannot locate the file $ifile\n";
}
if (! -e  "$jfile"){
  die "Error: Cannot locate the file $jfile\n";
}

my %muts        = ();
my $imuts       = 0;
my %uniqimuts   = ();
my %sample1muts = ();
open IN, "$ifile" or die;
while(<IN>){
  chomp($_);
  if ($_ =~ /\ /){
    die "Error: No white space permitted in the input files.\n";
  }
  if ($_ =~ /\t/){
    die "Error: No white space permitted in the input files.\n";
  }
  $muts{$_}++;
  $imuts++;
  $uniqimuts{$_}   = 1;
  $sample1muts{$_} = 1;
}
close IN;

my $jmuts = 0;
open IN, "$jfile" or die;
while(<IN>){
  chomp($_);
  if ($_ =~ /\ /){
    die "Error: No white space permitted in the input files.\n";
  }
  if ($_ =~ /\t/){
    die "Error: No white space permitted in the input files.\n";
  }
  $muts{$_}++;
  $jmuts++;
}
close IN;

my %recurrent = ();
my $intersect = 0;
foreach my $k1 (sort keys %muts){
  if ($muts{$k1} > 1){
    $recurrent{$k1} = 1;
    $intersect++;
  }
}

my %ref         = ();
my %numerator   = ();
my %denominator = ();
my $total  = 0;
my $refid  = $panel;
if ($panel eq "WEX" or $panel eq "Sample1"){
  $refid = "ALL";
}

# *******************************************************
# load MSI sets if needed
# *******************************************************
my %CUSTOM = ();
if ($type =~ /^(DIAZMSI|MOSAICMSI|DIAZMSIPOLE|COADMSS|UCECMSS|COADREADMSS)$/){
  open IN, "$abs_path/supporting-materials/background_dbs/$type.samples.txt" or die "Cannot find bg db list: $abs_path/supporting-materials/background_dbs/$type.samples.txt";
  while(<IN>){
    chomp($_);
    $CUSTOM{$_} = 1;  
  }
  close IN;
}


# *******************************************************
# load TCGA empirical control set
# *******************************************************
my %tcgaobs       = ();
my %tcgasamples   = ();
my %tcgamutcounts = (); 
open IN, "$abs_path/supporting-materials/background_dbs/stage2.output.$refid.txt" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_; 
  # TCGA-02-0003-01A-01D-1490-08	TCGA-02-0003	GBM	TP53	17	7577094	7577094	Missense_Mutation	G	A	60	46	43.4	c.844C>T	p.R282W
  # TCGA-02-0003-01A-01D-1490-08	TCGA-02-0003	GBM	TP53	17	7578396	7578396	Missense_Mutation	G	T	78	47	37.6	c.534C>A	p.H178Q
  my $id                   = $A[1];
  my $idtype               = $A[2];
  
  if (defined($exclude{$id})){ # exclude the sample?
    next;	
  } 

  if ($type =~ /^(DIAZMSI|MOSAICMSI|DIAZMSIPOLE|COADMSS|UCECMSS|COADREADMSS)$/){
    if (!defined($CUSTOM{$id})){ # skip it b/c it's not in the custom list
      next;	
    }
  }elsif($type eq "COADREAD" and $idtype !~ /^(COAD|READ)$/){
  	next;
  }elsif($type ne "ALLTCGA" and $type ne "COADREAD" and $idtype ne $type){ # skip this sample
    next;	
  }    
  
  my $mut                  = "$A[3]\_chr$A[4]_$A[5]-$A[6]\_$A[8]\_$A[9]";
  $tcgasamples{$id}        = 1;

  if ($panel eq "Sample1" and !defined($sample1muts{$mut})){
  	# the reference background is only the Sample1 mutations itself
    next;
  }

  $tcgaobs{$id}{$mut}      = 1;
  $tcgamutcounts{$mut}++;
}
close IN;
my @sortedtcgasamples = sort keys %tcgasamples;
$total = $#sortedtcgasamples+1;

# *******************************************************
# if we don't have any overlap return now
# *******************************************************
if ($intersect == 0){
  my $relatedness = sprintf("%3.3f",100*1/$total);
  print "BEGIN REPORT:\n";
  print "Input Cancer Type: $type \($types{$type}\)\n";
  print "Input Panel ID: $panel\n";
  print "Sample1: $ifile\n";
  print "Sample2: $jfile\n";
  print "Negative control samples: $total\n";
  print "Total mutations in Sample1: $imuts\n";
  print "Total mutations in Sample2: $jmuts\n";
  print "Total mutations Shared: $intersect\n";
  print "\n";
  print "No recurrent mutations found.\n"; 
  print "Probability samples are related < $relatedness\%\n";
  exit;
}

# otherwise, we have recurrent mutations to log
print "STATUS: Computing background frequencies of each Shared mutation...\n";
foreach my $k1 (sort keys %recurrent){
  if (defined($tcgamutcounts{$k1})){
    $ref{$k1}         = $tcgamutcounts{$k1}/$total;	
    $numerator{$k1}   = $tcgamutcounts{$k1};
    $denominator{$k1} = $total;
  }else{
    $ref{$k1}         = 1/$total;
    $numerator{$k1}   = 1;
    $denominator{$k1} = $total;  	
  }
}

# empirical assessment
print "STATUS: Empirical assessment comparing Sample1 to negative controls...\n";
my @simshared  = ();
my $greaterthanorequaltointersect = 0;
my $iterations = $total;
for my $ti (1 .. $iterations){
  my $kti      = $sortedtcgasamples[($ti-1)];
  
  # compute shared between Sample1 and this sample:
  my $ktishared = 0;
  foreach my $ktimut (sort keys %{$tcgaobs{$kti}}){
    if (defined($uniqimuts{$ktimut})){
      $ktishared++;  	
    }	
  }

  if ($ktishared >= $intersect){
    $greaterthanorequaltointersect++;  	
  }
  push @simshared, $ktishared;
}

print "\nEMPIRICAL RESULTS: $greaterthanorequaltointersect out of $iterations greater than or equal to observed value of $intersect Shared mutations...\n\n";

# for each recurrent mut we now have reference background frequencies
# BEGIN REPORT
print "BEGIN REPORT:\n";
print "Input Cancer Type: $type \($types{$type}\)\n";
print "Input Panel ID: $panel\n";
print "Sample1: $ifile\n";
print "Sample2: $jfile\n";
print "Negative control samples: $total\n";
print "Total mutations in Sample1: $imuts\n"; 
print "Total mutations in Sample2: $jmuts\n";
print "Total mutations Shared: $intersect\n";
print "\n";
print "SHARED MUTATION STATISTICS:\n";

my $format = "%-30s\t%-8s\t%s\n";
my $line   = sprintf($format, "Recurrent_Mutation", "General_Frequency_in_$type\(%)", "95%CI");
print "Shared_Mutation\tGeneral_Frequency_in_$type\(%)\t95%CI\n";

foreach my $k1 (sort {$ref{$b} <=> $ref{$a}} keys %ref){
  my $refprint  = sprintf("%3.3f", 100*$ref{$k1});
  my $p = $ref{$k1}; 
  my $n = $denominator{$k1};

  # Wilson score interval with continuity correction
  # Newcombe, R. G. (1998). "Two-sided confidence intervals for the single proportion: comparison of seven methods". Statistics in Medicine. 17 (8): 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. PMID 9595616.
  my ($wlow,$whi) = wilson($p, $n, 1.96);

  $wlow *= 100;
  $whi  *= 100;
  $whi   = sprintf("%3.3f", $whi);
  $wlow  = sprintf("%3.3f", $wlow);
  
  my $one    = $k1;
  my $two    = "$refprint\%";
  my $three  = "($wlow\%,$whi\%)";
  my $line   = sprintf($format, $one, $two, $three);
  print "$line";
  # print "$k1\t$refprint\% 95%CI:($wlow\%,$whi\%)\n";
}

if ($greaterthanorequaltointersect == 0){
  $greaterthanorequaltointersect = 1;	
}
my $type1errorestimate = sprintf("%3.3f", 100*$greaterthanorequaltointersect/$iterations);
my ($elow,$ehi) = wilson($greaterthanorequaltointersect/$iterations, $iterations, 1.96);
$elow *= 100;
$ehi  *= 100;
$ehi   = sprintf("%3.3f", $ehi);
$elow  = sprintf("%3.3f", $elow);
print "\nNEGATIVE CONTROL EMPIRICAL RESULTS:\n";
print "TYPE1 ERROR ESTIMATE: $type1errorestimate\% 95%CI:($elow\%,$ehi\%)\n";

my $ovp     = sprintf("%3.3f", 100-$type1errorestimate);
my $ovehi   = sprintf("%3.3f", 100-$elow);
my $ovelow  = sprintf("%3.3f", 100-$ehi);
print "OVERALL PROBABILITY OF RELATEDNESS: $ovp\% 95%CI:($ovelow\%,$ovehi\%)\n";

# Wilson score interval with continuity correction
# Newcombe, R. G. (1998). "Two-sided confidence intervals for the single proportion: comparison of seven methods". Statistics in Medicine. 17 (8): 857–872. doi:10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E. PMID 9595616.
sub wilson
{
  my ($p,$n,$z) = @_;

  my $wminus = 0;
  my $wplus  = 1;

  my $cminus = (2*$n*$p + $z**2 - ($z*sqrt($z**2 - (1/$n) + 4*$n*$p*(1-$p) + (4*$p-2)) + 1)) / (2*($n+$z**2));
  my $cplus  = (2*$n*$p + $z**2 + ($z*sqrt($z**2 - (1/$n) + 4*$n*$p*(1-$p) - (4*$p-2)) + 1)) / (2*($n+$z**2));

  if ($cminus > $wminus){
    $wminus = $cminus;
  }
  if ($cplus < $wplus){
    $wplus = $cplus;
  }
  return($wminus,$wplus);
}


