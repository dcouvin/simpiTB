#!/usr/bin/perl
use strict;
#use warnings;

#Examples of use: 
#perl simpiTB.pl <options> *.fasta
#perl simpiTB.pl <options> *.fastq
#perl simpiTB.pl <options> *.fastq.gz
#perl simpiTB.pl <options> list.txt   
#(please note that the txt file must finish by extension ".txt", and it must contain one accession perl line such as the following example)
#SRR13015794 
#SRR13015795
#accessions could also be separated by commas as follows: 
#perl simpiTB.pl SRR13015794,SRR13015795,SRR13015796

#variables
my $seqfasta = "";
my $seqfastq = "";
my $recap_total_seq = "reconciledTB_list.xls"; # global summary file
my $spollineages_entries = "spollineages.csv"; # input file for spollineages
my $phyloviz_entries = "phyloviz_profiles.tab"; # input file for phyloviz
my $grapetree_entries = "grapetree.tab"; # input file for grapetree
my %hSPOL = (); # hash map for all Spol sequences
my %hSPOLoctal = ();
my %hMIRU = (); # hash for 24-loci MIRU-VNTRs
my %h12miru = (); # hash for 12-loci MIRU-VNTRs
my %h15miru = (); # hash for 15-loci MIRU-VNTRs
my @tabFiles = ();
#my @tabRuns = ();
my %hResistance = (); 
my %hLineageMyk = (); # lineage from mykrobe
my %galruSpol = ();
my %spollineages = ();
my $isList = 0;
my $isFasta = 0;
my $isFastq = 0;
my $tbProf = "TBProfiler_result";
my $resultTBprof = "";
my %hLineageTBP = ();
my %hDrugTBP = ();
my %hDrugListTBP = ();
my %hFLC = ();  # results for fast-lineage-caller program
my %hSB =(); 
my $sbFile = "SBnumbers_Mbovis_240221.tab";
my $useMykrobe = 0;
my $useTBP = 0;
my $doMiru = 0;
my $useGrape = 0;
my $useRoary = 0;
my $fasttree = 0;
my $spades = 0;
my $sbPath = "";
my $outdir = "simpiTBResults";
my $removeS = 0;
my $pathToBin = ""; # "./bin/" replaced by ""
my $useResFinder = 0;
my %hReads = ();
my %hRF = ();
my $genCode = 11;
my %hGrape = ();
my %cladeInfo = ();

my $version = "1.0.1";

my $start = time();

print "##################################################################\n";
print "# --> Welcome to $0 (version $version)!\n";
print "##################################################################\n";

print "\n";

my $spadesProg = isProgInstalled("spades");
if($spadesProg){
  print "spades is...............OK \n";
}
my $prodigalProg = isProgInstalled("prodigal");
if($prodigalProg){
  print "prodigal is...............OK \n";
}
my $tbprofilerProg = isProgInstalled("tb-profiler");
if($tbprofilerProg){
  print "tb-profiler is...............OK \n";
}
my $roaryProg = isProgInstalled("roary");
if($roaryProg){
  print "roary is...............OK \n";
}
my $fasttreeProg = isProgInstalled("fasttree");
if($fasttreeProg){
  print "fasttree is...............OK \n";
}
my $grapetreeProg = isProgInstalled("grapetree");
if($grapetreeProg){
  print "grapetree is...............OK \n";
}
my $flcProg = isProgInstalled("fast-lineage-caller");
if($flcProg){
  print "fast-lineage-caller is...............OK \n";
}

print "\n";

if (@ARGV<1) {
  help_user_simple($0);
  exit 1;
}

if ($ARGV[0] eq "--help" || $ARGV[0] eq "-h") {
  help_user_advance($0);
  exit 0;
}

if ($ARGV[0] eq "--version" || $ARGV[0] eq "-v") {
  program_version($0);
  exit 0;	
}

#FASTA/Q files
if(@ARGV){
  for my $arg (@ARGV){
    if ($arg =~ m/.fasta/ or $arg =~ m/.fna/ ){  # or $arg =~ m/.fasta.gz/ or $arg =~ m/.fna.gz/
        print "FASTA file: $arg\n";
        #$seqfasta = $arg;

        push (@tabFiles,$arg);
        $isFasta = 1;
    }
    elsif ($arg =~ m/.fastq/ or $arg =~ m/.fq/ or $arg =~ m/.fq.gz/ or $arg =~ m/.fastq.gz/){
	print "FASTQ file: $arg\n";
        #$seqfastq = $arg;
        my $readName = $arg;
        if ($arg =~ m/_1.fastq.gz/ or $arg =~ m/_R1.fastq.gz/ or $arg =~ m/_2.fastq.gz/ or $arg =~ m/_R2.fastq.gz/ ){
          #my $newstring = $oldstring;
          $readName =~ s/_1.fastq.gz//g;
          $readName =~ s/_R1.fastq.gz//g;
          $readName =~ s/_2.fastq.gz//g;
          $readName =~ s/_R2.fastq.gz//g;
          #print "READ NAME is : $readName\n";
          $hReads{$readName} = $readName;
          @tabFiles = keys %hReads;  
          #push (@tabFiles,$readName);
          $isFastq = 0;
          $isList = 1;
        }
        else{
          push (@tabFiles,$arg);
          $isFastq = 1;
        }
        
    }
    elsif ($arg =~ m/.txt/){
        if (-e $arg){
            open my $handle, '<', $arg;
            chomp(@tabFiles = <$handle>);
            close $handle;
            $isList = 1;
        }
    }
    elsif ($arg =~ m/,/){
        @tabFiles = split (/,/, $arg) ;
            
        $isList = 1;
    }
  }
}  

##requirements
for (my $i = 0; $i <= $#ARGV; $i++) {
    if ($ARGV[$i]=~/--mykrobe/i or $ARGV[$i]=~/-mk/i) {
        $useMykrobe = 1;
    }
    elsif ($ARGV[$i]=~/--tbprofiler/i or $ARGV[$i]=~/-tp/i) {
        $useTBP = 1;
    }
    elsif ($ARGV[$i]=~/--miru/i or $ARGV[$i]=~/-mi/i) {
        $doMiru = 1;
    }
    elsif ($ARGV[$i]=~/--grape/i or $ARGV[$i]=~/-gr/i) {
        $useGrape = 1;
    }
    elsif ($ARGV[$i]=~/--roary/i or $ARGV[$i]=~/-ro/i) {
        $useRoary = 1;
    }
    elsif ($ARGV[$i]=~/--fasttree/i or $ARGV[$i]=~/-ft/i) {
        $fasttree = 1;
    }
    elsif ($ARGV[$i]=~/--spades/i or $ARGV[$i]=~/-sp/i) {
        $spades = 1;
    }
    elsif ($ARGV[$i]=~/--out/i or $ARGV[$i]=~/-o/i) {
        $outdir = $ARGV[$i+1];
    }
    elsif ($ARGV[$i]=~/--sb_path/i or $ARGV[$i]=~/--sbpath/i or $ARGV[$i]=~/-sb/i) {
        $sbPath = $ARGV[$i+1];
    }
    elsif ($ARGV[$i]=~/--remove/i or $ARGV[$i]=~/-rm/i) {
        $removeS = 1;
    }
    elsif ($ARGV[$i]=~/--resfinder/i or $ARGV[$i]=~/-rf/i) {
        $useResFinder = 1;
    }
    elsif ($ARGV[$i]=~/--pathToBin/i or $ARGV[$i]=~/--path/i or $ARGV[$i]=~/-p/i) {
        $pathToBin = $ARGV[$i+1];
    }
}

if(-d "GFF/") { system("rm -rf GFF/");  }
if(-d $outdir) { system("rm -rf $outdir");  }
mkdir($outdir);

    if(-e "${pathToBin}/SpoTyping-2.1/SpoTyping-v2.1-commandLine/SpoTyping.py"){
      print "spotyping is...............OK \n";
    }
    if(-e "${pathToBin}/MIRUReader/MIRUReader.py"){
      print "mirureader is...............OK \n";
    }
    if(-e "${pathToBin}/resfinder/run_resfinder.py"){
      print "resfinder is...............OK \n";
    }
    if(-e "${pathToBin}/SpolLineages/spollineages.jar"){
      print "spollineages is...............OK \n";
    }
    
foreach my $seq (@tabFiles){

    my $seqFolder = "spades".$seq;
    #$newstring = $oldstring;
    if ($seqFolder =~  m/\./) {
        #my @tabSeq = split (/\./, $seqFolder) ;
        #$seqFolder =~ $tabSeq[0];
        $seqFolder =~ s/.fastq.gz//g;
        $seqFolder =~ s/.fastq//g;
       print "SEQFOLDER IS : $seqFolder \n";
    }

    #check paired reads and consider as list


    #SpoTyping, MIRUReader, Galru, and MyKrobe and TBprofiler
    my $spotypingCmd = "";
    my $tbProfilerCmd = "tb-profiler profile --txt "; # DC removed $pathToBin.
    my $fastaTbProfilerCmd = "tb-profiler fasta_profile --txt ";
    #print "TBProfiler command: $tbProfilerCmd\n";
    my $resfinderCmd = "python3 ${pathToBin}/resfinder/run_resfinder.py -o RF -s \"Mycobacterium tuberculosis\" -l 0.6 -t 0.8 --acquired --point -db_res ${pathToBin}/resfinder/db_resfinder -db_point ${pathToBin}/resfinder/db_pointfinder -ifa";
    #RF/pheno_table.txt

    if(-e $seq and $isFasta){
        $spotypingCmd = "python2 ".$pathToBin."SpoTyping-2.1/SpoTyping-v2.1-commandLine/SpoTyping.py --noQuery -o spotyping --seq $seq ";
        #open GALRU, "galru $seq -t 16 |";
        if ($doMiru) { 
          print "MIRU analysis from FASTA input files... \n";
          open MIRU, "python3 ${pathToBin}MIRUReader/MIRUReader.py -r $seq -p mirus |"; 
          readMiru($seq);
        }
        else { $hMIRU{$seq} = "NA"; }
        if ($useMykrobe) {
          #open KROBE, 'sudo singularity exec -B $PWD docker://"quay.io/biocontainers/mykrobe:0.9.0--py37h13b99d1_2" mykrobe predict my_sample tb -f --format csv --seq '. "$seq |";
          open KROBE, "mykrobe predict -s my_sample --species tb -f --format csv --seq $seq |"; #-s SAMPLE --species tb 
	}
        else {  
          $hResistance{$seq} = "NA";
          $hLineageMyk{$seq} = "NA";
        }
        if($useResFinder) {
          $resfinderCmd .= " ".$seq;
          system ($resfinderCmd);
          resfinder_out($seq);
          if(-d "RF/") {
            system ("rm -rf RF/ "); #remove RF folder
          }
          #elsif(-d "RF/"){
          #  system ("mv -f RF/ $outdir/");
          #}
          #else{
          #}
        }
        if($useRoary){
          my $nameSeq = $seq;
          $nameSeq =~  s/\//_/ig;
          $nameSeq =~  s/.fasta//ig;
          $nameSeq =~  s/.fna//ig;
	  
	  #my $prokkaCmd = "prokka --force --outdir prokka --prefix result $seq";
          my $prokkaCmd = "prodigal -i $seq -c -m -g $genCode -f gff -o my_gff ";
	  system($prokkaCmd);
          my $gff = "my_gff";
	  system(" echo '##FASTA' >> $gff && cat $gff $seq > ${nameSeq}.gff ");
          if(! -d "GFF/") { mkdir("GFF/");  }
         
          if(-e "${nameSeq}.gff" and -d "GFF/") { system("mv ${nameSeq}.gff GFF/ "); }
          
        }

        $fastaTbProfilerCmd = $fastaTbProfilerCmd." -d ".$tbProf." ".$seq." ".$seq;
    }
    elsif(-e $seq and $isFastq){
        
        $spotypingCmd = "python2 ".$pathToBin."SpoTyping-2.1/SpoTyping-v2.1-commandLine/SpoTyping.py $seq --noQuery -o spotyping ";
        
        #open GALRU, "galru $seq -t 16 |";
        if ($doMiru) {
          if($spades) { 
            my $nameSeq = $seq;
            #$nameSeq =~  s/\//_/ig; #commented to replace by $seqFolder
            if(-d $seqFolder) { system("rm -rf spades$nameSeq");  } #spades$nameSeq replaced by $seqFolder
            my $spadesCmd = "spades --careful -s $seq -o $seqFolder "; #--only-assembler replaced by --careful
            my $pathContig = $seqFolder."/scaffolds.fasta";
            system($spadesCmd);

            if(-e $pathContig) { 
              open MIRU, "python3 ".$pathToBin."MIRUReader/MIRUReader.py -r $pathContig -p mirus |";  
              readMiru($seq);

              if($useResFinder) {
                $resfinderCmd .= " ".$pathContig;
                system ($resfinderCmd);
                resfinder_out($seq);

                if(-d "RF/") {
                  system ("rm -rf RF/ "); #remove RF folder
                }
                #elsif(-d "RF/"){
                #  system ("mv -f RF/ $outdir/");
                #}
                #else{
                #}
              }

            }
            if(-d $seqFolder)  { 
              if ($outdir =~  m/\//) { $outdir =~  s/\///i; }
              system ("cp -f $pathContig $outdir/scaffolds_$seq.fasta"); #removed scaffolds_$seq.fasta? 
              if($removeS) {  
                system ("rm -rf $seqFolder"); 
              }
              else { system ("mv -f $seqFolder $outdir/");  }
            }

          }
          else{
            open MIRU, "python3 ".$pathToBin."MIRUReader/MIRUReader.py -r $seq -p mirus --nofasta |"; 
            readMiru($seq);
          }
        }
        else { 
          $hMIRU{$seq} = "NA"; 
          $h15miru{$seq} = "NA";
          $h12miru{$seq} = "NA";
        }

        if ($useMykrobe) {
          #open KROBE, 'sudo singularity exec -B $PWD docker://"quay.io/biocontainers/mykrobe:0.9.0--py37h13b99d1_2" mykrobe predict my_sample tb --format csv --seq '. "$seq |";
          open KROBE, "mykrobe predict -s my_sample --species tb -f --format csv --seq $seq |";
        }
        else {  
          $hResistance{$seq} = "NA";
          $hLineageMyk{$seq} = "NA";
        }

        $tbProfilerCmd = $tbProfilerCmd." -1 ".$seq." -p ".$seq." -t 16 -d ".$tbProf;
    }
    elsif($isList){
      my $read1 = $seq."_1.fastq.gz";
      my $read2 = $seq."_2.fastq.gz";
      if(-e $seq."_R1.fastq.gz" and $seq."_R2.fastq.gz" ) {
        $read1 = $seq."_R1.fastq.gz";
        $read2 = $seq."_R2.fastq.gz";
      }

      if(-e $read1 and -e $read2){
        $spotypingCmd = "python2 ".$pathToBin."SpoTyping-2.1/SpoTyping-v2.1-commandLine/SpoTyping.py $read1 $read2 --noQuery -o spotyping ";
        #open GALRU, "galru $read1 -t 16 |";
        if ($doMiru) { 
          if($spades) {

            if(-d $seqFolder) { system("rm -rf $seqFolder");  } 
            my $spadesCmd = "spades --careful -1 $read1 -2 $read2 -o $seqFolder "; # --disable-rr  --isolate? default (last run was 2772 seconds) --only-assembler (2526 seconds) good mirus
            my $pathContig = $seqFolder."/scaffolds.fasta";
            system($spadesCmd);

            if(-e $pathContig) { 
              open MIRU, "python3 ".$pathToBin."MIRUReader/MIRUReader.py -r $pathContig -p mirus |";  
              readMiru($seq);

              if($useResFinder) {
                $resfinderCmd .= " ".$pathContig;
                system ($resfinderCmd);
                resfinder_out($seq);

                if(-d "RF/") {
                  system ("rm -rf RF/ "); #remove RF folder
                }
                #elsif(-d "RF/"){
                #  system ("mv -f RF/ $outdir/");
                #}
                #else{
                #}

              }

            }

            if(-d $seqFolder)  { 
              if ($outdir =~  m/\//) { $outdir =~  s/\///i; }
              system ("cp -f $pathContig $outdir/scaffolds_$seqFolder.fasta"); 
              if($removeS) {  
                system ("rm -rf $seqFolder"); 
              }
              else { system ("mv -f $seqFolder $outdir/");  }
            }
          }
          else {
            open MIRU, "python3 ".$pathToBin."MIRUReader/MIRUReader.py -r $read1 -p mirus --nofasta |"; 
            readMiru($seq);
          }
        }
        else { 
          $hMIRU{$seq} = "NA"; 
          $h15miru{$seq} = "NA";
          $h12miru{$seq} = "NA";
        }
        
        if ($useMykrobe) {
          #open KROBE, 'sudo singularity exec -B $PWD docker://"quay.io/biocontainers/mykrobe:0.9.0--py37h13b99d1_2" mykrobe predict my_sample tb --format csv --seq '. "$read1 $read2 |";
          open KROBE, "mykrobe predict -s my_sample --species tb -f --format csv --seq $read1 $read2 |";
        }
        else {  
          $hResistance{$seq} = "NA";
          $hLineageMyk{$seq} = "NA";
        }

        $tbProfilerCmd = $tbProfilerCmd." -t 16 -d ".$tbProf." -1 ".$read1." -2 ".$read2." -p ".$seq;
      }
      elsif(-e $read1 and ! -e $read2){
        $spotypingCmd = "python2 ".$pathToBin."SpoTyping-2.1/SpoTyping-v2.1-commandLine/SpoTyping.py $read1 --noQuery -o spotyping ";
        #open GALRU, "galru $read1 -t 16 |";
        if ($doMiru) {  
          if($spades) { 
            #my $seqFolder = "spades".$seq;
            #$newstring = $oldstring;

            my $spadesCmd = "spades --careful -1 $read1 -o $seqFolder "; #--only-assembler replaced by --careful
            my $pathContig = "${seqFolder}/scaffolds.fasta";
            system($spadesCmd);

            if(-e $pathContig) { 
              open MIRU, "python3 ".$pathToBin."MIRUReader/MIRUReader.py -r $pathContig -p mirus |";  
              readMiru($seq);

               if($useResFinder) {
                $resfinderCmd .= " ".$pathContig;
                system ($resfinderCmd);
                resfinder_out($seq);

                if(-d "RF/") {
                  system ("rm -rf RF/ "); #remove RF folder
                }
                #elsif(-d "RF/"){
                #  system ("mv -f RF/ $outdir/");
                #}
                #else{
                #}

              }
            }
            
            if(-d $seqFolder)  { 
              if ($outdir =~  m/\//) { $outdir =~  s/\///i; }
              system ("cp -f $pathContig $outdir/scaffolds_$seq.fasta"); 
              system ("mv -f spades$seq $outdir/"); 
            }
          }
          else {
            open MIRU, "python3 ".$pathToBin."MIRUReader/MIRUReader.py -r $read1 -p mirus --nofasta |"; 
            readMiru($seq);
          }
        }
        else { 
          $hMIRU{$seq} = "NA"; 
          $h15miru{$seq} = "NA";
          $h12miru{$seq} = "NA";
        }
        if ($useMykrobe) {
          #open KROBE, 'sudo singularity exec -B $PWD docker://"quay.io/biocontainers/mykrobe:0.9.0--py37h13b99d1_2" mykrobe predict my_sample tb --format csv --seq '. "$read1 |";
          open KROBE, "mykrobe predict -s my_sample --species tb --format csv --seq $read1 |";
        }
        else {  
          $hResistance{$seq} = "NA";
          $hLineageMyk{$seq} = "NA";
        }
        $tbProfilerCmd = $tbProfilerCmd." -t 16 -d ".$tbProf." -1 ".$read1." -p ".$seq;
      }
  
    }
    else{
        print "Input file does not exist or is not correct....please verify! \n";
        exit(1);
    }

    system ($spotypingCmd);

    #read SpoTyping/Galru results
    open (SPO, "<spotyping") or die "open : $!";
    while (<SPO>) {
        chomp();
	print "SPO file: $_\n";
	my @tab = split (/\t/, $_) ;
        $hSPOL{$seq} = $tab[1]."\t".$tab[2]; # $tab[0] replaced by $seq
        $hSPOLoctal{$seq} = $tab[2];
	#$hSPOL{$tab[0]} = $tab[1]."\t".$tab[2]; #$seq
    }
    close (SPO) or die "close file error : $!";
    system ("rm -f spotyping "); 
    system ("rm -f spotyping.log "); 
    system ("rm -f spotyping.SpoTyping.tmp.* ");
    system ("rm -rf *.fasta.* ");  # remove files generated when using MIRUReader
    system ("rm -rf *.fna.* ");
    #system ("rm -rf RF/ "); #remove RF folder
        
    #while (<GALRU>) {
     #   chomp();
     #   $galruSpol{$seq} = "$_";
    #}
    #close (GALRU) or die "close file error : $!";

    # gather SB numbers in hashTable taking into account SBnumbers_Mbovis_240221.tab 
    if(-e $sbPath){
      open (SB, "<$sbPath") or die "open : $!";
    }
    elsif (-e $sbFile) {
      open (SB, "<$sbFile") or die "open : $!";
    }

    if(-e $sbPath or -e  $sbFile) {
    while (<SB>) {
        chomp();
        if ($_ =~  m/^SB/) {
	  my @tabSB = split (/\t/, $_) ;
          if ($tabSB[2] ne '') { $hSB{$tabSB[2]} = $tabSB[0]; }  # $tabSB[2] replaced by $seq , and removed $tabSB[2] ne ''
	  #$hSPOL{$tab[0]} = $tab[1]."\t".$tab[2]; #$seq
        }
    }
    close (SB) or die "close file error : $!";
    }
 
    if ($useMykrobe) {
    while (<KROBE>) {
        chomp();
        #print "WE ARE IN KROBE RESULTS:\n";
        print $_."\n";
        if($_ =~  m/^\"my_sample\"/){
          #print "MY MAIN LINE IS: $_\n";
          my @tabK = split (/,/, $_) ;
          print "drug: $tabK[1]; susceptibility: $tabK[2]; lineage: $tabK[12]; lineage_per_covg: $tabK[15]; lineage_depth: $tabK[18]. \n\n";
          $tabK[2] =~  s/\"//g;
          $tabK[12] =~  s/\"//g;
          if($tabK[2] eq "R"){ $hResistance{$seq} .= "$tabK[1] ($tabK[2]) "; }
          #$hResistance{$seq} =~  s/\"//g;
          $hLineageMyk{$seq} = "$tabK[12] "; # replacing .= by =
          #$hLineageMyk{$seq} =~  s/\"//g;
        }
    }
    close (KROBE) or die "close file error : $!";
    }

    #TBprofiler results
    if($isFastq or $isList){
        system($tbProfilerCmd);

        #read results
        $resultTBprof = $tbProf."/results/".$seq.".results.txt";

	open (TBP, "<$resultTBprof") or die "open : $!";
        
        my $lineageTBP = "";
        my $drugTBP = "";
        my $drugListTBP ="";
        my @tabTBP = ();
        
        while (<TBP>) {
          chomp();
          if ($_ =~  m/^Strain/) {
            $_ =~ s/Strain: //;
            $lineageTBP = $_;
            
          }
          elsif ($_ =~  m/^Drug-resistance/) {
            $_ =~ s/Drug-resistance: //;
            $drugTBP = $_;
            
          }
          elsif ($_ =~ m/^Rifampicin/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") { $drugListTBP .= "Rifampicin (R) "; }
          }
          elsif ($_ =~ m/^Ethambutol/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") { $drugListTBP .= "Ethambutol (R) "; }
          }
          elsif ($_ =~ m/^Pyrazinamide/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Pyrazinamide (R) "; }
            print "IN TBProlier File: $drugListTBP\n";
          }
          elsif ($_ =~ m/^Streptomycin/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Streptomycin (R) "; }
          }
          elsif ($_ =~ m/^Fluoroquinolones/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Fluoroquinolones (R) "; }
          }
          elsif ($_ =~  m/^Amikacin/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Amikacin (R) "; }
          }
          elsif ($_ =~  m/^Capreomycin/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Capreomycin (R) "; }
          }
          elsif ($_ =~  m/^Kanamycin/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Kanamycin (R) "; }
          }
          elsif ($_ =~  m/^Cycloserine/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Cycloserine (R) "; }
          }
          elsif ($_ =~  m/^Ethionamide/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Ethionamide (R) "; }
          }
          elsif ($_ =~  m/^Clofazimine/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Clofazimine (R) "; }
          }
          elsif ($_ =~  m/^Para-aminosalicylic acid/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Para-aminosalicylic acid (R) "; }
          }
          elsif ($_ =~  m/^Delamanid/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Delamanid (R) "; }
          }
          elsif ($_ =~  m/^Bedaquiline/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Bedaquiline (R) "; }
          }
          elsif ($_ =~  m/^Linezolid/) {
            #$_ =~ s/Drug-resistance: //;  
            @tabTBP = split (/\t/, $_) ;
            if(defined($tabTBP[1]) and $tabTBP[1] eq "R") {$drugListTBP .= "Linezolid (R) "; }
          }
          #else{
           # $drugListTBP = "ND";
          #}
        }

        close (TBP) or die "close file error : $!";
        $hLineageTBP{$seq} = $lineageTBP ;
        $hDrugTBP{$seq} = $drugTBP;
        $hDrugListTBP{$seq} = $drugListTBP;
	
        #Fast-Lineage-Caller   (please note that this program is used when TBProf is used because it needs the .vcf file generated by TBP)
        # @ gzip vcf file from TBP directory (gzip -d 
        my $vcfTBprof = $tbProf."/vcf/".$seq.".targets.vcf.gz";   # .targets.csq.vcf.gz
        system ("gzip -d -f $vcfTBprof");
        my $newVCF = $tbProf."/vcf/".$seq.".targets.vcf"; 

        open FLC, "fast-lineage-caller --noheader $newVCF |";
        # read FLC results
        while (<FLC>) {
          chomp();
          my @flc = split (/\t/, $_) ;
          my $first = shift(@flc); 
          $hFLC{$seq} = join( "\t", @flc ); 
          #my $stringFLC = join( "\t", @flc ); 
        }
        close (FLC) or die "close file error : $!";
	
	if($removeS){
	  system("rm -rf $tbProf ");
	}
    }
     
}


#########################
#Summary files
open (RECAP,'>', $recap_total_seq) or die "could not open $!"; # removed Spoligotype (Galru)\t
print RECAP "#ID\tSpoligotype (SpoTyping)\tSpoligo (octal)\t24-loci MIRU-VNTR\t15-loci MIRU-VNTR\t12-loci MIRU-VNTR\tLineage (SpolLineages)\tSpoligotype International Type (SIT)\tSBnumber (Mbovis)\tResistance (Mykrobe)\tLineage (Mykrobe)\tResistance (TBProfiler)\tLineage (TBProfiler)\tResFinder\tColl2014\tFreschi2020\tLipworth2019\tShitikov2017\tStucki2016\n";

#SpolLineages and GrapeTree files
open (LINE,'>', $spollineages_entries) or die "could not open $!";

if($useGrape){
open (PHYLO,'>', $phyloviz_entries) or die "could not open $!";
if($doMiru){
  print PHYLO "#Strain\tGene_1\tGene_2\tGene_3\tGene_4\tGene_5\tGene_6\tGene_7\tGene_8\tGene_9\tGene_10\tGene_11\tGene_12\tGene_13\tGene_14\tGene_15\tGene_16\tGene_17\tGene_18\tGene_19\tGene_20\tGene_21\tGene_22\tGene_23\tGene_24\tGene_25\tGene_26\tGene_27\tGene_28\tGene_29\tGene_30\tGene_31\tGene_32\tGene_33\tGene_34\tGene_35\tGene_36\tGene_37\tGene_38\tGene_39\tGene_40\tGene_41\tGene_42\tGene_43\tGene_44\tGene_45\tGene_46\tGene_47\tGene_48\tGene_49\tGene_50\tGene_51\tGene_52\tGene_53\tGene_54\tGene_55\tGene_56\tGene_57\tGene_58\tGene_59\tGene_60\tGene_61\tGene_62\tGene_63\tGene_64\tGene_65\tGene_66\tGene_67\n";
}
else{
  print PHYLO "#Strain\tGene_1\tGene_2\tGene_3\tGene_4\tGene_5\tGene_6\tGene_7\tGene_8\tGene_9\tGene_10\tGene_11\tGene_12\tGene_13\tGene_14\tGene_15\tGene_16\tGene_17\tGene_18\tGene_19\tGene_20\tGene_21\tGene_22\tGene_23\tGene_24\tGene_25\tGene_26\tGene_27\tGene_28\tGene_29\tGene_30\tGene_31\tGene_32\tGene_33\tGene_34\tGene_35\tGene_36\tGene_37\tGene_38\tGene_39\tGene_40\tGene_41\tGene_42\tGene_43\n";
}

#my $count = 0;
}

#foreach my $key (keys %hSPOL)
foreach my $key (@tabFiles)
{
  my @tabSpo = split (/\t/, $hSPOL{$key}) ;
  print LINE "$key;$tabSpo[1];$hMIRU{$key}\n";

  if($useGrape){
    my @spo43 = split (//, $tabSpo[0]) ;
    my $spolString = join( "\t", @spo43 ); 

    my @miru24 = split (//, $hMIRU{$key}) ;
    my $miruString = join( "\t", @miru24 );
    
    $miruString =~  s/A/10/ig;
    $miruString =~  s/B/11/ig;
    $miruString =~  s/C/12/ig;
    $miruString =~  s/D/13/ig;
    $miruString =~  s/E/14/ig;
    $miruString =~  s/F/15/ig;
    #$miruString =~  s/-/0/ig;
    my @tabKey = split (/\./, $key) ;
    #my $id = $count++; # $tabKey[0]
	if($doMiru){
          print PHYLO "$tabKey[0]\t$spolString\t$miruString\n";
	  my $spoligomiru = $tabSpo[0].$hMIRU{$key};
	  $hGrape{$spoligomiru} .= $tabKey[0]."|";
	}
	else{
	  print PHYLO "$tabKey[0]\t$spolString\n";
	  my $spoligo = $tabSpo[0];
	  $hGrape{$spoligo} .= $tabKey[0]."|";
	}
  }
  
}

if($useGrape){ close (PHYLO) or die "close file error : $!"; }
close (LINE) or die "close file error : $!";  


#SpolLineages
my $spollineages_Cmd = "java -jar ".$pathToBin."SpolLineages/spollineages.jar -i $spollineages_entries -o outputSL.csv -cur -c -sit -uR -S ";
system($spollineages_Cmd);

      open (RESLINE, "<outputSL.csv") or die "open : $!";
      while (<RESLINE>) {
        chomp();
        if ($_ !~  m/^StrainID/) {
          my @tabRes = split (/;/, $_) ;
          $spollineages{$tabRes[0]} = "$tabRes[4]"."/"."$tabRes[5]\t$tabRes[12]";
	  my $currentSpol = $tabRes[1];
	  $currentSpol =~  s/n/1/ig;
	  $currentSpol =~  s/o/0/ig;
	  $cladeInfo{$currentSpol} = "$tabRes[4]"."/"."$tabRes[5]_SIT$tabRes[12]";
        }
      }
      close (RESLINE) or die "close file error : $!"; 

if(-e $spollineages_entries and -e "outputSL.csv") {   
  system("mv $spollineages_entries outputSL.csv $outdir");
}

#removed "and -e $grapetree_entries" in condition
if($useGrape){
#prepare entry-file
open (GRAPE,'>', $grapetree_entries) or die "could not open $!";
if($doMiru){
  print GRAPE "#Strain\tGene_1\tGene_2\tGene_3\tGene_4\tGene_5\tGene_6\tGene_7\tGene_8\tGene_9\tGene_10\tGene_11\tGene_12\tGene_13\tGene_14\tGene_15\tGene_16\tGene_17\tGene_18\tGene_19\tGene_20\tGene_21\tGene_22\tGene_23\tGene_24\tGene_25\tGene_26\tGene_27\tGene_28\tGene_29\tGene_30\tGene_31\tGene_32\tGene_33\tGene_34\tGene_35\tGene_36\tGene_37\tGene_38\tGene_39\tGene_40\tGene_41\tGene_42\tGene_43\tGene_44\tGene_45\tGene_46\tGene_47\tGene_48\tGene_49\tGene_50\tGene_51\tGene_52\tGene_53\tGene_54\tGene_55\tGene_56\tGene_57\tGene_58\tGene_59\tGene_60\tGene_61\tGene_62\tGene_63\tGene_64\tGene_65\tGene_66\tGene_67\n";
}
else{
  print GRAPE "#Strain\tGene_1\tGene_2\tGene_3\tGene_4\tGene_5\tGene_6\tGene_7\tGene_8\tGene_9\tGene_10\tGene_11\tGene_12\tGene_13\tGene_14\tGene_15\tGene_16\tGene_17\tGene_18\tGene_19\tGene_20\tGene_21\tGene_22\tGene_23\tGene_24\tGene_25\tGene_26\tGene_27\tGene_28\tGene_29\tGene_30\tGene_31\tGene_32\tGene_33\tGene_34\tGene_35\tGene_36\tGene_37\tGene_38\tGene_39\tGene_40\tGene_41\tGene_42\tGene_43\n";
}

foreach my $gkey (keys %hGrape){
 my $firstSpol_char = "";

 my @string_2 = split (//, $gkey) ;
 my $profiles_2 = join( "\t", @string_2 ); 
  
  if($doMiru){
    $firstSpol_char = substr($gkey, 0, 43); # 'first 43 chars correspond to spoligo'
    $gkey =~  s/A/10/ig;
    $gkey =~  s/B/11/ig;
    $gkey =~  s/C/12/ig;
    $gkey =~  s/D/13/ig;
    $gkey =~  s/E/14/ig;
    $gkey =~  s/F/15/ig; 
    
    print GRAPE "$hGrape{$gkey}".$cladeInfo{$firstSpol_char}."\t$profiles_2\n";
  }
  else{
    print GRAPE "$hGrape{$gkey}".$cladeInfo{$gkey}."\t$profiles_2\n";
  }
  
}
close (GRAPE) or die "close file error : $!"; 

#GrapeTree
if(-e $grapetree_entries) { 
my $grape_Cmd = "grapetree -p $grapetree_entries -m NJ > grapetreeNJ.nwk";
my $tabSize = @tabFiles;
if($tabSize >= 4){   # the NJ tree file will be generated only if $tabSize is >= 4
  system ($grape_Cmd);
}
system("mv $grapetree_entries $outdir");
}

if(-e "grapetreeNJ.nwk") { 
  system("mv grapetreeNJ.nwk $outdir"); 
}
}


#Final result   hSPOLoctal
foreach my $key2 (@tabFiles)
{

  if ($hDrugListTBP{$key2} eq '') { $hDrugListTBP{$key2}= "NA"; }
  if ($hResistance{$key2} eq '') { $hResistance{$key2}= "NA"; }
  if ($hLineageMyk{$key2} eq '') { $hLineageMyk{$key2} = "NA"; }
  if ($hDrugTBP{$key2} eq '') { $hDrugTBP{$key2} = "NA"; }
  if ($hLineageTBP{$key2} eq '') { $hLineageTBP{$key2} = "NA"; }
  if ($hSB{$hSPOLoctal{$key2}} eq '' ) { $hSB{$hSPOLoctal{$key2}} = "NA"; } 
  if ($hRF{$key2} eq '' ) { $hRF{$key2} = "NA"; }

  # else { $hSB{$hSPOLoctal{$key2}} = "NA"; }
  # print "  $key2 spoligotype is: $hSPOL{$key2} ; GALRU: $galruSpol{$key2} ;  24-loci MIRU value is: $hMIRU{$key2}\n"; # removed $galruSpol{$key2}\t
  print RECAP "$key2\t$hSPOL{$key2}\t$hMIRU{$key2}\t$h15miru{$key2}\t$h12miru{$key2}\t$spollineages{$key2}\t$hSB{$hSPOLoctal{$key2}}\t$hResistance{$key2}\t$hLineageMyk{$key2}\t$hDrugTBP{$key2}: $hDrugListTBP{$key2}\t$hLineageTBP{$key2}\t$hRF{$key2}\t$hFLC{$key2}\n";

}

close (RECAP) or die "close file error : $!";  


if($useRoary and -d "GFF/"){
  if(-d "roary") { system("rm -rf roary"); }
  my $roaryCmd = "roary -e --mafft -p 16 -g 100000 -f roary -v GFF/*.gff -cd 95 -r ";
  system($roaryCmd);

  #dnaDist & fastme / fasttree
  my $coreAln = "roary/core_gene_alignment.aln";
  if (-e $coreAln) { system ("cp $coreAln $outdir"); }
  
  if($fasttree){
    my $fasttreeCmd = "fasttree -nt -gtr < $coreAln > fasttree.nwk";
    system($fasttreeCmd);
    if (-e "fasttree.nwk") { system ("mv fasttree.nwk $outdir"); }
  }
  else {
    #my $dnaDistCmd = "/opt/dnaDist/build/dnaDist -b 100 $coreAln > cg_dnaDist.mat";
    #my $fastmeCmd = "/media/results_datacalcul/David/SLAMTB/fastme-2.1.5/binaries/fastme-2.1.5-linux64 -i cg_dnaDist.mat -o fastme_tree.nwk";

    #system($dnaDistCmd);
    #system($fastmeCmd);
    #if (-e "cg_dnaDist.mat" and -e "fastme_tree.nwk") { system ("mv cg_dnaDist.mat fastme_tree.nwk $outdir"); }
  }
}
#remove or move files generated and used by roary
if($removeS){
  system("rm -rf roary");
  system("rm -rf GFF/");
}
else{
  system("mv roary $outdir");
  system("mv GFF/ $outdir");
}
# Remove TBP repository
#if(-d $resultTBprof){
#  system("rm -rf $resultTBprof");
#  rmdir $resultTBprof;
#}

# Move files and directories to outdir
if (-d $tbProf) { system ("mv $tbProf $outdir"); }
if (-e $recap_total_seq) { system ("mv $recap_total_seq $outdir"); }

my $end = time();

my $total = $end - $start;
my $min = $total / 60;
my $hrs = $min / 60;

print "\n\n";
print "***** Thank you for using $0! \n";
print "***** Total time: $total seconds OR $min minutes OR $hrs hours *****\n";


### functions
#readMiru
sub readMiru {
# read MIRU results
  my $seq = shift @_;
    if ($doMiru) { 
    while (<MIRU>) {
        chomp();
        if ($_ =~  m/^mirus/) {
          print "24-loci $_\n";
          $_ =~  s/s//g;

          $_ =~  s/10/A/ig;
          $_ =~  s/11/B/ig;
          $_ =~  s/12/C/ig;
          $_ =~  s/13/D/ig;
          $_ =~  s/14/E/ig;
          $_ =~  s/15/F/ig;

          $_ =~  s/ND/-/ig;
          my @tm = split (/\t/, $_) ; # number (into brackets) indicates MIRU position according to SITVIT nomenclature 
	  
    #0154(1);0424(19);0577(15);0580(2);0802(12);0960(3);1644(4);1955(20);2059(5);2163b(16);2165(13);2347(21);2401(22);2461(14);2531(6);2687(7);2996(8);3007(9);3171(23);3192(10);3690;4052(17);4156(18);4348(11)

    $hMIRU{$seq} = $tm[1].$tm[4].$tm[6].$tm[7].$tm[9].$tm[15].$tm[16].$tm[17].$tm[18].$tm[20].$tm[24].$tm[5].$tm[11].$tm[14].$tm[3].$tm[10].$tm[22].$tm[23].$tm[2].$tm[8].$tm[12].$tm[13].$tm[19].$tm[21];

    $h12miru{$seq} = substr($hMIRU{$seq},0,12);  #=$tm[1].$tm[4].$tm[6].$tm[7].$tm[9].$tm[15].$tm[16].$tm[17].$tm[18].$tm[20].$tm[24].$tm[5];
    $h15miru{$seq} = $tm[4].$tm[6].$tm[7].$tm[17].$tm[20].$tm[5].$tm[11].$tm[3].$tm[10].$tm[22].$tm[23].$tm[2].$tm[8].$tm[13].$tm[21];

        }
    }  
    close (MIRU) or die "close file error : $!";
    }
}
# display global help document
sub help_user_simple {
	my $programme = shift @_;
	print STDERR  "Usage : perl $programme [options] *.fasta | *.fastq | accessions_list.txt | accession1,accession2...\n";
	print "Type perl $programme --version or perl $programme -v to get the current version\n";
	print "Type perl $programme --help or perl $programme -h to get full help\n";
}
#------------------------------------------------------------------------------
# display full help document
sub help_user_advance {
	print <<HEREDOC;
	
	Name: 
		$0
	
	Synopsis:
		A Perl pipeline to analyze WGS Mycobacterium tuberculosis complex (MTBC) data.
		
	Usage:
	  perl $0 [options] *.fasta | *.fastq | accessions_list.txt | accession1,accession2...
	  examples: 
       perl simpiTB.pl <options> *.fasta
       perl simpiTB.pl <options> *.fastq
       perl simpiTB.pl <options> *.fastq.gz
       perl simpiTB.pl <options> list.txt   
       (please note that the txt file must finish by extension ".txt", and it must contain one accession perl line such as the following example)
       #SRR13015794 
       #SRR13015795
       #accessions could also be separated by commas as follows: 
       perl simpiTB.pl SRR13015794,SRR13015795,SRR13015796
						 	
	General:
		--help or -h			displays this help 	
		--version or -v			displays the current version of the program
		
	Options ([XXX] represents the expected value):
		--out or -o [XXX]		allows to indicate the output repository containing generated files (default: $outdir)
		--path or -p [XXX]		allows to indicate the path to directory containing SpoTyping, MIRUReader, resfinder, and SpolLineages (default: "")
		--mykrobe or -mk		allows to use Mykrobe software to predict drug resistance profiles from sequencing reads
		--tbprofiler or -tp		allows to use TBProfiler software to predict drug resistance profiles from sequencing reads
		--resfinder or -rf		allows to use ResFinder tool to predict drug resistance from assembled genomes (FASTA)
		--miru or -mi			allows to use MIRUReader tool to predict 24-loci MIRU-VNTRs
		--grape or -gr			allows to use GrapeTree tool to infer a phylogeny from MIRU-VNTRs and Spoligotyping data
		--roary or -ro			allows to use Roary to perform a core gene alignment (pangenome analysis)
		--fasttree or -ft		allows to use FastTree to perform a phylogenetic analysis from alignment file
		--spades or -sp			allows to use SPAdes for genome assembly (when input files are FASTQ sequence reads)
		--sbpath or -sb [XXX]	allows to use SB numbers annotation collected from Mbovis.org database (an example file is provided in simpiTB repository)
		--remove or -rm			allows to remove supplementary output files/folders
		
HEREDOC
}
#------------------------------------------------------------------------------
# display program version 
sub program_version {
	my $programme = shift @_;
	print "\n $programme, version : $version\n";
	print "\n A Perl pipeline to analyze WGS Mycobacterium tuberculosis complex (MTBC) data. \n";
}
#------------------------------------------------------------------------------
sub resfinder_out {
  
  my $seq = shift @_;
  my @drugs = ('rifampicin', 'ethambutol', 'pyrazinamide', 'streptomycin', 'fluoroquinolone', 'amikacin', 'capreomycin', 'kanamycin', 'd-cycloserine', 'ethionamide', 'clofazimine', 'para-aminosalicylic acid', 'delamanid', 'bedaquiline', 'linezolid');
  my $res = "";
  my $out = "RF/pheno_table.txt";

        if(-e $out){
          for my $elem (@drugs){
            my $line = `grep '$elem.*Resistant' '$out'`;
            if($line){
              #my @tabdrug = split(/\t/,$line);
              $res .= "$elem (R), "; #"$tabdrug[0] (R), ";
            }
          }
          $hRF{$seq} = $res;
        }
        else{
          hRF{$seq} = "NA";
        }
}
#------------------------------------------------------------------------------

# Test if Prog is installed
sub isProgInstalled {
 
	my $program = shift;
 	my @PATH=split(":","$ENV{'PATH'}");     #Get users Path
	my $found = 0;
	foreach my $path (@PATH) 
	{
   	  if (-x "$path/$program") 
   	  {
   		  $found= 1;
   		  #last;
    	  }
	}
	
    return $found;
}

#------------------------------------------------------------------------------
