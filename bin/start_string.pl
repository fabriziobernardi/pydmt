#!/usr/bin/perl


 use strict;
 use Env;
 use Cwd;
 my $pa = getcwd;

 print "$pa\n";

 my $string = "pytdmt.py --model \/usr\/local\/pytdmtV2.1\/earth-models\/Aquila ";

 opendir(DIR,"."); my @dir=readdir(DIR); close(DIR);

 foreach $_ (@dir) {
   chomp($_);
   next if "$_" eq ".";
   next if "$_" eq "..";
   my @f = split(/\./,$_);
   if("$f[$#f]" eq "fseed") {
     $string .= "--fseed $_ ";
     last;
   }
 }

 #20130104_075006_37873_14722_0101

 my @d = split(/\//,$pa);

 my @l = split(/_/,$d[$#d]);
 my $Y = substr($l[0],0,4);
 my $M = substr($l[0],4,2);
 my $D = substr($l[0],6,2);
 my $h = substr($l[1],0,2);
 my $m = substr($l[1],2,2);
 my $s = substr($l[1],4,2);

 my $la = sprintf("%.3f",$l[2]/1000);
 my $lo = sprintf("%.3f",$l[3]/1000);
 my $pd = sprintf("%.1f",$l[4]/10);

 my $data = "--ori ".$Y."-".$M."-".$D."T".$h.":".$m.":".$s;
 my $epi  = "--epi \"$la $lo $pd\" ";

 my $titl = "--title $d[$#d] --outdir inv --pre 100 --len 200 --lowpass \"2 0.05\" --highpass \"2 0.02\" --range \"50 300\"";

 

 print "$string $data $epi $titl --pltname auto \n";
