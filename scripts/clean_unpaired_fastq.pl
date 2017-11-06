#!/usr/bin/env perl

# usage: $PSHOME/scripts/clean_unpaired_fastq.pl NA07048.1.fastq

$file = $ARGV[0];
$outFile = $file.".clean";

open FILE, $file;
open OUTFILE, ">$outFile";

while($line1 = <FILE>){
	$line2 = <FILE>;
	$line3 = <FILE>;
	$line4 = <FILE>;
	if($line1=~/\//){
		print OUTFILE "$line1$line2$line3$line4";
	}else{
		#print "$line1$line2$line3$line4";
	}


}

close FILE;
close OUTFILE;

system("mv $outFile $file");

