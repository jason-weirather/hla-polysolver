#!/usr/bin/perl

# usage: /cga/wu/sachet/hla/hla_caller/keep_only_read_pairs.pl hla_c_07_01_01_01.1960.temp.sam temp.hla_c_07_01_01_01.1960.temp.sam
# program assumes that the input sam file is sorted by read and has no header

$inFile = $ARGV[0];
$outFile = $ARGV[1];

open INFILE, $inFile || die "Cannot open $inFile\n";
open OUTFILE, ">$outFile" || die "Cannot open $outFile\n";

$count = 0;
while(1){
	$line1 = <INFILE>;
	$line2 = <INFILE>;
	#print "####\n$line1$line2";
	@f1 = split(/\t/,$line1);
	@f2 = split(/\t/,$line2);
	if($f1[0] eq $f2[0]){
		print OUTFILE $line1;
		print OUTFILE $line2;
	} else {
		seek(INFILE, -length($line2), 1);
	}
	if(eof){
		last;
	}
}


#while($line = <INFILE>){
#	chomp($line);
#	$count++;
#	@f = split(/\t/,$line);
#	$curRead = $f[0];
#	if($count == 1){
#		$prevRead = $curRead;
#		$prevLine = $line;
#		$curCount = 1;
#		next;
#	}
#	if($curRead eq $prevRead){
#		print OUTFILE $prevLine,"\n",$line,"\n";
#		$curCount++;
#	}
#	
#
#}


close INFILE;
close OUTFILE;




