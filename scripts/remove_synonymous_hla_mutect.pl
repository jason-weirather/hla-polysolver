#!/usr/bin/perl

#usage: ./remove_synonymous_hla_mutect.pl HNSC-TCGA-BA-6873/HNSC-TCGA-BA-6873.mutect.filtered.annotated HNSC-TCGA-BA-6873 HNSC-TCGA-BA-6873

$inFile = $ARGV[0];
$indiv = $ARGV[1];
$dir = $ARGV[2];
$outFileNonsyn = $dir."/".$indiv.".mutect.filtered.nonsyn.annotated";
$outFileSyn =  $dir."/".$indiv.".mutect.filtered.syn.annotated"; 

#@remove = ("p.R339*","p.T238I","p.A269A","p.T214T","p.T238T");
@remove = ("p.A269A","p.T214T","p.T238T");

%removeHash = ();

for($i = 0; $i <= $#remove; $i++){
	$removeHash{$remove[$i]} = 1;
	#print "remove: $remove[$i]\n";
}


open INFILE, $inFile || die "Cannot open $inFile\n";
open OUTFILENON, ">$outFileNonsyn" || die "Cannot open $outFileNonsyn\n";
open OUTFILESYN, ">$outFileSyn" || die "Cannot open $outFileSyn\n";

$header = <INFILE>;
chomp($header);

print OUTFILENON "$header\n";
print OUTFILESYN "$header\n";

while($line = <INFILE>){
	chomp($line);
	@f = split(/\t/,$line);
	$peptideChange = $f[-1];
	#print "$peptideChange\n";
	if($removeHash{$peptideChange} == 1){
		next;
	}
	$peptideChange = substr($peptideChange,2);
	$_ = $peptideChange;
	($a1,$pos,$a2) = /(\D+)(\d+)(\D+)/;
	#print "breakdown: $a1 $pos $a2\n";	
	if($a1 eq $a2){
		print OUTFILESYN "$line\n";
		next;
	}
	print OUTFILENON "$line\n";
}


close INFILE;
close OUTFILENON;
close OUTFILESYN;



