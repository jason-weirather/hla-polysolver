#!/usr/bin/env perl

# usage: $PSHOME/scripts/convertToStandard.pl winners.hla.txt inferred.alleles.txt $PSHOME

$inFile = $ARGV[0];
$outFile = $ARGV[1];
$PSHOME = $ARGV[2];

$lib = $PSHOME."/scripts/common_functions_typing.pl";
require $lib;

open INFILE, $inFile || die "Cannot open $inFile\n";
open OUTFILE, ">$outFile" || die "Cannot open $outFile\n";

while($line = <INFILE>){
	chomp($line);
	@f = split(/\t/,$line);
	$allele4digit = get_4digit($f[1]);	
	$allele4digit =~ s/hla_a_/A*/;
	$allele4digit =~ s/hla_b_/B*/;
	$allele4digit =~ s/hla_c_/C*/;
	$allele4digit =~ s/_/:/;
	print OUTFILE "$allele4digit\n";
        $allele4digit = get_4digit($f[2]);
        $allele4digit =~ s/hla_a_/A*/;
        $allele4digit =~ s/hla_b_/B*/;
        $allele4digit =~ s/hla_c_/C*/;
        $allele4digit =~ s/_/:/;
	print OUTFILE "$allele4digit\n";
}

close INFILE;
close OUTFILE;
