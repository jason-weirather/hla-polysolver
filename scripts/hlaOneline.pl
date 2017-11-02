#!/usr/bin/perl

#usage: /cga/wu/sachet/hla/hla_caller/capture/code_base_041713/hla_oneline_fh_050713.pl winners.R0k6.scale23.hla hla.intervals

$inFile = $ARGV[0];
$outFile = $ARGV[1];

open INFILE, $inFile || die "Cannot open $inFile";
open OUTFILE, ">$outFile" || die "Cannot open $outFile";

$line = <INFILE>;
chomp($line);
@f = split(/\t/,$line);
print OUTFILE $f[1],"\n",$f[2],"\n";

$line = <INFILE>;
chomp($line);
@f = split(/\t/,$line);
print OUTFILE $f[1],"\n",$f[2],"\n";

$line = <INFILE>;
chomp($line);
@f = split(/\t/,$line);
print OUTFILE $f[1],"\n",$f[2],"\n";

close INFILE;
close OUTFILE;
