#!/usr/bin/env perl

# usage: $PSHOME/first_allele_calculations_fork.pl $PSHOME/data/ids $PSHOME $SAMTOOLS 8 $race $iFile $outDir

use Parallel::ForkManager;
use POSIX;

$idsFile = $ARGV[0];
$PSHOME = $ARGV[1];
$SAMTOOLS = $ARGV[2];
$nFork = $ARGV[3]; # should be a multiple of 4, usually 4 or 8
$race = $ARGV[4];
$iFile = $ARGV[5];
$outDir = $ARGV[6];

$n1=`wc -l $idsFile`;
chomp($n1);
@f = split(/ /,$n1);
$n1 = $f[0];

print "n1=$n1\tnFork=$nFork\n";

$nPerFile = $n1/$nFork;
print "nPerFile=$nPerFile\n";
$ceil = ceil($nPerFile);
if($ceil != $nPerFile){
        $nPerFile = ceil($nPerFile);
}

print "nPerFile=$nPerFile\n";

while($nPerFile % 4 != 0){
        $nPerFile = $nPerFile + 1;
}

print "n1=$n1\tnPerFile=$nPerFile\tnFork=$nFork\n";

open IDSFILE, $idsFile || die "Cannot open $idsFile\n";

@ids = ();
while($line = <IDSFILE>){
	chomp($line);
	push @ids, $line;
}

close IDSFILE;

$index = 0;
for($j = 1; $j <= $nFork; $j++){
	if($index >= $n1){
		last;
	}
        $outFile = $outDir."/"."ids_$j";
        open OUTFILE, ">$outFile" || die "Cannot open OUTFILE\n";
	for($k = 0; $k < $nPerFile; $k++){
                if($index >= $n1){
                        last;
                }
                print "$index $ids[$index]\n";
                print OUTFILE $ids[$index],"\n";
                $index++;
        }
	close OUTFILE;
}	

	
$pm = Parallel::ForkManager->new($nFork);

for($j = 1; $j <= $nFork; $j++){
	$pid = $pm->start and next;
	$outFile = $outDir."/"."ids_$j";
	`$PSHOME/scripts/shell_first_allele_calculations $outFile $race $PSHOME $SAMTOOLS $iFile $outDir`;
	$pm->finish;
}

$pm->wait_all_children;

