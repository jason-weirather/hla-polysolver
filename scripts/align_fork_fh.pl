#!/usr/bin/env perl

# usage: $PSHOME/scripts/align_fork_fh.pl $outDir/merged.1.fastq $outDir/merged.2.fastq $NUM_THREADS $format $PSHOME/data/abc_complete.nix $outDir/nv.complete.chr6region.R0k6.sam 0 $NOVOALIGN_DIR


use Parallel::ForkManager;
use POSIX;

$f1 = $ARGV[0];
$f2 = $ARGV[1];
$nFork = $ARGV[2]; # should be a multiple of 4, usually 4 or 8
$format = $ARGV[3];
$nixFile = $ARGV[4];
$samPrefix = $ARGV[5];
$softClip = $ARGV[6];
$NOVOALIGN_DIR = $ARGV[7];
$digits = 2;

$n1=`wc -l $f1`;
chomp($n1);
@f = split(/ /,$n1);
$n1 = $f[0];

$n2=`wc -l $f2`;
chomp($n2);
@f = split(/ /,$n2);
$n2 = $f[0];

if($n1 != $n2){
	die "n1=$n1 is not equal to n2=$n2\n";
}


print "n1=$n1\tn2=$n2\n";

$nPerFile = $n1/$nFork;
$ceil = ceil($nPerFile);
if($ceil != $nPerFile){
	$nPerFile = ceil($nPerFile);
}

while($nPerFile % 4 != 0){
	$nPerFile = $nPerFile + 1;
}

print "n1=$n1\tn2=$n2\tnPerFile=$nPerFile\tceil=$ceil\n";
	
`split -l $nPerFile -d -a 2 $f1 $f1`;
`split -l $nPerFile -d -a 2 $f2 $f2`;

my $pm = Parallel::ForkManager->new($nFork);

for($i = 0; $i < $nFork; $i++){
	my $pid = $pm->start and next;
	$suffix = $i;
        if(length($i) < $digits){
                $diff = $digits - length($i);
                for($j = 0; $j < $diff; $j++){
                        $suffix = "0".$suffix;
                }
        }

	$novoalignPath = $NOVOALIGN_DIR."/novoalign";
	print "novoalignPath=$novoalignPath\n";

	if($softClip == 0){
		`time $novoalignPath -d $nixFile -f $f1$suffix $f2$suffix -F $format -R 0 -r all -o SAM -o FullNW | grep -P '\thla' > $samPrefix$suffix`;
	} else{
		`time $novoalignPath -d $nixFile -f $f1$suffix $f2$suffix -F $format -R 0 -r all -o SAM -g 20 -x 3 | grep -P '\thla' > $samPrefix$suffix`;
	}
	
	$pm->finish;
}

$pm->wait_all_children;

for($i = 0; $i < $nFork; $i++){
	$suffix = $i;
        if(length($i) < $digits){
                $diff = $digits - length($i);
                for($j = 0; $j < $diff; $j++){
                        $suffix = "0".$suffix;
                }
        }

	`cat $samPrefix$suffix >> $samPrefix`;

	`rm -rf $samPrefix$suffix $f1$suffix $f2$suffix`;
}


