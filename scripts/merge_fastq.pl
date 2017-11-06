#!/usr/bin/env perl

#usage: $PSHOME/scripts/merge_fastq.pl $outDir/tag $outDir/chr6region $outDir/merged

# this program assumes that the fastq files have already been cleaned by clean_unpaired_fastq.pl

use lib "/cga/wu/sachet/software/lib/perl5/x86_64-linux";
use List::MoreUtils qw(uniq);
use List::Util 'max';
use List::Util 'min';
use Array::Utils qw(:all);

use Bio::DB::Sam;
use Dumpvalue;
use Data::Dumper;
use List::MoreUtils 'true';
use List::MoreUtils 'indexes';

$f1 = $ARGV[0];
$f2 = $ARGV[1];
$f3 = $ARGV[2];

# get all the first fastq reads

$file1 = $f1.".1.fastq";
$file2 = $f1.".2.fastq";

open FILE1, $file1 || die "Cannot open $file1\n";
open FILE2, $file2 || die "Cannot open $file2\n";

%hash1 = ();
@p1 = ();
while($line1 = <FILE1>){
	#print "line=$line1";
	$index = index($line1,"/1");
	$head = substr($line1,0,$index);
	push @p1,$head;
	$seq1 = <FILE1>;
	$third1 = <FILE1>;
	$qual1 = <FILE1>;
	$line2  = <FILE2>;
        $seq2 = <FILE2>;
        $third2 = <FILE2>;
        $qual2 = <FILE2>;
	$hash1{$head}->{seq1} = $seq1;
        $hash1{$head}->{third1} = $third1;
        $hash1{$head}->{qual1} = $qual1;
        $hash1{$head}->{seq2} = $seq2;
        $hash1{$head}->{third2} = $third2;
        $hash1{$head}->{qual2} = $qual2;
}

close FILE1;
close FILE2;

# get all the second fastq reads

$file1 = $f2.".1.fastq";
$file2 = $f2.".2.fastq";

open FILE1, $file1 || die "Cannot open $file1\n";
open FILE2, $file2 || die "Cannot open $file2\n";

%hash2 = ();
@p2 = ();
while($line1 = <FILE1>){
        $index = index($line1,"/1");
        $head = substr($line1,0,$index);
	push @p2,$head;
        #print "head=$head\n";
        $seq1 = <FILE1>;
        $third1 = <FILE1>;
        $qual1 = <FILE1>;
        $line2  = <FILE2>;
        $seq2 = <FILE2>;
        $third2 = <FILE2>;
        $qual2 = <FILE2>;
        $hash2{$head}->{seq1} = $seq1;
        $hash2{$head}->{third1} = $third1;
        $hash2{$head}->{qual1} = $qual1;
        $hash2{$head}->{seq2} = $seq2;
        $hash2{$head}->{third2} = $third2;
        $hash2{$head}->{qual2} = $qual2;
}

close FILE1;
close FILE2;

# get all the unique reads

@pIntersect = intersect(@p1,@p2);
@pDiff = array_diff(@p1,@p2);
@pBoth = (@pIntersect,@pDiff);

#@p = (@p1,@p2);
#@p = sort @p;
#@pBoth = ();
#for($i = 0; $i <= $#p; $i++){
#	if($i == 0){
#		$prev = $p[$i];
#		push @pBoth, $prev;
#		next;
#	}
#	if($p[$i] eq $prev){
#		next;
#	} else{
#		push @pBoth, $p[$i];
#		$prev = $p[$i];
#	}
#	
#}


# print the output
	
$file1 = $f3.".1.fastq";
$file2 = $f3.".2.fastq";

open FILE1, ">$file1" || die "Cannot open $file1\n";
open FILE2, ">$file2" || die "Cannot open $file2\n";

for($i = 0;$i <= $#pBoth; $i++){
	$read = $pBoth[$i];
	if(defined($hash1{$read})){
		print FILE1 "$read/1\n$hash1{$read}->{seq1}$hash1{$read}->{third1}$hash1{$read}->{qual1}"; 
		print FILE2 "$read/2\n$hash1{$read}->{seq2}$hash1{$read}->{third2}$hash1{$read}->{qual2}";
	
	}elsif(defined($hash2{$read})){
                print FILE1 "$read/1\n$hash2{$read}->{seq1}$hash2{$read}->{third1}$hash2{$read}->{qual1}";
                print FILE2 "$read/2\n$hash2{$read}->{seq2}$hash2{$read}->{third2}$hash2{$read}->{qual2}";
	}else{
		die "$read is undefined in both files?\n";
	}
}
