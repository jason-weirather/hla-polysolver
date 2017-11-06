#!/usr/bin/env perl

use Math::BaseCalc;

# usage: ./remove_sam_tag.pl LUSC-37-5819.hla_a_03_01_01_01.RG.bam LUSC-37-5819.hla_a_03_01_01_01.RG.rm.tag.bam 1

# change NotPrimaryAlignment tag in flags field  (flag at 2^8 position)
# change mapping quality to 70 (any non-zero value should be fine) field 5


$inFile = $ARGV[0];
$outFile = $ARGV[1];
$eventCountThreshold = $ARGV[2]; # reject reads where # events > this threshold. If you want to set this filter off, set it to a high number like 10

$calc = new Math::BaseCalc(digits => [0,1]); #Binary
$bin_string = $calc->to_base('7'); # Convert 465 to binary
$string = $calc->from_base('011');
#print "$bin_string\t$string\n";

$tempSam = $inFile.".temp.sam";

open FILE, ">$tempSam";

@lines = `/cga/wu/sachet/software/samtools/samtools view -h $inFile`;

for($i = 0; $i <= $#lines; $i++){
	#print "i=$i\n";
	$line = $lines[$i];
	chomp($line);
	if($line =~ /^\@/){
		print FILE "$line\n";
	}else{

		@f = split(/\t/,$line);
		$cigar = $f[5];
		($insertionCount,$deletionCount) = @{process_cigar($cigar)};
		#remove read if the two ends don't match to the same allele
		if(!($f[6] eq "=")){
			next;
		}

		# remove if the read has more than one event (mismatch, insertion, deletion)

		 $_ = $line;
                /(NM:i:\d+)/;
                $nm = $1;
                $insertionEventCount = @{[$cigar =~ /I/g]};
                $deletionEventCount = @{[$cigar =~ /D/g]};
                #$mismatchCount = @{[$md =~ /[A-Za-z]/g]};
                $editCount = substr($nm,5);
		$mismatchEventCount = $editCount - $insertionCount - $deletionCount;
		$eventCount = $insertionEventCount + $deletionEventCount + $mismatchEventCount;

		$_=$line;
		/(MD:Z:\w+)/;
		$md = $1;
		$md = substr($md,5);
		@mdArr = split(//,$md);

		#### CHECK ####

		#$eventCount = $editCount;
		
		#print "$line\neventCount=$eventCount\tinsertionEventCount=$insertionEventCount\tinsertionCount=$insertionCount\tdeletionEventCount=$deletionEventCount\tdeletionCount=$deletionCount\tmismatchEventCount=$mismatchEventCount\tnm=$nm\tcigar=$cigar\n";

		#### CHECK ####

		
		#### CHANGE ####
		#print "$eventCount\n";
		#next;
		#### CHANGE ####

		#print "$line\n";
		#print "$md\t$nm\t$cigar\t$insertionCount\t$deletionCount\t$mismatchCount\t$eventCount\n";

		if($eventCount > $eventCountThreshold){
			next;
		}

		$bin_string = $calc->to_base($f[1]);
		#print "f[1]=$f[1]\t$bin_string\n";
		@flags = split(//,$bin_string);
		#print "flags1 = @flags\n";
		if($#flags >= 8){
			$flags[-9] = 0;
		}
		#print "flags2 = @flags\n";
		$bin_string = join("",@flags);
		#print "$bin_string\n";
		$string = $calc->from_base($bin_string);
		#print "$f[1]\t$string\n";
		$f[1] = $string;
		$f[4] = 70;
		$line2 = join("\t",@f);
		print FILE "$line2\n";
	}

}

close FILE;

system("/cga/wu/sachet/software/samtools/samtools view -bS -o $outFile $tempSam");
system("rm -f $tempSam");

########

sub process_cigar { # start sub process_cigar

	my ($cigar) = @_;

	my ($i,@f,$iCount,$dCount,$c,$number,$type,@arr);

	# the cigar string is an alternating string of numbers and characters. read them in order.

	@f = split(//,$cigar);
	
	$iCount = 0;
	
	$dCount = 0;

	$number = ();

	for($i = 0; $i <= $#f; $i++){

		$c = $f[$i];

		if($c =~ /\d/){
			
			$number = $number.$c;

		} else{

			$type = $c;

			if($type eq "I"){

				$iCount = $iCount + $number;
			}

			if($type eq "D"){

				$dCount = $dCount + $number;

			}

			$number = ();

		}


	}	

	@arr = ($iCount,$dCount);

	return \@arr;

} # end sub process_cigar
