#!/usr/bin/perl
# This script convert blastpgp output (m = 0 flag) to block of aligned texts

use Switch;
use strict;
use List::Util qw[min max];
use List::Util qw[sum];

my $targetfile = $ARGV[0];
my $outputfile = $ARGV[1];
my $rounds     = $ARGV[2]; 	# number of rounds of psibllast performed.

my $line_start = 0;
my $tgt_seq = "";
my $temp_seq = "";
my $temp_seq_clean = "";
my $start = -1;
my $tgt_start = 0;
my $tgt_length = 0;
our @Profile = ();

open (OUT, ">$outputfile") or die "can't create output file";
open (PRF, "$targetfile") or die "Can't open $targetfile file";
while (my $line = <PRF>){
	chomp $line;
	if ($line =~ /Results from round $rounds/){
		$line_start = 1;
	} # end if

	if ($line_start == 1){
		if ($line =~ /Identities = /){
			if ($start != -1){
				# clean up gaps 
				for (my $i=0; $i<length($tgt_seq); $i++){
					if (substr($tgt_seq, $i, 1) ne "-"){
						$temp_seq_clean = $temp_seq_clean.substr($temp_seq, $i, 1);
					} # end if
				} # end for

				# add dashes for nonaligned portions in the start and end
				for (my $i=1; $i<$tgt_start; $i++){
					$temp_seq_clean = "-".$temp_seq_clean;
				} # end for
				push (@Profile, $temp_seq_clean);		

				$temp_seq_clean = "";
				$tgt_seq = "";
				$temp_seq = "";
				$start = 0;
			} elsif ($start == -1){
				$start = 0;
			} # end if
		} # end if
		
		if ($line =~ /^Query/){
			my @columns = split (/\s+/, $line);
			$tgt_seq = $tgt_seq.$columns[2];
			if ($start == 0){
				$tgt_start = $columns[1];
				$start = 1;
			} # end if
		} elsif ($line =~ /^Sbjct/){
			my @columns = split(/\s+/, $line);
			$temp_seq = $temp_seq.$columns[2];
		}  # end if

		# get length of tgt seq
		if ($line =~ /^length of query/){
			my @columns = split(/:/, $line);
			$tgt_length = $columns[1];
		} # end if
	} # end if
} # end while
close PRF;

# clean up gaps for last alignment
for (my $i=0; $i<length($tgt_seq); $i++){
	if (substr($tgt_seq, $i, 1) ne "-"){
		$temp_seq_clean = $temp_seq_clean.substr($temp_seq, $i, 1);
	} # end if
} # end for

# add dashes for nonaligned portions in the start
for (my $i=1; $i<$tgt_start; $i++){
	$temp_seq_clean = "-".$temp_seq_clean;
} # end for
push (@Profile, $temp_seq_clean);               

# add dashes for nonaligned portions in the end
for (my $j=0; $j <= $#Profile; $j++){
	my $seq = $Profile[$j];
	for (my $i=length($Profile[$j]); $i< $tgt_length; $i++){
		$seq = $seq."-";
	} # end for
	$Profile[$j] = $seq;
	print OUT "$seq\n"
} # end for
close OUT;

