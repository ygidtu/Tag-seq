#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<5;

sub usage {
        my $usage = << "USAGE";

        Remove ODN tag for read2, paired-end.
        Author: zhoujj2013\@gmail.com
        Usage: $0 <R1.fq> <R2.fq> <ODN_tag> <outdir> <prefix>

        Fix: strip /1 /2 suffix when matching R1/R2 read IDs.
USAGE
print "$usage";
exit(1);
};


my $r1_f = shift;
my $r2_f = shift;
my $odn = shift;
my $outdir = shift;
my $prefix = shift;

-f $r1_f or die "[ERROR] R1 file not found: $r1_f\n";
-f $r2_f or die "[ERROR] R2 file not found: $r2_f\n";
-d $outdir or die "[ERROR] Output dir not found: $outdir\n";
length($odn) > 0 or die "[ERROR] ODN tag is empty\n";

open OUT1,">","$outdir/$prefix.rmODN.R1.fq" or die "Cannot write R1 output: $!";
open OUT2,">","$outdir/$prefix.rmODN.R2.fq" or die "Cannot write R2 output: $!";
open OUT3,">","$outdir/$prefix.rmODN.stat" or die "Cannot write stat: $!";

my %kept_pairs;  # key = base_id, value = 1 if both R1 and R2 were kept
my %filtered_R1; # key = base_id, value = 1 if R1 was filtered out
my %filtered_R2; # key = base_id, value = 1 if R2 was filtered out

print STDERR "[INFO] Processing R2: $r2_f\n";
open IN,"$r2_f" or die $!;
while(<IN>){
	my $id1 = $_;
	my $raw_id = $1 if(/(\S+)/);
	(my $base_id = $raw_id) =~ s/\/\d+$//;  # strip /1, /2 suffix
	my $seq = <IN>;
	my $id2 = <IN>;
	my $qual = <IN>;

	# Process R2 and determine if it should be kept
	my $start = 0;
	while($seq =~ /$odn/g){
		$start = pos($seq);
	}
	
	chomp($seq);
	chomp($qual);
	my $len = length($seq);
	my $new_seq = substr($seq,$start,$len-1);
	my $new_qual = substr($qual,$start,$len-1);

	if($start == 0){
		# ODN not found in sequence
		$filtered_R2{$base_id} = 1;
		next;
	}elsif(length($new_seq) > 0){
		(my $hdr = $id1) =~ s|/\d+$||;
		print OUT2 "$hdr";
		print OUT2 "$new_seq\n";
		print OUT2 "$id2";
		print OUT2 "$new_qual\n";
		$kept_pairs{$base_id} = 1;  # both R1 and R2 are kept
		$filtered_R2{$base_id} = 0;  # R2 was kept
	}else{
		$filtered_R2{$base_id} = 1;
		next;
	}
}
close IN;

my $read_total = 0;
my $read_with_odn = 0;

print STDERR "[INFO] Processing R1: $r1_f\n";
open IN,"$r1_f" or die $!;
while(<IN>){
	$read_total++;
	my $id1 = $_;
	my $raw_id = $1 if(/(\S+)/);
	(my $base_id = $raw_id) =~ s|/\d+$||;
	my $seq = <IN>;
	my $id2 = <IN>;
	my $qual = <IN>;

	if(exists $filtered_R2{$base_id} && $filtered_R2{$base_id} == 1){
		# R2 was filtered, skip R1 too
		$filtered_R1{$base_id} = 1;
		next;
	}else{
		$read_with_odn++;
		(my $hdr = $id1) =~ s|/\d+$||;
		print OUT1 "$hdr";
		print OUT1 "$seq";
		print OUT1 "$id2";
		print OUT1 "$qual";
	}
}
close IN;

print OUT3 "Raw flagment count: $read_total\n";
print OUT3 "Read count with ODN: $read_with_odn\n";
print OUT3 "R1 filtered (R2 removed): " . (scalar keys %filtered_R1) . "\n";
print OUT3 "R2 filtered (ODN not found or too short): " . (scalar keys %filtered_R2) . "\n";
print OUT3 "Kept pairs: " . (scalar keys %kept_pairs) . "\n";

close OUT1;
close OUT2;
close OUT3;

my $pct = $read_total > 0 ? sprintf("%.2f", $read_with_odn / $read_total * 100) : 0;
print STDERR "[INFO] Done: $read_with_odn / $read_total reads passed ODN filter ($pct%)\n";
