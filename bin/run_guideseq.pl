#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<2;

sub usage {
        my $usage = << "USAGE";

        This script is designed for running Tag-seq analysis.
        Author: zhoujj2013\@gmail.com 2020/3/25
        Last update: 2020/9/29
        Usage: $0 <config.txt> <all|create_makefile|align|find_target>
               $0 <config.txt> resume [all|create_makefile|align|find_target]

        Use 'resume' to skip already-completed steps (safe to re-run after failure).

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my $step=shift;

my $resume = 0;
if($step eq "resume"){
    $resume = 1;
    $step = shift || "all";
}

my %conf;
&load_conf($conf, \%conf);

$conf{OUTDIR} = abs_path($conf{OUTDIR});
$conf{FORWARD_LIB_R1} = abs_path($conf{FORWARD_LIB_R1});
$conf{FORWARD_LIB_R2} = abs_path($conf{FORWARD_LIB_R2});
$conf{REVERSE_LIB_R1} = abs_path($conf{REVERSE_LIB_R1});
$conf{REVERSE_LIB_R2} = abs_path($conf{REVERSE_LIB_R2});
$conf{CTRL} = abs_path($conf{CTRL});
$conf{BIN} = abs_path($Bin);
$conf{GRNA} = abs_path($conf{GRNA});

#print Dumper(\%conf);

mkdir $conf{OUTDIR} unless(-d "$conf{OUTDIR}");

if($step eq "all" || $step eq "create_makefile"){

# Check if already done
my $mk1 = "$conf{OUTDIR}/$conf{PREFIX}_plus/makefile";
my $mk2 = "$conf{OUTDIR}/$conf{PREFIX}_minus/makefile";
if($resume && -f $mk1 && -f $mk2){
    print STDERR "[INFO] Makefiles already exist, skip creation (use 'create_makefile' to force re-generate).\n";
}else{
# Validate config keys
foreach my $k (qw(PREFIX OUTDIR FORWARD_LIB_R1 FORWARD_LIB_R2 FORWARD_LIB_TAG REVERSE_LIB_R1 REVERSE_LIB_R2 REVERSE_LIB_TAG GRNA INDEX REF CHROMSIZE MINLEN READLEN MAXINS THREAD ADAPTER BIN FASTQC STAR BEDTOOLS SAMTOOLS PICARD AdapterRemoval umi_tools water bedops)){
    defined $conf{$k} or die "[ERROR] Missing config key: $k\n";
}

# Validate input files
foreach my $k (qw(FORWARD_LIB_R1 FORWARD_LIB_R2 REVERSE_LIB_R1 REVERSE_LIB_R2)){
    -f $conf{$k} or die "[ERROR] $k file not found: $conf{$k}\n";
}
# create sample.lst
open OUT,">","$conf{OUTDIR}/sample.lst" || die $!;
print OUT "$conf{PREFIX}_plus\t$conf{FORWARD_LIB_R1}\t$conf{FORWARD_LIB_R2}\tforward\t$conf{FORWARD_LIB_TAG}\n";
print OUT "$conf{PREFIX}_minus\t$conf{REVERSE_LIB_R1}\t$conf{REVERSE_LIB_R2}\treverse\t$conf{REVERSE_LIB_TAG}\n";
close OUT;


# create config file
open OUT,">","$conf{OUTDIR}/config.txt" || die $!;
print OUT "#sample\n";
print OUT "SAMPLE\t$conf{OUTDIR}/sample.lst\n";
print OUT "OUTDIR\t$conf{OUTDIR}\n";

print OUT "# parameter for detecting potential cutting sites\n";
print OUT "MinSupportReadCount\t$conf{MinSupportReadCount}\n";
print OUT "MinCuttingEventCount\t$conf{MinCuttingEventCount}\n";

print OUT "# parameter for detecting off-target sites\n";
print OUT "MaxMismatch\t$conf{MaxMismatch}\n";
print OUT "MaxGap\t$conf{MaxGap}\n";
print OUT "MaxGapMismatch\t$conf{MaxGapMismatch}\n";

print OUT "# parameter\n";
print OUT "MINLEN\t$conf{MINLEN}\n";
print OUT "READLEN\t$conf{READLEN}\n";
print OUT "MAXINS\t$conf{MAXINS}\n";
print OUT "THREAD\t$conf{THREAD}\n";
print OUT "ADAPTER\t$conf{ADAPTER}\n";

print OUT "#program\n";
print OUT "BIN\t$conf{BIN}\n";
print OUT "FASTQC\t$conf{FASTQC}\n";
print OUT "STAR\t$conf{STAR}\n";
print OUT "BEDTOOLS\t$conf{BEDTOOLS}\n";
print OUT "SAMTOOLS\t$conf{SAMTOOLS}\n";
print OUT "PICARD\t$conf{PICARD}\n";
print OUT "AdapterRemoval\t$conf{AdapterRemoval}\n";
print OUT "umi_tools\t$conf{umi_tools}\n";
print OUT "water\t$conf{water}\n";
print OUT "bedops\t$conf{bedops}\n";

print OUT "# species\n";
print OUT "INDEX\t$conf{INDEX}\n";
print OUT "GENOME\t$conf{GENOME}\n";
print OUT "REF\t$conf{REF}\n";
print OUT "CHROMSIZE\t$conf{CHROMSIZE}\n";
close OUT;

	# run guide-seq
	chdir $conf{OUTDIR};
	print STDERR "[INFO] Running guideseq.v5.pl to generate Makefiles ...\n";
	my $guide_ret = system("perl $conf{BIN}/guideseq.v5.pl config.txt > $conf{PREFIX}.log 2>$conf{PREFIX}.err");
	if($guide_ret != 0){
		die "[ERROR] guideseq.v5.pl failed (exit=$guide_ret). Check:\n  tail -50 $conf{PREFIX}.err\n  tail -50 $conf{OUTDIR}/$conf{PREFIX}_plus/err\n";
	}
	print STDERR "[INFO] Makefiles generated successfully.\n";
	}
} # create makefile END. 

if($step eq "all" || $step eq "align"){
	chdir $conf{OUTDIR};
	
	foreach my $sample ("$conf{PREFIX}_plus", "$conf{PREFIX}_minus"){
		my $dir = "$conf{OUTDIR}/$sample";
		if(!-d $dir){
			print STDERR "[WARN] Sample dir not found: $dir (skip)\n";
			next;
		}
		if(!-f "$dir/makefile"){
			print STDERR "[WARN] Makefile not found: $dir/makefile (skip)\n";
			next;
		}
		
		# Resume: check if all steps already complete
		if($resume && -f "$dir/03visual.finished"){
			print STDERR "[INFO] $sample: all steps already complete, skip.\n";
			next;
		}
		
		print STDERR "[INFO] Running make for $sample ...\n";
		my $ret = system("cd $dir && make > log 2>err");
		if($ret != 0){
			print STDERR "[ERROR] make failed for $sample (exit=$ret). Check:\n";
			print STDERR "  tail -50 $dir/err\n";
			print STDERR "  tail -50 $dir/log\n";
			# Continue to allow other samples to complete
		}else{
			print STDERR "[INFO] make completed for $sample.\n";
		}
	}
	
	# get stat
	print STDERR "[INFO] Generating statistics ...\n";
	if(-d "$conf{OUTDIR}/$conf{PREFIX}_plus" && -d "$conf{OUTDIR}/$conf{PREFIX}_minus"){
		my $stat_ret = system("perl $conf{BIN}/getMatrics.pl $conf{PREFIX}_plus $conf{PREFIX}_minus > stat.txt 2>>$conf{PREFIX}.err");
		if($stat_ret != 0){
			print STDERR "[WARN] getMatrics.pl failed (exit=$stat_ret), some files may be missing.\n";
		}else{
			print STDERR "[INFO] Statistics written to $conf{OUTDIR}/stat.txt\n";
		}
	}
} # alignment END


if($step eq "all" || $step eq "find_target"){
	# detect off-targets for each sgRNAs
	open IN,"$conf{GRNA}" || die $!;
	while(<IN>){
		chomp;
		my @t = split /\s+/;
		my $id = $t[0];
		my $grna_seq = $t[1];
		my $pam = $t[2];
		$conf{PAM} = $pam;
		
		print STDERR "## Detect off-target sites for $id, sgRNA: $grna_seq, PAM:$pam. Start ... ##\n";

		# Resume: skip if already completed
		my $done_marker = "$conf{OUTDIR}/$id.find.target/$id.parsing_water_for_visualization.offtarget.bed";
		if($resume && -f $done_marker){
			print STDERR "[INFO] $id: already done (found $done_marker), skip.\n";
			print STDERR "## Detect off-target sites for $id, sgRNA: $grna_seq, PAM:$pam. End  ... ##\n";
			next;
		}

		mkdir "$conf{OUTDIR}/$id.find.target" unless(-d "$conf{OUTDIR}/$id.find.target");
		chdir "$conf{OUTDIR}/$id.find.target";

		if($conf{CTRL} eq "none"){
			`touch blank.bed`;
			$conf{CTRL} = "./blank.bed";
		}

		# create gRNA.fa
		open OUT,">","$conf{OUTDIR}/$id.find.target/grna.fa" || die $!;
		print OUT ">grna\n";
		print OUT "$grna_seq\n";
		close OUT;

		open OUT, ">", "run.sh" || die $!;
		print OUT "perl $conf{BIN}/find_targetsite.v2.pl ../$conf{PREFIX}\_plus/02potentialTargets/$conf{PREFIX}\_plus.plus.proximal.sorted ../$conf{PREFIX}\_plus/02potentialTargets/$conf{PREFIX}\_plus.minus.proximal.sorted ../$conf{PREFIX}\_minus/02potentialTargets/$conf{PREFIX}\_minus.plus.proximal.sorted ../$conf{PREFIX}\_minus/02potentialTargets/$conf{PREFIX}\_minus.minus.proximal.sorted $conf{CTRL} $conf{OUTDIR}/$id.find.target/grna.fa $conf{PAM} $conf{BIN}/../data/$conf{GENOME}.blacklist.bed $conf{REF} $conf{CHROMSIZE} $id $conf{MinSupportReadCount} $conf{MinCuttingEventCount} $conf{MaxMismatch} $conf{MaxGap} $conf{MaxGapMismatch}\n";
		close OUT;
		
		# Check prerequisite files exist before running
		my @prereqs = (
			"../$conf{PREFIX}_plus/02potentialTargets/$conf{PREFIX}_plus.plus.proximal.sorted",
			"../$conf{PREFIX}_plus/02potentialTargets/$conf{PREFIX}_plus.minus.proximal.sorted",
			"../$conf{PREFIX}_minus/02potentialTargets/$conf{PREFIX}_minus.plus.proximal.sorted",
			"../$conf{PREFIX}_minus/02potentialTargets/$conf{PREFIX}_minus.minus.proximal.sorted",
			$conf{REF},
			$conf{CHROMSIZE},
			"$conf{BIN}/../data/$conf{GENOME}.blacklist.bed"
		);
		my $all_ok = 1;
		foreach my $pf (@prereqs){
			if(!-f $pf && !-e $pf){
				print STDERR "[WARN] Prerequisite not found: $pf\n";
				$all_ok = 0;
			}elsif(-f $pf && !-s $pf){
				print STDERR "[WARN] Prerequisite is empty (0 targets): $pf\n";
				$all_ok = 0;
			}
		}
		unless($all_ok){
			print STDERR "[ERROR] Skipping $id: prerequisite files missing. Run 'align' step first.\n";
			next;
		}
		
		# Check if any targets exist - skip if empty
		my $total_targets = 0;
		foreach my $pf (@prereqs[0..3]){
			if(-f $pf && -s $pf){ $total_targets++; }
		}
		if($total_targets == 0){
			print STDERR "[WARN] $id: no potential targets found in alignment output, skip.\n";
			next;
		}
		
		print STDERR "[INFO] Running find_targetsite for $id ...\n";
		my $ft_ret = system("perl $conf{BIN}/find_targetsite.v2.pl " .
			"../$conf{PREFIX}_plus/02potentialTargets/$conf{PREFIX}_plus.plus.proximal.sorted " .
			"../$conf{PREFIX}_plus/02potentialTargets/$conf{PREFIX}_plus.minus.proximal.sorted " .
			"../$conf{PREFIX}_minus/02potentialTargets/$conf{PREFIX}_minus.plus.proximal.sorted " .
			"../$conf{PREFIX}_minus/02potentialTargets/$conf{PREFIX}_minus.minus.proximal.sorted " .
			"$conf{CTRL} $conf{OUTDIR}/$id.find.target/grna.fa $conf{PAM} " .
			"$conf{BIN}/../data/$conf{GENOME}.blacklist.bed $conf{REF} $conf{CHROMSIZE} " .
			"$id $conf{MinSupportReadCount} $conf{MinCuttingEventCount} " .
			"$conf{MaxMismatch} $conf{MaxGap} $conf{MaxGapMismatch} " .
			"2>> $conf{OUTDIR}/$id.find.target/run.err");
		if($ft_ret != 0){
			print STDERR "[ERROR] find_targetsite.v2.pl failed for $id (exit=$ft_ret)\n";
			print STDERR "  Check: $conf{OUTDIR}/$id.find.target/run.err\n";
			next;
		}
		
		# check sites
		my $offtarget_bed = "../$id.parsing_water_for_visualization.offtarget.bed";
		unless(-f $offtarget_bed){
			print STDERR "[WARN] Offtarget BED not found for $id (may have no hits)\n";
		}
		
		mkdir "sites" unless(-d "sites");
		chdir "sites";
		my $mon_ret = system("perl $conf{BIN}/monitoring_offtargets.pl " .
			"../$id.parsing_water_for_visualization.offtarget.bed " .
			"../../$conf{PREFIX}_plus/03visual/$conf{PREFIX}_plus.plus.bdg " .
			"../../$conf{PREFIX}_plus/03visual/$conf{PREFIX}_plus.minus.bdg " .
			"../../$conf{PREFIX}_minus/03visual/$conf{PREFIX}_minus.plus.bdg " .
			"../../$conf{PREFIX}_minus/03visual/$conf{PREFIX}_minus.minus.bdg " .
			"> monitoring_offtargets.log 2>monitoring_offtargets.err");
		if($mon_ret != 0){
			print STDERR "[WARN] monitoring_offtargets.pl failed for $id (exit=$mon_ret)\n";
		}
		chdir "..";
		
		# check global site distribution (remove this part, 2020.09.29)
		#mkdir "global" unless(-d "global");
		#chdir "global";
		#if($conf{GENOME} eq "hg19"){
		#	`cp $conf{BIN}/RIdeogram/data/human_karyotype.hg19.txt ./karyotype.txt`;
		#	`cp $conf{BIN}/RIdeogram/data/human_genedensity.hg19.txt ./genedensity.txt`;
		#}elsif($conf{GENOME} eq "hg38"){
		#	`cp $conf{BIN}/RIdeogram/data/human_karyotype.hg38.txt ./karyotype.txt`;
		#	`cp $conf{BIN}/RIdeogram/data/human_genedensity.hg38.txt ./genedensity.txt`;
		#}
		#`perl $conf{BIN}/RIdeogram/prepare_target_offtarget_for_Rideogram.pl ../$id.all.sites.merged.confirmed ../$id.parsing_water_for_visualization.offtarget.bed > ./offtargets.txt`;
		#
		#open UL,">","./run.sh" || die $!;
		#print UL "ulimit -s 16384\n";
		#print UL "Rscript $conf{BIN}/RIdeogram//run_RIdeogram.R offtargets.txt $id > RIdeogram.log 2>RIdeogram.err\n";
		#close UL;
		#`sh ./run.sh`;
		#chdir "..";
		print STDERR "## Detect off-target sites for $id, sgRNA: $grna_seq, PAM:$pam. End  ... ##\n";
	}
}# find target END

chdir "..";
chdir "..";

#########################
sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}
