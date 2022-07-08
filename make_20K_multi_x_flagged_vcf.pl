#!/usr/bin/perl

use strict;
use 5.010;
#use Statistics::R;

if($#ARGV != 0){
	print "QC-version 07/08/2022 10:29am by Mike Schmidt, mschmidt\@med.miami.edu\n";
	die "args: <control file>\n";
	
}

=begin
VFLAGS
0	Passed
1	Variant failed preliminary QC: GATK “FILTER” ≠ “PASS” or is in tranche ≥ 99.7%
2	Variant failed preliminary QC: All genotypes have DP<10 and/or GQ<20
3	Monomorphic
4	Call Rate ≤ 80%
5	Average mean depth > 500 reads
7	male-het > <user defined>
11	outside caputre region
=cut

my($ctrlfile, $temp, $names, $vcf, $shortest_seg);
my($ref1, $ref2, $ref3, $ref4, $ref5, $ref6, $ref7, $ref8, $ref9, $ref10, $ref11, $ref12, $ref13, $ref14, $ref15, $ref16);				# references returned from functions
my($qc_grp_cnt, $race_grp_cnt);																	# number of QC/RACE-groups
my(@vcfs);																						# holds VCF input file name(s)
my(@race_groups);																				# groups listed in column RACE of the .fam file (DutchIsolate, Hispanic, Caucasian, ...)
my(@qc_groups);																					# groups listed in column QC_GROUP (ADNI, FAM, CC)
my(@sex);																						# sex listed in column 6 of the .fam file. 0 = male, 1 = female. all set to female chr1-22
my(@real_sex);																					# unmodified
my(@sample_ids);																				# The ID-names from column 2 of .fam file and in the VCF file
my(@neworder);																					# numeric odering of IDs
my(@keep);																						# keep (=1) for QC or discard (=0)
my(@par_index);																					# [][0] = father index, [][1] = mother index
my(%index);																						# numeric index for IDs in VCF
my(%race_assign);																				# key = RACE (Hipanic, DutchIsolate, ...) value = number 0, ..., N
my(%qc_assign);																					# key = QC_GROUP (ADNI, CC, FAM) value = number 0, 1, 2
my(%parameters);																				# parameters from control file and defau;t values
my($i, $j);
my(@race_names);																				# ordered list of RACE names
my(@qc_names);																					# ordered list of QC names
my(@aff_stat);																					# affection statu. 0 = unaff, 1 = aff, less 0 is unknown
my(@usehwe);
my(%capture);																					# capture kit name of key= qc_group value = map name
my(@qc_group_map);																				# which map-vlue to use for the qc_lookup
$ctrlfile = $ARGV[0];
$ref1 = get_data($ctrlfile);																	# read the control file and fill in defaults as needed
%parameters = %$ref1;
my($debug_info) = "dbg_" . $parameters{"vcf"};
our($dbg_handle);
open($dbg_handle,">$debug_info");
display_parameters(\%parameters);																# print to stdout what was entered or used as default
$vcf =  $parameters{"vcf"};																		# the raw VCF data file
my($id_file) = $parameters{"id_file"};															# the .fam file. It is expected to have the same IDs in the exact same order as the VCF


print stderr "reading $id_file\n";

($ref1, $ref2, $ref3, $ref4, $ref5, $ref6, $ref7, $ref8, $ref9, $ref10, $ref11, $ref12, $race_grp_cnt, $qc_grp_cnt, $ref13, $ref14, $ref15, $ref16) = read_pedfile($id_file, $vcf, \%parameters); 
																								# see function body for details
@sex = @$ref1;
@par_index = @$ref2;
@keep = @$ref3;
@race_groups = @$ref4;
@sample_ids = @$ref5;																			# the SAMPLE_IDs as they occur in the VCF file
%index = %$ref6;																				# key = SAMPLE_ID, value = position = 0, 1, ..., N
%race_assign = %$ref7;
@qc_groups = @$ref8;
%qc_assign = %$ref9;
@race_names = @$ref10;
@qc_names = @$ref11;
@aff_stat = @$ref12;
@usehwe = @$ref13;
%capture = %$ref14;
my(@male_thresholds) = @$ref15;
@real_sex = $ref16;
my($key, $value);
my($target_size, $chr, $bp_cnt);
my(%bp_hash);
my(@bp_lst);
my(@endpoints);
if($parameters{"read_capture"}){
	print stderr "reading $vcf 1st pass\n";
	($ref1, $ref2, $ref3, $bp_cnt) = get_vcf_positions($vcf);
	print $dbg_handle "\n$bp_cnt positions\n\n";
	%bp_hash = %$ref1;
	@bp_lst = @$ref2;
	@endpoints = @$ref3;
}
while(($key, $value) = each(%qc_assign)){
	$qc_names[$value] = $key; 
}

my(%capture_lookup);
my(@map_lst);
if($parameters{"read_capture"}){
	($ref1, $ref2) = get_capture_relations(\%capture, \%qc_assign, $parameters{"capture_relations"});
	@qc_group_map = @$ref1;
	@map_lst = @$ref2;

	$ref1 = read_capture_maps(\@bp_lst, \@endpoints, \%bp_hash, \@map_lst, $parameters{"chr"});

	%capture_lookup = %$ref1;
}

while(($key, $value) = each(%race_assign)){
	$race_names[$value] = $key; 
}


$ref1 = qc_group_cnt(\@keep, \@race_groups, \@qc_groups, \@race_names, \@qc_names, \@sex);
print "\nQC-subsets:\n";
for($i = 0; $i < @qc_names; $i++){
	print "$i.\t$qc_names[$i]\n";
}

print "\nEthnicities:\n";
for($i = 0; $i < @race_names; $i++){
	print "$i.\t$race_names[$i]\n";
}

print stderr "reading $vcf final pass\n";
read_vcf($vcf, \%parameters, \@keep, \@sample_ids, \@sex, \@race_groups, \@qc_groups, \@par_index, \%index, \@race_names, \@qc_names, $ref1, $race_grp_cnt, $qc_grp_cnt, \@aff_stat, \@usehwe, \%capture_lookup, \@qc_group_map, $target_size, $shortest_seg, \@male_thresholds, @real_sex);
close($dbg_handle);


##################################
# this function reads the VCF ad does all the QC. It's a monster
sub read_vcf
{
	my($infile) = $_[0];
	my($ref1) = $_[1];
	my($ref2) = $_[2];
	my($ref3) = $_[3];
	my($ref4) = $_[4];
	my($ref5) = $_[5];
	my($ref6) = $_[6];
	my($ref7) = $_[7];
	my($ref8) = $_[8];
	my($ref9) = $_[9];
	my($ref10) = $_[10];
	my($ref11) = $_[11];
	my($race_grp_cnt) = $_[12];
	my($qc_grp_cnt) = $_[13];
	my($ref12) = $_[14];
	my($ref13) = $_[15];
	my($ref14) = $_[16];
	my($ref15) = $_[17];
	my($target_size) = $_[18];
	my($shortest_seg) = $_[19];
	my($ref16) = $_[20];
	my($ref17) = $_[21];
	open(dbg,">debug.txt");
	my(@male_thresholds) = @$ref16;													# the MAX # of male hets has to be LESS-Equal than this number
	my($key, $val);
	my(%parameters) = %$ref1;														# program global parameters
	my(@keep) = @$ref2;																# the IDs we want to keep (=1) versus those to be skipped over (=0)
	my(@sample_ids) = @$ref3;														# sample IDs as found in the VCF
	my(@sex) = @$ref4;																# 0 = male, 1 = female
	my(@race_groups) = @$ref5;														# groups listed in column RACE of the .fam file (DutchIsolate, Hispanic, Caucasian, ...) but as numbers
	my(@qc_groups) = @$ref6;														# groups listed in column QC_GROUP (ADNI, FAM, CC) but as numbers
	my(@par_index) = @$ref7;														# parent_index [][0] = father-index [][1] = mother-index
	my(%index) = %$ref8;															# key = SAMPLE_ID value = position 0, 1, ..., N
	my(@race_names) = @$ref9;														# ordered list of RACE names
	my(@qc_names) = @$ref10;														# ordered list of QC names	
	my(@qc_group_cnt) = @$ref11;													# count of FEMALE races in the QC-groups. Used for hetz & pHWE calculations
	my(@affstat) = @$ref12;															# aff.status. 0 = unaff, 1 = aff, "< 0" = unkown/other dementia, etc.
	my(@usehwe) = @$ref13;															# use for HWE & zHet calculations
	my(%capture_lookup) = %$ref14;													# xor together over the various target maps, 2^0 + 2^1 + ... + 2^(N-1)
	my(@qc_group_map) = @$ref15;													# holds the 2^i value for the qc_group used in the capture_lookup
	my(@real_sex) = @$ref17;														# unmodified
	my(@genotypes);																	# PASSed genotypes for current SNV the genotype on a per sample basis. -1= missing; -10= failed; 0 = ref.ref; 1= ref,alt; 2 = alt,alt
	my(@race_gt_cnt_m, @race_gt_cnt_f);												# genotype counts on per RACE basis [0]= pass_RR, [1]= pass_RA, [2]= pass_AA, [3]= fail_RR, [4]= fail_RA, [5]= fail_AA, [6]= missing
	my(@race_gt_cnt_con_m, @race_gt_cnt_con_f);										# as above but only counts controls/unaff
	my(@qc_gt_cnt_m, @qc_gt_cnt_f);													# genotype counts on per QC_GROUP basis, ditto
	my($qc_id);																		# numeric value of QC cohort, -1 if cohort not listed in ctrl file
	my($race_id);																	# numeric value of RACE cohort, -1 if cohort not listed in ctrl file
	my($temp);																		# temp variable
	my($t1, $t2, $t3, $tt, $aaf);													# some temp variables for calculationg MAFs
	my($line);																		# line extracted from infile
	my($head, $info, $fields);														# leading 9 columns of SNV data
	my($i, $j, $k, $l, $m, $n, $o);													# loop counters
	my($ref, $var);																	# reference and variant allele
	my($bp, $gt, $dp, $ad, $gtnum);													# individual values to be extracted for each SNV and individual
	my($sum_gt, $p1gt, $p2gt);														# temp variables for mend. incons
	my($str);																		# geenric string variable
	my($endpoint, $curpt, $ontarget);												# indel to see if in capture region
	my($dbg1, $dbg2);																# debug variables
	my(@lform);																		# array that holds descriptor (8th column)
	my(@indivmaster);																# individual master array that keeps track of various counts over ALL positions
	my(@indivslave);																# individual master array that keeps track of various counts over CURRENT position. added to master depending on VFLAG
	my(@position_pass_type);														# to facilitate final indiv mess. (0,1,2) is pass, VFLAG!=11, fail
	my(@idindex);																	# IDs of individuals in order as they appear
	my(@het_cnt);																	# count od genotypes used for het-Z
	my(@headers);																	# the header in the companion file(s)
	my(@an_lst);																	# see which alleles are left AFTER QC, some weaker ones might drop out, update ALT column
	my($an);																		# allele count for INFO field AFTER QC
	my($indivmastersize);															# size of @indivmaster
	my($indexgt);																	# index for GT and DP and AD
	my($indexdp);																	# position of DP in the format field. Sometimes those chage for different positions
	my($indexad);																	# as above for AD
	my($indexgq);																	# as above for GQ
	my($pass_snv);																	# if entire SNV passed
	my($findex, $mindex);															# parent index
	my($ogt, $fgt, $mgt, $pgt);														# genotypes of : offspring, father, mother, parent (single parent case)
	my($fid, $mid);																	# father & mother index
	my($head, $tail);																# to split string with regex
	my($bp);																		# SNV base pair
	my($chr);																		# chromosome number
	my($qual);																		# SNV QUAL-score
	my($vtype);																		# variant type. "SNV", Insertion" or "Deletion"
	my(@lst);																		# array that holds strings of a line extracted from infile
	my($total);																		# total genotype calls
	my($called);																	# total called genotype calls
	my($altfail);																	# total 0/1 or 1/1 genotypes failed
	my($reffail);																	# total 0/0 genotypes fail
	my($passed);																	# genotypes PASSed
	my($num_alleles);																# variant allele count
	my(@mendincon);																	# mendelian incons. [0]= single-parent-off pairs with goog GT [1]= single-parent-off pairs incon [2]= single-parent-off pairs with non identical genotypes
																					# [3] = triads with good GT [4]= triads with incons
	my($missing);																	# missing calls per SNV
	my(@vflag);																		# flag returned, one for each QC_GROUP
	my(@snphwe);																	# one value per QC_group and RACE
	my(@hetz);																		# one value per QC_group and RACE
	my(@str_hwe, @str_hetz, @str_clean, @str_abh);									# formatted strings for summary file
	my($badcall);																	# exceeds missingness threshold
	my($hidepth);																	# exceeds depth threshold
	my(@variant_alleles);
	my($pass_str, $fail_str,$maf_str, $abh_str);									# for printing putposes to create comma separated string of pass/fail/MAF values 
	my(@badgtcnt);																	# FAILed genotypes count for current SNV
	my(@goodgtcnt);																	# PASSed genotypes count for current SNV
	my(@indiv);																		# holds individual call fields as read from VCF
	my($flagged_line);																# the tail we will print as output to flagged VCF
	my($mindp_female) = $parameters{"minDP_female"};								# min depth for females
	my($mindp_male) = $parameters{"minDP_male"};									# min depth for males
	my($mindp);																		# minDP for current sample
	my($mingq) = $parameters{"minGQ"};												# minimum PL-score
	my($minstat)  = $parameters{"minStat"};											# minimum to calculate HWE & hetz stats
	my($printflg) = $parameters{"printFLG"};										# print the flagged VCF (sometimes might only want SNV-summary)
	my($isx) = $parameters{"isX"};													# is chrX 9=1) or chr 1-22 (=0)
	my($use_phwex) = $parameters{"pHWEx"};											# 0 | 1 to use the R-HardyWeinberg on chrX
	my($use_capture) = $parameters{"read_capture"};									# 1 if targeted map, 0 for whole exome without target maps
	my($start_bp) = $parameters{"start_bp"};										# in case we don't want the entire VCF
	my($end_bp) = $parameters{"end_bp"};											# in case we don't want the entire VCF
	my($output_prefix) = $parameters{"output_prefix"};								# location where to write output. Same dir is default
	my($file_size) = -s $infile;													# size of input VCF
	my(@passes);																	# Does SNV pass or fail for any reason? One for each QC_GROUP
	my($hom1, $hom2, $het);
	my($targetnum);																	# bp specific target num with all QC-groups xor together
	my($currenttargetnum);															# QC-group specific target we want to lookup (power of 2)
	my($is_on_target);																# is the current bp on target for specific QC-group
	my($gq);																		# GQ-score
	my($ms);																		# MS score
	my($prematch, $postmatch);														# temp variables for string manipulation
	my($overall_maf, $overall_ac);													# string for overall MAF & AC calculation
	my($agt, $agt_m, $agt_f, $aa, $normal, $index);									# count of alt. genotypes. $normal = 1 -> MAF <= 0.5, $normal = 0 -> MAF > 0.5
	my($rgt_m, $rgt_f);																# count of ref. genotypes
	my(@maf);																		# MAF[qc_grp_cnt][allele], just raw counts of passing genotypes
	my(@mafp);																		# as above but as proportion
	my(@rmaf);																		# race & QC-group specific MAF
	my(@altgt);																		# ATL alleles for detection mon SNVs
	my($ref);																		# reference
	my(@gt);																		# valid genotype count
	my($z_het);																		# excess het. z-score
	my($good);																		# 1 if GT passed, 0 otherwise
	my($chr_name);																	# name of chr = {1,2,...,22,X,Y,...}
	my($infofix) = 0;																# need to add a bunch of stuff to the header relating to changes in INFO field
	my($pos);																		# position: <chr> <bp>
	my(@depth_sum);																	# sums depth of non-missing calls by QC_groups
	my(@depth_cnt);																	# count how many contribute to sum
	my(@dbg_sum);																	# sums depth of non-missing calls by QC_groups
	my(@dbg_cnt);																	# count how many contribute to sum
	my(@lstt);																		# temp list	
	my(@abhet_cnt);																	# [qc_names][allele_cnt][0]= ref [1]=alt
	my(@ad_lst);																	# allelic depth (AD) for ABHet calculation
	my($depth);																		# the quotient of sum / cnt
	my($passed);																	# total of PASSed calls
	my($cnt) = 0;																	# line count
	my($notpassed);																	# if the SNV as a whole has a "PASS" or "VQSRTrancheSNP..." which is FAIL
	my($validparents);																# number of valid parents
	my($flg);																		# handle for flagged VCF
	my(@pindex);																	# parent index array
	my(@lstgt);																		# genotypes listed. Only used for missingness calculation.
	my(@malehets, @validgt);														# male hets over the various QC-subgroups, and non "./." anysex-genotypes
	my($is_multal);																	# is multi. allelic?
	my($is_mono);																	# Is the SNV monomorphic? 1= yes, 0 = no
	my($sumcalls);																	# sum of all non-missing individual calls
	my($callrate);																	# = calls / (calls + missing)
	my($missthres) = $parameters{"min_call_rate"};									# we need to be > $missthres
	my($hidepththres) = $parameters{"maxDP"};										# maximum depth
	my($min_pass) = $parameters{"maxTranche"};										# 7th column
	my($debug) = $parameters{"debug"};												# debug will print genotypes on per line to qc_group specific .dbg file
	my(@dbg_out);																	# debug outfiles
	my($dbg_main);																	# all calls as listing
	my($filtered);																	# combo of BIPASS, Mono, depth, Badcall
	my($non_standard) = 0;															# anything that is not a single Ref and single ALt allele of exactly one base. (indels are non-standard and anything multi-allelic)
	my($hetz_maf);																	# MAF minus exclusions
	my($tranche);																	# 7th column of VCF
	my($non_missing);																# calls that are non-miising genotypes (not necessarily passing)
	my(@outstream);																	# SNV summary output
	my(@outnames);																	# names of companion files
	my($stream);																	# input from VCF
	my($usegz) = $parameters{"usegz"};												# use compressed or ASCII
	my($max_gt_cnt);																# theoretically possible max GT count
	my($allele_cnt);																# allele count (REF + ATLs) as listed in VCF file
	my(@alleles);																	# the alleles as integers, like 0 1
	my($missing_gt_index);															# index where to store missing GT counts
	my($mi_cnt);																	# count of SNVs MIs
	my($total_mi, $prop_mi);														# total MI-pairs, prop_mi = mi_cnt / total_mi
	my($het_obs, $het_exp, $incr);													# needed for zHet calculation
	my($rsid, $fout, $hidp);														# for SNV-summary file, rsID, flitered-out, high-depth
	my($gtf, $mono, $crate, $cb, $gpass, $itr, $miss);								# various silly stats for the SNV-summary file
	my($num, $den, $clean_str);														# numerator and denominator for callrate calculation
	my($sample_cnt) = 0;															# total samples
	my($ontarget_cnt);																# how many samples are on target at particular position
	my(@subset_cnt);																# counts of samples in each subset
	my($const_max_gt) = 78;															# allow for 78 (0..77) genotypes (= 12 alleles). Otherwise group into [.][78]
	my($m_ref, $m_alt);																# chrX male genotype counts, REF versus ALT (combined if more than 1 ALT)
	my($chr23) = 0;																	# pull out PAR and relabel chr23 into seperate VCF
	my($line_cnt) = 0;																# to give progress indicator
	$start_bp += 0;
	$end_bp += 0;

=begin
	my($R) = Statistics::R->new();													# to perform R calculations within perl
	if($use_phwex){
		$R->startR;																	# start R
		$R->send('library(HardyWeinberg);');										# load relevant R libraries. Those need to be installed on the system, otherwise it will fail.
	}
=cut
	my(@ps_auto);																	# 2 pseudo-autosomal regions to be excluded

	if($infile =~ m/gz$/){															# open .gz or uncompressed text file.vcf
		open($stream, "-|", "gzip -dc $infile") or die "Can not open infile $infile\n";
		$usegz = 1;
	}
	else{
		open($stream,"<$infile") or die "Can not open infile $infile\n";
	}
	if($isx){
		@ps_auto = split(/\s+/, $parameters{"PAR"});
		for($i = 0; $i < @ps_auto; $i++){
			$ps_auto[$i] += 0;
		}
		undef($chr23);
		if($usegz){
			unlink("chr23.vcf.gz");
			open($chr23, "| bgzip -c > chr23.vcf.gz");
		}
		else{
			open($chr23, ">chr23.vcf");
		}
	}
	my(%titv);																		# holds all TITv combinations. 5 = Ti and 6 = Tv
	$titv{"AG"} = "5";
	$titv{"GA"} = "5";
	$titv{"TC"} = "5";
	$titv{"CT"} = "5";
	$titv{"TG"} = "6";
	$titv{"GT"} = "6";
	$titv{"GC"} = "6";
	$titv{"CG"} = "6";
	$titv{"TA"} = "6";
	$titv{"AT"} = "6";
	$titv{"CA"} = "6";
	$titv{"AC"} = "6";
	
	
	
	$hidepththres += 0;
	$missthres += 0;
	$minstat += 0;
	for($i = 0; $i < @qc_names; $i ++){
		$subset_cnt[$i] = 0;
		for($j = 0; $j < @race_names; $j++){
			$subset_cnt[$i] += $qc_group_cnt[$i][$j];
			$sample_cnt += $qc_group_cnt[$i][$j];
		}
	}
	my(%hz_excl);
	$temp = $parameters{"z_het_excl"};
	@lst = split(/\s+/, $temp);
	for($i = 0; $i < @lst; $i++){
		if($lst[$i] eq "0"){
			last;
		}
		$hz_excl{$lst[$i]} = 1;
	}
	
	$mindp += 0;	
	$mingq += 0;



	print STDERR "\nprocessing $infile\n";
	$infile =~ s/\.gz$//;
	print $dbg_handle "bp\tREF\tALT\tVTYPE\t@qc_names\n";	
	@lst = split(/\//, $infile);
	$info = $lst[$#lst];
	for($i = 0; $i < @qc_names; $i++){												# setup various headers needed for summary files
		if($isx){
			$headers[$i] = "CHR\tPOS\tPASS\tFAIL\tMissing\tGT_Failed\tClean\tMono\tCallRate\tCallBad\tGATKPass\tMAF\tMeanDepth\tHiDepth\tABHet\tMend_Incon\tMend_pairs\tpropMI\tMultiAllele\tFilteredOut\tVFLAGS\trsID\tRefAllele\tAltAllele\tQUAL\tFILTER\tVTYPE\tInTargetRegion\tMaleHet";
		}
		else{
			$headers[$i] = "CHR\tPOS\tPASS\tFAIL\tMissing\tGT_Failed\tClean\tMono\tCallRate\tCallBad\tGATKPass\tMAF\tMeanDepth\tHiDepth\tABHet\tMend_Incon\tMend_pairs\tpropMI\tMultiAllele\tFilteredOut\tVFLAGS\trsID\tRefAllele\tAltAllele\tQUAL\tFILTER\tVTYPE\tInTargetRegion";		
		}
		for($j = 0; $j < @race_names; $j++){
			if($qc_group_cnt[$i][$j] >= $minstat){
				$headers[$i] .= "\t" . "nClean_" . $race_names[$j];
			}
		}
		$headers[$i] .= "\t";
		for($j = 0; $j < @race_names; $j++){
			if($qc_group_cnt[$i][$j] >= $minstat){
				$headers[$i] .= "Zhet_" . $race_names[$j] . ",";
			}			
		}
		$headers[$i] =~ s/,$//;
		$headers[$i] .= "\t";
		for($j = 0; $j < @race_names; $j++){
		
			if($qc_group_cnt[$i][$j] >= $minstat){
				$headers[$i] .= "pHWE_" . $race_names[$j] . ",";
			}
		}	
		$headers[$i] =~ s/,$//;
		$outnames[$i] = $output_prefix . $qc_names[$i] . $info . "." . $parameters{"snv_summary"};
		open($outstream[$i],">$outnames[$i]");
		$temp = $outstream[$i];
		print $temp "$headers[$i]\n";
		if($debug){
			$temp = $qc_names[$i] . ".dbg";
			open($dbg_out[$i],">$temp");
		}
	} # for($i = 0; $i < @qc_names; $i++){
	
	$temp =  $output_prefix . "flagged_" . $info;													# outfile.vcf, either compressed or uncompressed just like input.vcf
#	die "output= $temp\ninfo= $info\n";
	if($printflg){
		if($usegz){
			open($flg, "| bgzip -c > $temp.gz");
		}
		else{
			open($flg,">$temp");
		}
	}

	while(chomp($line = <$stream>)){
		if($line =~ m/^#/){
			if($chr23){
				print $chr23 "$line\n";
			}
			if($line =~ m/^#CHROM/){
				$indivmastersize = 0;
				$line =~ s/^(\S+\s+){9,9}//;
				$head = $&;
				
				@lst = split(/\s+/, $line);

				@idindex = split(/\s+/, $names);
				$names = "";
				$m = 0;
				$n = 0;
# all the individual stats are not used in the multi-allelic case but probably will at some point in the future. So they stay in the code for now.
				for($j = 0; $j < @sample_ids; $j++){									# init individual array
					$indivmaster[$j][0] = $sample_ids[$j];								# name
					$indivmaster[$j][1] = 0;											# missing genotypes
					$indivmaster[$j][2] = 0;											# singleton
					$indivmaster[$j][3] = 0;											# private doubleton
					$indivmaster[$j][4] = 0;											# doubleton
					$indivmaster[$j][5] = 0;											# Ti
					$indivmaster[$j][6] = 0;											# Tv
					$indivmaster[$j][7] = 0;											# het. genotypes
					$indivmaster[$j][8] = 0;											# alt.hom. genotypes
					$indivmaster[$j][9] = 0;											# sum depth of "PASS" genotypes
					$indivmaster[$j][10] = 0;											# parent count
					$indivmaster[$j][11] = 0;											# 1P MI count
					$indivmaster[$j][12] = 0;											# 2P MI count
					$indivmaster[$j][13] = 0;											# valid 1P counts
					$indivmaster[$j][14] = 0;											# valid 2P counts
					$indivmaster[$j][15] = 0;											# indel count
					$indivmaster[$j][16] = 0;											# valid genotype count
					$indivmaster[$j][17] = 0;											# pre-QC 0/0
					$indivmaster[$j][18] = 0;											# pre-QC 0/alt
					$indivmaster[$j][19] = 0;											# pre-QC alt/alt
					$indivmaster[$j][20] = 0;											# pre-QC 0/0 set ./.
					$indivmaster[$j][21] = 0;											# pre-QC 0/alt set ./.
					$indivmaster[$j][22] = 0;											# pre-QC alt/alt set ./.

					$indivslave[$j][0] = $sample_ids[$j];								# name
					$indivslave[$j][1] = 0;												# missing genotypes
					$indivslave[$j][2] = 0;												# singleton
					$indivslave[$j][3] = 0;												# private doubleton
					$indivslave[$j][4] = 0;												# doubleton
					$indivslave[$j][5] = 0;												# Ti
					$indivslave[$j][6] = 0;												# Tv
					$indivslave[$j][7] = 0;												# het. genotypes
					$indivslave[$j][8] = 0;												# hom. genotypes
					$indivslave[$j][9] = 0;												# sum depth of "PASS" genotypes
					$indivslave[$j][10] = 0;											# parent count
					$indivslave[$j][11] = 0;											# 1P MI count
					$indivslave[$j][12] = 0;											# 2P MI count
					$indivslave[$j][13] = 0;											# valid 1P counts
					$indivslave[$j][14] = 0;											# valid 2P counts
					$indivslave[$j][15] = 0;											# indel count
					$indivslave[$j][16] = "";											# genotype	
					
					$indivmastersize++;
					if($keep[$j]){
						$names .= $sample_ids[$j] . "\t";
					}
				}
				$names =~ s/\s+$//;
				if($printflg){
					print $flg "$head$names\n";
				}
			}
			else{
				if($printflg){
					
					if(($infofix == 0) && ($line =~ m/^##INFO=/)){
						$str = fix_header(\@qc_names);
						print $flg "$str";
						$infofix = 1;
					}
					print $flg "$line\n";
				}
			}
			next;
		} # if($line =~ m/^#/)
		if(($start_bp > 0) && ($usegz == 0)){
			$stream = find_start_bp($stream, $start_bp, $file_size, $usegz);
			$start_bp = 0;
		}
		elsif(($start_bp > 0) && ($usegz == 1)){
			$line =~ m/^\S+\s+\S+/;
			$temp = $&;
			$temp =~ m/\S+$/;
			$bp = $&;
			$cnt++;
	
			if($start_bp > $bp){
				next;
			}
		}
#		$line =~ s/^chr//;																
		$line =~ m/^\S+/;
		$chr_name = $&;
		if($parameters{"chr"} ne $chr_name){	
			print "unexpected chr_name= $chr_name\n";
			next;
		}
		$line =~ s/^(\S+\s+){7,7}//;
		$head = $&;
		@lst = split(/\s+/, $head);
		$chr = $lst[0];

		$bp = $lst[1];
		$rsid = $lst[2];
		$ref = $lst[3];
		$var = $lst[4];
		$qual = $lst[5];
		$tranche = $lst[6];
		$non_standard = 0;
		$bp += 0;
		if($end_bp < $bp){
			last;
		}
		if($isx){
			$temp = 0;
			for($i = 0; $i < @ps_auto; $i += 2){
				if(($bp >= $ps_auto[$i]) && ($bp <= $ps_auto[$i+1])){
					$temp = 1;
					last;
				}
			}
			if($temp){																					# is in PAR, next position
				$temp = $head;
				$temp =~ s/^\S+//;
				$temp = "chrX" . $temp . $line;
				print $chr23 "$temp\n";
				next;
			}
		}
#		if($isx && ($ps_auto[0] <= $bp && $ps_auto[1] >= $bp) || ($ps_auto[2] <= $bp && $ps_auto[3] >= $bp) || ($ps_auto[4] <= $bp && $ps_auto[5] >= $bp)){
#			next;
#		}		
		$temp = $var;
		$temp =~ s/,//g;
		$allele_cnt = length($var) - length($temp) + 2;													# figure out how many ALT alleles
		$max_gt_cnt = $allele_cnt * ($allele_cnt+1) / 2;												# calculate possibl genotype counts
		$missing_gt_index = 2 * $max_gt_cnt;
		$ms = 1;

		$line =~ s/^\S+//;
		$info = $&;
		$line =~ s/^\s+//;
		$line =~ s/^\S+//;
		$fields = $&;
		@lform = split(/:/, $fields);
		$line =~ s/^\s+//;
		@lst = split(/\s+/, $line);	

		$indexgt = -1;																					# index for GT and DP and AD
		$indexdp = -1;
		$indexad = -1;
		$indexgq = -1;
		for($i = 0; $i < @lform; $i++){																	# fields are not always in same order or some are extra/missing, so have to test each position
			if($lform[$i] eq "GT"){
				$indexgt = $i;
			}
			elsif($lform[$i] eq "DP"){
				$indexdp = $i;
			}
			elsif($lform[$i] eq "AD"){
				$indexad = $i;
			}
			elsif($lform[$i] eq "GQ"){
				$indexgq = $i;
			}
		} # for($i = 0; $i < @lform; $i++)	

		$called = 0;
		$altfail = 0;
		$reffail = 0;
		$passed = 0;																					# total valid genotype calls

		for($i = 0; $i < 5; $i++){
			$mendincon[$i] = 0;
		}
		$missing = 0;															
		$non_missing = 0;
		$flagged_line = "";
			
		undef(@race_gt_cnt_m);																			# reset some arrays just in case
		undef(@race_gt_cnt_con_m);
		undef(@qc_gt_cnt_m);
		undef(@race_gt_cnt_f);																			# reset some arrays just in case
		undef(@race_gt_cnt_con_f);
		undef(@qc_gt_cnt_f);		
		undef(@abhet_cnt);
		for($i = 0; $i < $qc_grp_cnt; $i++){															# init various arrays that count qc_subgroup specific stuff
			for($m = 0; $m < ((2*$max_gt_cnt)+1); $m++){												# for: PASS, FAIL, missing								
				$qc_gt_cnt_m[$i][$m] = 0;
				$qc_gt_cnt_f[$i][$m] = 0;
			}
			$depth_sum[$i] = 0;
			$depth_cnt[$i] = 0;
#			$dbg_sum[$i] = 0;
#			$dbg_cnt[$i] = 0;			
			$altgt[$i] = 0;
			$malehets[$i] = 0;
			$validgt[$i] = 0;
			for($k = 0; $k < $race_grp_cnt; $k++){
				for($m = 0; $m < ((2*$max_gt_cnt)+1); $m++){
					$race_gt_cnt_m[$i][$k][$m] = 0;														# overall count
					$race_gt_cnt_f[$i][$k][$m] = 0;
				}
				for($m = 0; $m < $max_gt_cnt; $m++){
					$race_gt_cnt_con_m[$i][$k][$m] = 0;													# controls-only count (used in pHWE)
					$race_gt_cnt_con_f[$i][$k][$m] = 0;
				}
			}
			for($m = 0; $m < $allele_cnt; $m++){
				$abhet_cnt[$i][$m][0] = 0;
				$abhet_cnt[$i][$m][1] = 0;
			}
		} #for($i = 0; $i < $qc_grp_cnt; $i++)
		for($i = 0; $i < @sample_ids; $i++){
			$indivslave[$i][1] = 0;																		# missing genotypes
			$indivslave[$i][2] = 0;																		# singleton
			$indivslave[$i][3] = 0;																		# private doubleton
			$indivslave[$i][4] = 0;																		# doubleton
			$indivslave[$i][5] = 0;																		# Ti
			$indivslave[$i][6] = 0;																		# Tv
			$indivslave[$i][7] = 0;																		# het. genotypes
			$indivslave[$i][8] = 0;																		# hom. genotypes
			$indivslave[$i][9] = 0;																		# sum depth of "PASS" genotypes
			$indivslave[$i][10] = 0;																	# depth SNV count
			$indivslave[$i][11] = 0;																	# 1P MI count
			$indivslave[$i][12] = 0;																	# 2P MI count
			$indivslave[$i][13] = 0;																	# valid 1P counts
			$indivslave[$i][14] = 0;																	# valid 2P counts
			$indivslave[$i][15] = 0;																	# indel count
			$indivslave[$i][16] = -999;																	# genotype	
			$indivslave[$i][17] = -1;																	# allele 1
			$indivslave[$i][18] = -1;																	# allele 2
		}

		$agt = 0;																						# Alternative GenoTypes, everything valid that is not 0/0
		$agt_m = 0;
		$agt_f = 0;
		$rgt_m = 0;
		$rgt_f = 0;
		for($j = 0; $j < @lst; $j++){																	# now process one individual at a time
			$genotypes[$j] = -1000;																		# missing genotype
			$race_id = $race_groups[$j];
			$qc_id = $qc_groups[$j];
			if($keep[$j] == 0){
				next;
			}
			if($lst[$j] =~ m/^\./){																		# missing genotype
				$indivmaster[$j][1]++;																	# keep track of missing genotypes on a per individual basis
				if($sex[$j] == 0){
					$qc_gt_cnt_m[$qc_id][$missing_gt_index]++;
					$race_gt_cnt_m[$qc_id][$race_id][$missing_gt_index]++;
				}
				else{
					$qc_gt_cnt_f[$qc_id][$missing_gt_index]++;
					$race_gt_cnt_f[$qc_id][$race_id][$missing_gt_index]++;				
				}
				$flagged_line .=  "\t" . $lst[$j];														# 3 = missing genotype
				
				if($debug){
					$temp = $dbg_out[$qc_id];
					$dp = 0;
					$gq = 0;
					$gt = "./.";
				}
				next;
			}
			
			$called++;
			@indiv = split(/:/, $lst[$j]);
			
			$gt = $indiv[$indexgt];
			$dp = $indiv[$indexdp];
			$ad = $indiv[$indexad];
			$gq = $indiv[$indexgq];
			if($dp !~ m/^\d+/){
				$dp = 0;
			}
			$depth_sum[$qc_id] += $dp;
			$depth_cnt[$qc_id]++;
			$indivmaster[$j][9] += $dp;
			$indivmaster[$j][10]++;
			
			@alleles = split(/\/|\|/, $gt);																# split the standard 0/1 or phased 0|1 into 0 and 1 (or whatever alleles you have)
			$alleles[0] += 0;
			$alleles[1] += 0;
			if($sex[$j] == 0){
				$mindp = $mindp_male;
			}
			else{
				$mindp = $mindp_female;
			}

			$gtnum = ($alleles[0] * $allele_cnt) + $alleles[1] - (($alleles[0] * ($alleles[0] + 1)) / 2);							# turn $gt into numeric value
																										# Example with 5 alleles {0,1,2,3,4}. Require: allele[0] <= allele[1]
																										# 	0	1	2	3	4
																										#		5	6	7	8
																										#			9	10	11
																										#				12	13
																										#					14
																										#	Where (0,0) -> 0, (0,1) -> 1, ..., (1,1) -> 5, ...,(4,4) -> 14
																										#

			if(($sex[$j] == 0) && ($alleles[0] != $alleles[1]) && ($dp >= $mindp) && ($gq >= $mingq)){	# male hets -> FAIL
				if($alleles[0] == 0){
					$indivmaster[$j][21]++;
				}
				else{
					$indivmaster[$j][22]++;
				}
				$malehets[$qc_id]++;
				$validgt[$qc_id]++;
				$temp = $max_gt_cnt + $gtnum;
				$qc_gt_cnt_m[$qc_id][$temp]++;
				$race_gt_cnt_m[$qc_id][$race_id][$temp]++;
				$genotypes[$j] = -1*$gtnum - 1;															# -999 = male het on X
				$gt = "./.";
			}
			elsif(($dp >= $mindp) && ($gq >= $mingq)){													### PASS
				$indivslave[$j][17] = $alleles[0];															# might need those later for MI investigation
				$indivslave[$j][18] = $alleles[1];

				$indivmaster[$j][16]++;
				if(($alleles[0] == 0) && ($alleles[1] == 0)){
					$indivmaster[$j][17]++;
				}elsif($alleles[0] == 0){
					$indivmaster[$j][18]++;
				}
				else{
					$indivmaster[$j][19]++;
				}
				$validgt[$qc_id]++;
				$genotypes[$j] = $gtnum;																# update genotype
				if($sex[$j] == 0){
					$qc_gt_cnt_m[$qc_id][$gtnum]++;
					$race_gt_cnt_m[$qc_id][$race_id][$gtnum]++;
					if($alleles[0] > 0){
						$agt_m++;
					}
					else{
						$rgt_m++;
					}
				}
				elsif($sex[$j] == 1){
					$qc_gt_cnt_f[$qc_id][$gtnum]++;
					$race_gt_cnt_f[$qc_id][$race_id][$gtnum]++;	
					if($alleles[0] > 0){
						$agt_f++;
					}
					else{
						$rgt_f++;
					}
					if($alleles[1] > 0){
						$agt_f++;
					}
					else{
						$rgt_f++;
					}					
				}
				if(($affstat[$j] == 1) && ($usehwe[$j] == 1)){
					if($sex[$j] == 0){
						$race_gt_cnt_con_m[$qc_id][$race_id][$gtnum]++;
					}
					elsif($sex[$j] == 1){
						$race_gt_cnt_con_f[$qc_id][$race_id][$gtnum]++;
					}
				}
				if($gtnum > 0){
					$altgt[$qc_id]++;
					if($gtnum > $allele_cnt){
						$altgt[$qc_id]++;
					}
				}
				if($alleles[0] != $alleles[1]){															# ABHet
					@ad_lst = split(/,/, $ad);
					$m = $ad_lst[$alleles[0]];
					$n = $ad_lst[$alleles[1]];
					$abhet_cnt[$qc_id][$alleles[0]][0] += $m;
					$abhet_cnt[$qc_id][$alleles[0]][1] += $n;
					$abhet_cnt[$qc_id][$alleles[1]][0] += $n;
					$abhet_cnt[$qc_id][$alleles[1]][1] += $m;
					
				}
				if(($alleles[0] == 0) && ($alleles[1] > 0)){
					$indivmaster[$j][7]++;																# indiv het
				}
				elsif($alleles[0] != 0){
					$indivmaster[$j][8]++;																# indiv hom.
				}

			}
			else{																						### FAIL
				if(($alleles[0] == 0) && ($alleles[1] == 0)){
					$indivmaster[$j][20]++;
				}elsif($alleles[0] == 0){
					$indivmaster[$j][21]++;
				}
				else{
					$indivmaster[$j][22]++;
				}
				$validgt[$qc_id]++;
				$temp = $max_gt_cnt + $gtnum;
				if($sex[$j] == 0){
					$qc_gt_cnt_m[$qc_id][$temp]++;
					$race_gt_cnt_m[$qc_id][$race_id][$temp]++;
				}
				else{
					$qc_gt_cnt_f[$qc_id][$temp]++;
					$race_gt_cnt_f[$qc_id][$race_id][$temp]++;				
				}
				$genotypes[$j] = -1*$gtnum - 1;															# failing GTs 

				$gt = "./.";
			}
			$lst[$j] =~ s/^\S{3}//;
			$lst[$j] = $gt . $lst[$j];
			$flagged_line .= "\t" . $lst[$j];
		} # for($j = 0; $j < @lst; $j++){
		if(($agt_m + $agt_f) <= ($rgt_m + $rgt_f)){
			$agt = $agt_m + $agt_f;
			$normal = 1;
		}
		else{
			$agt = $rgt_m + $rgt_f;
			$normal = 0;
		}	
		$endpoint = $bp - 1 + length($ref);																# get the endpoints of the REF segment
		$curpt = $bp;
		if($use_capture){
			$targetnum = $capture_lookup{$bp};
		}

		### MAF for each QC-group and allele ###
		for($i = 0; $i <= @qc_names; $i++){																# init counts for MAF calculations
			for($j = 0; $j < $allele_cnt; $j++){
				$maf[$i][$j] = 0;
				$mafp[$i][$j] = 0;
			}
		}
		$n = $#qc_names + 1;																			# storage for the overall MAF/AC counts
		$l = 0;
		for($i = 0; $i < $allele_cnt; $i++){															# to get an overall allele count across all QC-subgroups
			$an_lst[$i] = 0;
		}

		for($i = 0; $i < @qc_names; $i++){																# step through all QC-subgroups for MAF purposes
			$total = 0;
			for($j = 0; $j < $max_gt_cnt; $j++){
				$total += 2 * $qc_gt_cnt_f[$i][$j];														# females x 2
				$total += $qc_gt_cnt_m[$i][$j];															# males
			}

			$index = 0;
			for($m = 0; $m < $allele_cnt; $m++){														# count the females in standard fashion
				for($o = $m; $o < $allele_cnt; $o++){													# note that the females get counted twice, once for each allele
					$maf[$i][$m] += $qc_gt_cnt_f[$i][$index];
					$maf[$i][$o] += $qc_gt_cnt_f[$i][$index];
					$maf[$n][$m] += $qc_gt_cnt_f[$i][$index];
					$maf[$n][$o] += $qc_gt_cnt_f[$i][$index];					
					$index++;
				}
			}
			if($isx){
				$index = 0;
				for($m = 0, $o = $allele_cnt; $m < $allele_cnt; $m++, $o--){							# only count the males once down th main diagonal
					$maf[$i][$m] += $qc_gt_cnt_m[$i][$index];
					$maf[$n][$m] += $qc_gt_cnt_m[$i][$index];
					$index += $o;
				}
			}			
			$l += $total;
			if($total > 0){
				for($j = 0; $j < $allele_cnt; $j++){													# MAF for the various alleles, including REF (=first)
					$an_lst[$j] += $maf[$i][$j];				
					$mafp[$i][$j] = $maf[$i][$j] / $total;
				}
			}
		}
		$an = 0;
		$temp = "";																						# this will be the new $var (ALT column) 
		@lstt = split(/,/, $var);
		for($i = 0; $i < $allele_cnt; $i++){															# look at AC (=allele count) and update ALT column
			if($an_lst[$i]){
				$an += $an_lst[$i];																		# 
			}
		}		
		
		$overall_maf = "";																				# strings that hold AF & AV to be inserted into info field
		$overall_ac = "";
		for($j = 1; $j < $allele_cnt; $j++){
			if($l > 0){
				$overall_ac .= "$maf[$n][$j],";
				$temp = $maf[$n][$j] / $l;
				$overall_maf .= "$temp,";
			}
			else{
				$overall_ac .= "0,";
				$overall_maf .= "0,";
			}
		}
		# AF = allele frequency
		# AC = allele count for each minor allele
		# AN = total allele count
		$overall_ac =~ s/,$//;																			# delete trailing ,
		$overall_maf =~ s/,$//;
		$info =~ m/AF=/;																				# break string to inset newly calculated AF & AC strings
		$prematch = $` . "AF=";
		$postmatch = $';
		$postmatch =~ m/;/;
		$postmatch = $';
		$info = $prematch . $overall_maf . ";" . $postmatch;											# patch the whole mess back together
		$info =~ m/AC=/;																				# ditto for AC
		$prematch = $` . "AC=";
		$postmatch = $';
		$postmatch =~ m/;/;
		$postmatch = $';
		$info = $prematch . $overall_ac . ";" . $postmatch;
		$info =~ m/AN=/;																				# ditto for AC
		$prematch = $` . "AN=";
		$postmatch = $';
		$postmatch =~ m/;/;
		$postmatch = $';
		$info = $prematch . $an . ";" . $postmatch;
		
		### end QC-group MAF
		$ontarget_cnt = 0;
		for($i = 0; $i < $qc_grp_cnt; $i++){															# step through all QC-groups for VFLAG purposes
			if($use_capture){
				$ontarget = $targetnum & $qc_group_map[$i];
			}
			else{
				$ontarget = 1;																			# without capture map everything is on target
			}
			if($ontarget){
				$ontarget_cnt += $subset_cnt[$i];
			}
			$depth = 0;
			if($depth_cnt[$i]){
				$depth = $depth_sum[$i] / $depth_cnt[$i];
			}
			for($m = 0; $m < ((2*$max_gt_cnt)+1); $m++){												# copy for the more simple minded get_vflag function	
				$lstgt[$m] = $qc_gt_cnt_m[$i][$m] + $qc_gt_cnt_f[$i][$m];								# only needed to calculate missing rate
			}
			($vflag[$i], $passes[$i]) = get_vflag(\@lstgt, $depth, \%parameters, $tranche, $altgt[$i], $ontarget, $allele_cnt, $malehets[$i], $male_thresholds[$i]);	# see what VFLAG that yields
			$str_hwe[$i] = "";																			# init various strings for the summary file
			$str_clean[$i] = "";
			$str_hetz[$i] = "";

			for($k = 0; $k < $race_grp_cnt; $k++){														# step through the races inside each QC-group
				if($qc_group_cnt[$i][$k] >= $minstat){													# if it has the min number of FEMALE samples
					$hom1 = $race_gt_cnt_con_f[$i][$k][0];												# 0/0 genotype
					$hom2 = 0;
					$het = 0;
					for($m = 1; $m < $allele_cnt; $m++){												# 0/X genotypes, X={1, 2, ..., max_ALT_alleles}
						$het += $race_gt_cnt_con_f[$i][$k][$m];
					}
					for($m = $allele_cnt; $m < $max_gt_cnt; $m++){
						$hom2 += $race_gt_cnt_con_f[$i][$k][$m];
					}
					if($isx && $use_phwex){
						$m_ref = $race_gt_cnt_con_m[$i][$k][0];											# male REF genotype on chrX
						$m_alt = 0;
=begin
						for($m = $allele_cnt, $n = $allele_cnt-1; $m < $max_gt_cnt; $m += $n, $n--){
							$m_alt += $race_gt_cnt_con_m[$i][$k][$m];
						}
						$temp = $m_alt + $m_ref + $het + $hom1 + $hom2;
						if($temp  >= $minstat){
							$R->set('a', $m_alt);
							$R->set('b', $m_ref);
							$R->set('aa', $hom2);
							$R->set('ab', $het);
							$R->set('bb', $hom1);
							$R->send('snv <- matrix(c(A = a, B = b, AA = aa, AB = ab, BB = bb), nrow=1)');
							$R->send('res <- HWExactStats(snv,x.linked=TRUE,verbose=TRUE)');
							$R->send('print(res)');
							$temp = $R->read;
							$temp =~ s/\n/\t/g;
							$temp =~ s/^\S+\s+//;
							$snphwe[$i][$k] = $temp;
						}
						else{
							$snphwe[$i][$k] = -1;
						}
=cut
					}
					else{
						$temp = $het + $hom1 + $hom2;
						if($temp >= $minstat){
							$snphwe[$i][$k] = snphwe($het, $hom1, $hom2);									# calculate pHWE
						}
						else{
							$snphwe[$i][$k] = -1;
						}
					}
					$het_exp = 1;
					$total = 0;
					for($m = 0; $m < $max_gt_cnt; $m++){
						$total += $race_gt_cnt_f[$i][$k][$m] + $race_gt_cnt_con_f[$i][$k][$m];
					}				
					$het_obs = $total;
					$incr = $allele_cnt;
					$m = 0;
					for($m = 0, $incr = $allele_cnt; $m < $max_gt_cnt; $m += $incr, $incr--){			# count everything on the main diagonal
						$het_obs -= ($race_gt_cnt_f[$i][$k][$m] + $race_gt_cnt_con_f[$i][$k][$m]);		# and subtract from total = sum of all hets
					}				
					for($m = 0; $m < $allele_cnt; $m++){
						$rmaf[$m] = 0;
					}

					$index = 0;
					for($m = 0; $m < $allele_cnt; $m++){												# MAF counts to calculate race specific chisq
						for($o = $m; $o < $allele_cnt; $o++){
							$rmaf[$m] += $race_gt_cnt_f[$i][$k][$index] + $race_gt_cnt_con_f[$i][$k][$index];
							$rmaf[$o] += $race_gt_cnt_f[$i][$k][$index] + $race_gt_cnt_con_f[$i][$k][$index];
							$index++;
						}
					}
					if($total > 0){
						for($m = 0; $m < $allele_cnt; $m++){
							$rmaf[$m] /= (2*$total);
							$het_exp -= $rmaf[$m] * $rmaf[$m];											# hets expected based on MAFs of all alleles
						}
					}
					$hetz[$i][$k] = 888888;																# init with something, should get overwritten
					if($total == $het_obs){
						$hetz[$i][$k] = 999999;															# the "all het" case
					}
					elsif($het_obs == 0){																# the "no het" and all missing case
						$hetz[$i][$k] = -999999;
					}
					elsif($total > 0){																	# the interesting case
						$het_obs /= $total;
						$temp = $het_exp - $het_obs;
						$hetz[$i][$k] = $total * $temp * $temp / ($het_obs * (1 - $het_obs));
						$hetz[$i][$k] = sqrt($hetz[$i][$k]);
						if($het_exp > $het_obs){
							$hetz[$i][$k] *= -1;
						}
					}
				
					$str_hetz[$i] .= sprintf("%7.5f,\t", $hetz[$i][$k]);								# add to output string
					$str_hwe[$i] .= sprintf("%g,\t", $snphwe[$i][$k]);


					if($isx){
						$str_clean[$i] .= sprintf("%i,",$race_gt_cnt_m[$i][$k][0]);
						for($m = $allele_cnt, $n = $allele_cnt-1; $m < $max_gt_cnt; $m += $n, $n--){
							$str_clean[$i] .= sprintf("%i,",$race_gt_cnt_m[$i][$k][$m]);
						}					
					}
					for($m = 0; $m < $max_gt_cnt; $m++){
						$str_clean[$i] .= sprintf("%i,", $race_gt_cnt_f[$i][$k][$m]);
					}
					$str_clean[$i] =~ s/,$//;
					$str_clean[$i] .= ";";
					if($isx){
						$str_clean[$i] .= sprintf("%i,",$race_gt_cnt_con_m[$i][$k][0]);
						for($m = $allele_cnt, $n = $allele_cnt-1; $m < $max_gt_cnt; $m += $n, $n--){
							$str_clean[$i] .= sprintf("%i,",$race_gt_cnt_con_m[$i][$k][$m]);
						}					
					}
					for($m = 0; $m < $max_gt_cnt; $m++){
						$str_clean[$i] .= sprintf("%i,", $race_gt_cnt_con_f[$i][$k][$m]);
					}
					$str_clean[$i] =~ s/,$//;					
				}
				$str_clean[$i] .= "\t";
			}
			$str_hetz[$i] =~ s/\s+//g;																	# clean up some unwanted characters in output strings
			$str_clean[$i] =~ s/\s+$//;
			$str_hwe[$i] =~ s/\s+//g;
			$str_hetz[$i] =~ s/,$//;
			$str_hwe[$i] =~ s/,$//;
			
		} #for($i = 0; $i < $qc_grp_cnt; $i++)

		

		$vtype = get_vtype($ref, $var);
		
		print $dbg_handle "$bp\t$ref\t$var\t$vtype\t@vflag\n";
		for($i = 0; $i < @qc_names; $i++){
			$abh_str = "ABHet_" . $qc_names[$i] . "=";
			$str_abh[$i] = "";
			for($j = 0; $j < $allele_cnt; $j++){
				$temp = ".";
				$den = $abhet_cnt[$i][$j][0] + $abhet_cnt[$i][$j][1];
				if($den > 0){
					$temp = $abhet_cnt[$i][$j][0] / $den;
				}
				$abh_str .= "$temp" . ",";
				$str_abh[$i]  .= "$temp" . ",";
			}
			$abh_str =~ s/,$//;
			$str_abh[$i] =~ s/,$//;
			$info = $abh_str . ";" . $info;			
		}
		for($i = 0; $i < @qc_names; $i++){
			$temp = "VFLAGS_" . $qc_names[$i] . "=" . $vflag[$i] . ";";
			$info = $temp . $info;
		}

		$info =~ s/\s+$//;
=begin
		$info .= ";VariantInTargetFraction=" . "$ontarget_cnt" . "/" . "$sample_cnt";
		$temp = -1;
		if($sample_cnt > 0){
			$temp = $ontarget_cnt / $sample_cnt;
		}
		$info .= ";VariantInTargetRatio=" . "$temp";
=cut
		if($printflg){																					# this prints the flagged.vcf
			$head =~ s/\s+$//;																			# make sure we don't have any unwanted white space
			$info =~ s/^\s+//;
			$info =~ s/\s+$//;
			$fields =~ s/^\s+//;
			$fields =~ s/\s+$//;
			$flagged_line =~ s/^\s+//;
			$flagged_line =~ s/\s+$//;
			if($head =~ m/^chr/){
				print $flg "$head\t$info\t$fields\t$flagged_line\n";
			}
			else{
				print $flg "chr$head\t$info\t$fields\t$flagged_line\n";
			}
		}
		$head =~ m/^(\S+\s+){2,2}/;
		$pos = $&;
		$head =~ s/^\S+\s+\S+\s+//;
		if($agt == 1){
			if($normal){
				for($k = 0; $k < @genotypes; $k++){													# step through all genotypes
					if(($genotypes[$k] > 0) && $keep[$k]){
						$indivmaster[$k][2]++;														# singleton on ALT
						last;
					}
				}
			}
			else{
				for($k = 0; $k < @genotypes; $k++){
					if(($genotypes[$k] > 0) && ($genotypes[$k] < $allele_cnt) && $keep[$k]){
						$indivmaster[$k][2]++;														# singleton on REF
						last;
					}
				}
			}
		}
		elsif($agt == 2){
			$temp = 0;
			if($normal){
				for($k = 0; $k < @genotypes; $k++){													# step through all genotypes
					if($genotypes[$k] > 0 && $keep[$k]){
						if($genotypes[$k] >= $allele_cnt){
							$indivmaster[$k][3]++;													# private doubleton, 2 alt alleles
							$temp += 2;
						}
						else{
							$indivmaster[$k][4]++;													# doubleton, 1 alt allele
							$temp++;
						}
					}
					if($temp == 2){
						last;
					}				
				}
			}
			else{																					# REF is minor allele
				for($k = 0; $k < @genotypes; $k++){													# step through all genotypes
					if($genotypes[$k] < $allele_cnt && $keep[$k]){
						if($genotypes[$k] == 0){													# private doubleton, 2 ref alleles
							$indivmaster[$k][3]++;
							$temp += 2;
						}
						elsif($genotypes[$k] > 0){
							$indivmaster[$k][4]++;							
							$temp++;
						}
					}
					if($temp == 2){
						last;
					}	
				}	
			}
		}
		
		if($isx == 0){
			$mi_cnt = 0;
			$total_mi = 0; 
			$prop_mi = -1;
			for($k = 0; $k < @genotypes; $k++){														# mend. incons. and transfer of tempdata to permanent data				
				if(($genotypes[$k] >= 0) && $keep[$k]){							
					$fid = $par_index[$k][0];														# father index
					$mid = $par_index[$k][1];														# mother index
					$fgt = -1;
					$mgt = -1;
					if($fid >= 0){
						$fgt = $genotypes[$fid];
						if($fgt >= 0){
							$total_mi++;
						}
					}
					if($mid >= 0){
						$mgt = $genotypes[$mid];
						if($mgt >= 0){
							$total_mi++;
						}
					}
					if(($fgt >= 0) && ($mgt >= 0)){
						$indivmaster[$k][14]++;														# valid 2P MI-check
						if((($indivslave[$k][17] == $indivslave[$fid][17]) || ($indivslave[$k][17] == $indivslave[$fid][18])) && (($indivslave[$k][18] == $indivslave[$mid][17]) || ($indivslave[$k][18] == $indivslave[$mid][18]))){
							next;
						}
						if((($indivslave[$k][18] == $indivslave[$fid][17]) || ($indivslave[$k][18] == $indivslave[$fid][18])) && (($indivslave[$k][17] == $indivslave[$mid][17]) || ($indivslave[$k][17] == $indivslave[$mid][18]))){
							next;
						}							
						$indivmaster[$k][12]++;
						$mi_cnt++;
					}
					elsif($fgt >= 0){
						$indivmaster[$k][13]++;														# valid 1P MI-check
						if(($indivslave[$k][17] != $indivslave[$fid][17]) && ($indivslave[$k][17] != $indivslave[$fid][18]) && ($indivslave[$k][18] != $indivslave[$fid][17]) && ($indivslave[$k][18] != $indivslave[$fid][18])){
							$indivmaster[$k][11]++;							# 1P MI
							$mi_cnt++;
						}
					}
					elsif($mgt >= 0){
						$indivmaster[$k][13]++;														# valid 1P MI-check
						if(($indivslave[$k][17] != $indivslave[$mid][17]) && ($indivslave[$k][17] != $indivslave[$mid][18]) && ($indivslave[$k][18] != $indivslave[$mid][17]) && ($indivslave[$k][18] != $indivslave[$mid][18])){
							$indivmaster[$k][11]++;													# 1P MI
							$mi_cnt++;
						}						
					}
				}
			}
			if($total_mi > 0){
				$prop_mi = $mi_cnt / $total_mi;
			}
		}
		else{
			$mi_cnt = ".";
		}
		for($i = 0; $i < @qc_names; $i++){
			$temp = $outstream[$i];
			$depth = 0;
			if($depth_cnt[$i] > 0){
				$depth = $depth_sum[$i] / $depth_cnt[$i];
			}
			$pass_str = "";
			$fail_str = "";
			$maf_str = "";
			$num = 0;
			$gtf = 0;
			$den = 0;
			if($isx){
				$den = $qc_gt_cnt_m[$i][$missing_gt_index];												# those are the ./.
				for($j = 0; $j < $max_gt_cnt; $j++){
					$pass_str .= sprintf("%i,",$qc_gt_cnt_m[$i][$j]);									# list of passing genotypes
					$num += $qc_gt_cnt_m[$i][$j];
					$index = $j + $max_gt_cnt;
					$fail_str .= sprintf("%i,",$qc_gt_cnt_m[$i][$index]);								# list of failed genotypes
					$gtf += $qc_gt_cnt_m[$i][$index]
				}
				$pass_str =~ s/\S$//;
				$fail_str =~ s/\S$//;
				$pass_str .= ";";
				$fail_str .= ";";
			}
			$den += $qc_gt_cnt_f[$i][$missing_gt_index];												# those are the ./.

			for($j = 0; $j < $max_gt_cnt; $j++){
				$pass_str .= sprintf("%i,",$qc_gt_cnt_f[$i][$j]);										# list of passing genotypes
				$num += $qc_gt_cnt_f[$i][$j];
				$index = $j + $max_gt_cnt;
				$fail_str .= sprintf("%i,",$qc_gt_cnt_f[$i][$index]);									# list of failed genotypes
				$gtf += $qc_gt_cnt_f[$i][$index]
			}
			$den += $num + $gtf;

			$pass_str =~ s/\S$//;
			$fail_str =~ s/\S$//;
	
	
			for($j = 0; $j < $allele_cnt; $j++){
				$maf_str .= sprintf("%8.5f,", $mafp[$i][$j]);											# build MAF string for summary file
			}
			$maf_str =~ s/\S$//;
			$maf_str =~ s/\s+//g;
			$hidp = 0;
			if($vflag[$i] =~ m/5/){
				$hidp = 1;
			}
			$fout = 0;
			if($vflag[$i] ne "0"){
				$fout = 1;
			}
			$clean_str = "0";
			
			$gpass = 1;
			if(($vflag[$i] =~ m/1/) && ($vflag[$i] !~ m/11/)){											# all those redundant indicators directly derived from VFLAG
				$gpass = 0;
			}
			$mono = 0;
			if($vflag[$i] =~ m/3/){
				$mono = 1;
			}
			$crate = 0;
			if($den > 0){
				$crate = $num / $den;
			}
			$cb = 0;
			if($vflag[$i] =~ m/4/){
				$cb = 1;
			}
			$itr = 1;
			if($vflag[$i] =~ m/11/){
				$itr = 0;
			}
			$miss = $qc_gt_cnt_m[$i][$missing_gt_index] + $qc_gt_cnt_f[$i][$missing_gt_index];
			if($isx){
				printf($temp "%s\t%i\t%s\t%s\t%i\t%i\t%s\t%i\t%7.5f\t%i\t%i\t%s\t%7.5f\t%i\t%s\t.\t.\t.\t.\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n",$chr, $bp, $pass_str, $fail_str, $miss, $gtf, $clean_str, $mono, $crate, $cb, $gpass, $maf_str, $depth, $hidp, $str_abh[$i], $fout, $vflag[$i], $rsid, $ref, $var, $qual, $tranche, $vtype, $itr, $malehets[$i], $str_clean[$i], $str_hetz[$i], $str_hwe[$i]);
			}
			else{
				printf($temp "%s\t%i\t%s\t%s\t%i\t%i\t%s\t%i\t%7.5f\t%i\t%i\t%s\t%7.5f\t%i\t%s\t%i\t%i\t%7.5f\t.\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\n",$chr, $bp, $pass_str, $fail_str, $miss, $gtf, $clean_str, $mono, $crate, $cb, $gpass, $maf_str, $depth, $hidp, $str_abh[$i], $mi_cnt, $total_mi, $prop_mi, $fout, $vflag[$i], $rsid, $ref, $var, $qual, $tranche, $vtype, $itr, $str_clean[$i], $str_hetz[$i], $str_hwe[$i]);			
			}
		} # for($i = 0; $i < @qc_names; $i++)

		$k = $titv{"$ref$var"};
		if($k){
			$k += 0;																					# $j is either 5 or 6, the correct index for the sample counting array
			for($i = 0; $i < @genotypes; $i++){
				if($genotypes[$i] > 0){
					$indivmaster[$i][$k]++;
				}
			}
		}
		$line_cnt++;
		if(0 == ($line_cnt % 50)){
			print STDERR ".";
			if(0 == ($line_cnt % 5000)){
				print STDERR "\tbp= $bp\tcnt= $line_cnt\n";
			}
		}
	} # while(chomp($line = <$stream>))
	for($i = 0; $i < @qc_names; $i++){
		$temp = $outstream[$i];
		close($temp);
		if($debug){
			$temp = $dbg_out[$i];
			close($temp);
		}
	}
	if($debug){
		close($dbg_main);
	}
	if($printflg){
		close($flg);
	}
	close($stream);
	close(coh);
	print_indiv_master(\%parameters, \@indivmaster, $const_max_gt, \@keep, \@real_sex);						
	print "\n";
	close(dbg);
	if($chr23){
		close($chr23);
	}
=begin
	if($use_phwex){
		$R->stopR();
	}
=cut
}

##################################
# fix header in flagged VCF with stuff that's added to INFO column

sub fix_header
{
	my($ref) = $_[0];
	my(@qc_names) = @$ref;
	my($i);
	my($str) = "";
	for($i = 0; $i < @qc_names; $i++){
		$str .= "##INFO=<ID=";
		$str .= "VFLAG_" . $qc_names[$i];
		$str .= ",Number=.,Type=Integer,Description=\"Pipeline-specific QC variant flags\">";
		$str .= "\n";
	}
	for($i = 0; $i < @qc_names; $i++){
		$str .= "##INFO=<ID=";
		$str .= "ABHet_" . $qc_names[$i];
		$str .= ",Number=1,Type=Float,Description=\"Allelic Read Ratio\">";
		$str .= "\n";
	}
#	$str .= "##INFO=<ID=VariantInTargetFraction,Number=1,Type=String,Description=\"(on target) / total\">\n";
#	$str .= "##INFO=<ID=VariantInTargetRatio,Number=1,Type=Float,Description=\"(on target) / total as float\">\n";
	return($str);
}

##################################
# (\@mafp, \@vflag, 1+$#qc_names, $allele_cnt, $ref, $var)
sub get_vtype
{
	my($ref, $var) = @_;
	my($vtype) = "";
	
	my($lr) = length($ref);
	my($lv) = length($var);
	my($i, $l, $min, $max);
	my(@vars) = split(/,/, $var);
	if(($lr == 1) && ($lv == 1) && ($var ne "*")){
		$vtype = "SNV";
	}
	elsif($var =~ m/,/){
		$vtype = "Mult";
	}
	else{
		$vtype = "indel";
	}
	return($vtype);
}

##################################
# is_mono($lstgt[0], $lstgt[1], $lstgt[2])

sub is_mono{
	my($rr, $rv, $vv) = @_;
	my($result) = 1;
	if($rv > 0){
		$result = 0;
	}
	elsif(($rr > 0) && ($vv > 0)){
		$result = 0;
	}
	return($result);
}

##################################

#					$indivmaster[$j][0] = $sample_ids[$j];								# name
#					$indivmaster[$j][1] = 0;											# missing genotypes
#					$indivmaster[$j][2] = 0;											# singleton
#					$indivmaster[$j][3] = 0;											# private doubleton
#					$indivmaster[$j][4] = 0;											# doubleton
#					$indivmaster[$j][5] = 0;											# Ti
#					$indivmaster[$j][6] = 0;											# Tv
#					$indivmaster[$j][7] = 0;											# het. genotypes
#					$indivmaster[$j][8] = 0;											# hom. genotypes
#					$indivmaster[$j][9] = 0;											# sum depth of "PASS" genotypes
#					$indivmaster[$j][10] = 0;											# parent count
#					$indivmaster[$j][11] = 0;											# 1P MI count
#					$indivmaster[$j][12] = 0;											# 2P MI count
#					$indivmaster[$j][13] = 0;											# valid 1P counts
#					$indivmaster[$j][14] = 0;											# valid 2P counts
#					$indivmaster[$j][15] = 0;											# indel count
#					$indivmaster[$j][16] = 0;											# valdi gt cnt.
#					$indivmaster[$j][17] = 0;											# pre-QC 0/0
#					$indivmaster[$j][18] = 0;											# pre-QC 0/alt
#					$indivmaster[$j][19] = 0;											# pre-QC alt/alt
#					$indivmaster[$j][20] = 0;											# pre-QC 0/0 set ./.
#					$indivmaster[$j][21] = 0;											# pre-QC 0/alt set ./.
#					$indivmaster[$j][22] = 0;											# pre-QC alt/alt set ./.
#
#					failing genotypes deconvolut: 1-*($gt +1)
# Not needed for this version, but usre to be revived in the future
# print_indiv_master(\%parameters, \@indivmaster, $const_max_gt, \@keep, \@sex, \@indivfailcnt, \@indivpasscnt);
sub print_indiv_master
{
	my($ref1) = $_[0];
	my($ref2) =$_[1];
	my($const_max_gt) = $_[2];
	my($ref3) = $_[3];
	my($ref4) = $_[4];
	my(@keep) = @$ref3;
	my(@sex) = @$ref4;
	my(%parameters) = %$ref1;
	my(@indivmaster) = @$ref2;

	my(@gt_lst);
	my($i);
	my($temp);
	my($sum);
	my($total_gt);
	my($ratio);
	my($outfile) = $parameters{"output_prefix"};
	my($missing);
	my($titv);
	my($t1, $t2, $t3, $t4);
	my($passstr, $failstr);
	$outfile .= "_" . $parameters{"chr"} . "_";
	$outfile .= $parameters{"indiv_summary"};
	open(out,">$outfile");
	print out "SampleID\tSex\tPass\tFail\tMissing\tSingleton\tPrivate_Doubleton\tDoubleton\tHetHom\tIndMeanDepth\t1P_MI\t2P_MI\tMI_pairs\tTi\tTv\tTiTvRatio\n";
	for($i = 0; $i < @keep; $i++){
		if($keep[$i] == 0){
			next;
		}
		$passstr = "";
		$failstr = "";
		$total_gt = $indivmaster[$i][10];
		$passstr = $indivmaster[$i][17] . "," . $indivmaster[$i][18] . "," . $indivmaster[$i][19];
		$failstr = $indivmaster[$i][20] . "," . $indivmaster[$i][21] . "," . $indivmaster[$i][22];

		$passstr =~ s/(0,)+$//;
		$failstr =~ s/(0,)+$//;
		$passstr =~ s/,$//;
		$failstr =~ s/,$//;		
		$titv = -1;
		if($indivmaster[$i][6] > 0){
			$titv = $indivmaster[$i][5] / $indivmaster[$i][6];
		}
		if($passstr eq ""){
			$passstr = "0";
		}
		if($failstr eq ""){
			$failstr = "0";
		}		
		$t1 = $sex[$i] + 1;
		$t2 = -1;
		if($indivmaster[$i][8] > 0){
			$t2 = $indivmaster[$i][7]/ $indivmaster[$i][8];
		}
		$t3 = -1;
		if($total_gt > 0){
			$t3 = $indivmaster[$i][9] / $total_gt;
		}
		$t4 = 0;
		$t4 = $indivmaster[$i][13] + $indivmaster[$i][14];
		printf(out "%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%7.5f\t%8.5f\t%i\t%i\t%i\t%i\t%i\t%7.5f\n", $indivmaster[$i][0], $t1, $passstr, $failstr, $indivmaster[$i][1], $indivmaster[$i][2], $indivmaster[$i][3],$indivmaster[$i][4], $t2, $t3, $indivmaster[$i][11], $indivmaster[$i][12], $t4, $indivmaster[$i][5], $indivmaster[$i][6], $titv);
	}
	close(out);
}

#####################################
# get an inventory of all bp from the VCF file. This is needed to determine the mapping. 
# The initial value in the HASH is 0 for all positions

sub get_vcf_positions
{
	my($infile) = $_[0];
	my($stream);
	my($line, $bp);
	my(%hash);
	my(@lst);
	my(@bplst);
	my(@endpoints);
	my($i) = 0;
	if($infile =~ m/gz$/){
		open($stream, "-|", "gzip -dc $infile") or die "Can not open infile $infile\n";
	}
	else{
		open($stream,"<$infile") or die "Can not open infile $infile\n";
	}
	while(chomp($line = <$stream>)){
		if($line =~ m/^#/){
			next;
		}
		$line =~ s/^\S+\s+\S+\s+\S+\s+\S+//;
		$line = $&;
		@lst = split(/\s+/, $line);
		$bp = $lst[1];
		$hash{$bp} = "0";
		$bp += 0;
		$bplst[$i] = $bp;
		$endpoints[$i] = $bp + length($lst[3]) - 1;
		$i++;
		if(0 == ($i % 5000)){
			print "SNV cnt = $i\n";
		}
	}
	close($stream);
	print STDERR "total position = $i\n";
	return(\%hash, \@bplst, \@endpoints, $i);
}

#####################################
# the old hetZ calculation routine. 
# Leave as it might get revived in the future

=begin
sub excess_het_gt_only
{
	my($hom1) = $_[0];
	my($het) = $_[1];
	my($hom2) = $_[2];

	
	my($hetexp, $hetobs, $chisq, $maf);													# excess heterozygot estimate
	my($temp);
	my($x);
	my($result) = -999999;																# the default value for het_Z-score in case it can't be calculated due to missing heterozygots
	if(($hom1 == 0) && ($hom2 == 0)){
		return(999999);
	}
	my($total) = 2 * ($hom1 + $het + $hom2);
	if($total > 0){
		$maf = ($het + 2 * $hom2) / $total;
		if($maf > 0.5){																	# MAF in [0, 0.5] ... not that it really matters
			$maf = 1 - $maf;
		}
	}
	else{
		return(888888);
	}
	if(($maf > 0) && ($het > 0)){														# the actual calculation
		$total = $hom1 + $het + $hom2;
		$hetexp = 2 * $maf * (1 - $maf);
		$hetobs = $het / $total;
		$temp = $hetobs - $hetexp;
		$chisq = ($total * $temp * $temp) / ($hetobs * (1 - $hetobs));
		$result = sqrt($chisq);
		if($temp < 0){
			$result *= -1;
		}
	}
	elsif(($maf > 0) && ($het == 0)){													# case with no hets and MAF > 0			
		$result = 888888;																# we keep SNVs with het_z = 888888
	}
	elsif($maf == 0){																	# minor allele was in Dutch or it gets thrown out b/c monomorphic
		$result = 888888;
	}

	return($result);
}
=cut

##################################
##################################
# reads the control file and returns various parameters specified in a hash:
# 
# pedfile = name of file with pedigree info in standard 6 column format
# vcf =  name of VCF input file
# maxDP = max individual call depth
# maxDP = maximum mean depth of SNV
# min_call_rate = minimum call rate. Rate has to be > min_call_rate
# minDP = minimum depth, applies only to 0/0 genotypes. Depth must be > minDP
# etc.

sub get_data
{
	my($infile) = $_[0];
	my(%hash);
	my($line);
	my($found);
	my($temp);
	my($i, $j);
	my($nl);
	my(@lst, @lst2);
	open(inf,"<$infile") or die "Can not open control file $infile\n";

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^vcf:/){
			$line =~ s/^\S+\s+//;
			if($line =~ m/#/){
				$line = $`;
			}
			@lst = split(/\s+/, $line);
			$temp = "";
			for($i = 0; $i < @lst; $i++){
				$temp .= $lst[$i] . " ";
			}
			$temp =~ s/\s+$//;
			$hash{"vcf"} = $temp;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find vcf: Need one or more VCF input files to continue. Bye!\n";
		die;
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^maxDP:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"maxDP"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find maxDP: default = 500\n";
		$hash{"maxDP"} = 500;
	}

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^min_call_rate:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"min_call_rate"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find min_call_rate: default = 0.8\n";
		$hash{"min_call_rate"} = 0.8;
	}

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^snv_summary:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"snv_summary"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find snv_summary: default = snv_summary.out\n";
		$hash{"snv_summary"} = "snv_summary.txt";
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^indiv_summary:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"indiv_summary"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find indiv_summary: default = indiv_summary.out\n";
		$hash{"indiv_summary"} = "indiv_summary.out";
	}
=begin
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^debug:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"debug"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find debug: default = 0 (no debug)\n";
		$hash{"debug"} = 0;
	}
=cut
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^printFLG:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"printFLG"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find printFLG: default = 1\n";
		$hash{"printFLG"} = 1;
	}

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^minDP_female:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"minDP_female"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find minDP_female: default = 10\n";
		$hash{"minDP_female"} = 10;
	}

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^minDP_male:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"minDP_male"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find minDP_male: default = 10\n";
		$hash{"minDP_male"} = 10;
	}
		
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^minGQ:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"minGQ"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find minGQ: default = 20\n";
		$hash{"minGQ"} = 20;
	}		

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^chr:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"chr"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find chr: Expecting same as 1st column in capture.bed files, something like: chr17\n";
	}	

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^start_bp:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"start_bp"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find start_bp: default = 0\n";
		$hash{"start_bp"} = 0;	
	}

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^end_bp:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"end_bp"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find start_bp: default = 999999999\n";
		$hash{"end_bp"} = 999999999;	
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^add2mapsegments:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"add2mapsegments"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find add2mapsegments: default = 7\n";
		$hash{"add2mapsegments"} = 7;
	}	

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^isX:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"isX"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find isX to determine chr 1-22 or X: default = 0 (chr 1-22)\n";
		$hash{"isX"} = 0;
	}		
	
	if($hash{"isX"}){
		seek(inf, 0, 0);
		$found = 0;
		while(chomp($line = <inf>)){
			$line =~ s/^\s+//;
			if($line =~ m/^PAR:/){
				$line =~ s/^\S+\s+//;
				@lst2 = split(/\s+/, $line);
				$j = 0;
				for($i = 0; $i < @lst2; $i++){
					if($lst2[$i] =~ m/^[0-9]+$/){
						$lst[$j] = $lst2[$i];
						$j++;
					}
				}
				@lst = sort{$a <=> $b}@lst;
				$nl = "";
				for($i = 0; $i < @lst; $i ++){
					$nl .= "$lst[$i]\t";
				}
				$nl =~ s/\s+$//;
				$hash{"PAR"} = $nl;
				$found = 1;
				last;
			}
		}
		if(1 != ($#lst % 2)){
			die "PAR: need an even number of bp for regions, like: 1000\t2000\t40000000\t41000000\n";
		}
		if($found == 0){
			print STDERR "Did not find PAR for chrX: default = \"0\t0\"\n";
			$hash{"PAR"} = "0\t0";
		}	
	}

	if($hash{"isX"}){
		seek(inf, 0, 0);
		$found = 0;
		while(chomp($line = <inf>)){
			$line =~ s/^\s+//;	
			if($line =~ m/^pHWEx:/){
				$line =~ s/^\S+\s+//;
				$line =~ s/^\S+//;
				$hash{"pHWEx"} = $&;
				$found = 1;
				last;			
			}
		}
		if($found == 0){
			$hash{"pHWEx"} = 0;
		}
	}
	else{
		$hash{"pHWEx"} = 0;
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^minStat:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"minStat"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find minStat to calculate pop. stats: default = 6\n";
		$hash{"minStat"} = 6;
	}		
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^id_file:/){
			$line =~ s/^\S+\s+//;
			if($line =~ m/#/){
				$line = $`;
			}
			@lst = split(/\s+/, $line);
			$temp = "";
			for($i = 0; $i < @lst; $i++){
				$temp .= $lst[$i] . " ";
			}
			$temp =~ s/\s+$//;
			$hash{"id_file"} = $temp;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find id_file. Bye!\n";
		die;
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^capture_relations:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"capture_relations"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find name of capture_relations info file. Bye!\n";
		die;
	}	

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^usegz:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"usegz"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find usegz. Default= 0\n";
		$hash{"usegz"} = 0;
	}	
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^output_prefix:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$temp = $&;
			if((($temp =~ m/^\//) || ($temp =~ m/^\./)) && ($temp !~ m/\/$/)){
				$temp .= "/";
			}
			$hash{"output_prefix"} = $temp;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find name of output_prefix info file. Default is empty string.\n";
		$hash{"output_prefix"} = "";
	}	
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^id_col:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"id_col"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find id_col: default = 0\n";
		$hash{"id_col"} = 0;
	}

	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^race_subset_col:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"race_subset_col"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find race_subset_col: default = 6\n";
		$hash{"race_subset_col"} = 6;
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^error_rate:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"error_rate"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find error_rate: default = 0.0001\n";
		$hash{"error_rate"} = "0.0001";
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^error_threshold:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"error_threshold"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find error_threshold: default = 0.0001\n";
		$hash{"error_threshold"} = "0.0001";
	}
	
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^read_capture:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"read_capture"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find read_capture: default = 0\n";
		$hash{"read_capture"} = 0;
	}

	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^qc_subset_col:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"qc_subset_col"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not find qc_subset_col: default = 5\n";
		$hash{"qc_subset_col"} = 5;
	}
	
	seek(inf, 0, 0);
	$found = 0;
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		if($line =~ m/^maxTranche:/){
			$line =~ s/^\S+\s+//;
			$line =~ s/^\S+//;
			$hash{"maxTranche"} = $&;
			$found = 1;
			last;
		}
	}
	if($found == 0){
		print STDERR "Did not maxTranche: default = 99.7\n";
		$hash{"maxTranche"} = "99.7";
	}
	

	
	return(\%hash);
}

##################################
# print to screen the parameters in use

sub display_parameters
{
	my($ref) = $_[0];
	my(%hash) = %$ref;
	my($key, $value);
	my($i) = 0;
	my(@unsorted);
	my(@sorted);
	while(($key, $value) = each(%hash)){												# get all key-value pairs
		$unsorted[$i][0] = $key;
		$unsorted[$i][1] = $value;
		$i++;
	}
	@sorted = sort{$a->[0] cmp $b->[0]} @unsorted;										# sort by key
#	open(out,">parameters.txt");
	
	for($i = 0; $i < @sorted; $i++){													# print
		printf("# %20s = %30s\n", $sorted[$i][0], $sorted[$i][1]);
		printf($dbg_handle "# %20s = %30s\n", $sorted[$i][0], $sorted[$i][1]);
	}
	print $dbg_handle "\n";
#	close(out);
}

##################################
# beat the current FAM file into what is described under read_pedfile()
# simply permute the relevant columns into the needed position
sub fix_infile{
	my($infile) = $_[0];
	my($line);
	my(@lst);
	
	open(inf,"<$infile") or die "Can not open $infile\n";
	open(out,">temp_input.fam");
	print out "SampleID	FamID\tFatherID\tMotherID\tSex\tQC-subset\tRace\tuseSample\taff.stat\tuseHWE\tcaptureMap\n";	# this is the column order we want
	$line = <inf>;
	while(chomp($line = <inf>)){
		@lst = split(/\s+/, $line);
		print out "$lst[0]\t$lst[1]\t$lst[2]\t$lst[3]\t$lst[4]\t$lst[6]\t$lst[7]\t$lst[8]\t$lst[5]\t$lst[9]\t$lst[10]\t\n";
	}
	close(inf);
	close(out);
}

##################################
# expectes .fam like this, with header line:
# Since the input historically has the tendency to change it calls fix_infile() first.
# This function needs to be adjusted to deal with the FAM file format du jour.
# The result needs to conform exactly to what's shown below.
# 
#	0.	sample_ID (as in VCF)
#	1.	fam_ID (not used)
#	2.	father_ID (as in VCF)
#	3.	mother_ID (as in VCF)
#	4.	sex
#	5.	QC-subset
#	6.	race
#	7.	use sample (=1) or ignore sample (=0)
#	8.	aff.stat (0 = unknown/other, 1 = unaff, 2 = aff)
#	9.	use pHWE (0 = not, 1 = do use)
#  10.  capture map
#
#	Note that the sampleID MUST be unique.

sub read_pedfile
{
	my($infile) = $_[0];																# pedfile
	my($vcffile) = $_[1];
	my($ref) = $_[2];
	my(%parameters) = %$ref;
	my($id_col) = $parameters{"id_col"};												# column number (0,1,2,...) that holds ID we find in the VCF header line
	my($race_subset_col) = $parameters{"race_subset_col"};								# RACE subset-name in this column
	my($qc_subset_col) = $parameters{"qc_subset_col"};									# QC subset-name in this column
	my($isx) = $parameters{"isX"};														# is chrX (=1) or chr 1-22 (=0)
	my(@lst);																			# temp array to split line on white space
	my($line);																			# each line from $infile
	my(%vcf_ids1);																		# hash that holds IDs in VCF
	my(%parent_ids);																	# the parent IDs don't conform to the normal ID patterns
	my(%subj_sample);																	# key= subject_ID, value= sample_ID
	my(%race_index);																	# key = qc_group, value = index
	my($info);
	my($temp);
	my($fatalerror) = 0;																# counts IDs in VCF that are not listed in .fam file
	my($key, $value);
	my($p1, $p2);
	my(@sampleids);																		# 2-dim array [][0] = SampleID, [][1] = keep (0= delete, 1= keep)
	my(@race_groups);																	# DutchIsolate, Hispanic, etc.
	my(@qc_groups);																		# FAM, ADNI, CC
	my(@male_cnt);																		# number of males per QC-group. Used to calculate male-het cutoffs
	my(@parent_index);																	# index of "father, mother"
	my(@race_names);																	# actual names in orders as numeric substitudes
	my(@qc_names);																		# actual names in orders as numeric substitudes
	my(@sex);																			# -1 = unknown, 0 = male, 1 = female, set to all femal for chr1-22
	my(@real_sex);																		# -1 = unknown, 0 = male, 1 = female, stays this way for chr1-22
	my(@keep);																			# 0 or 1 if to keep this sample
	my(@gr_lst);																		# QC group listing
	my(@affstat);
	my(@usehwe);																		# 0 or 1 if used for zHet & HWE calculations
	my(@dummy);
	my(%capture);																		# capture kit and map used
	my(%race_gr_hash);																	# to numerically sort in the various RACEs like: DutchIsolate, Hispanic, etc.
	my(%qc_gr_hash);																	# to numerically sort in the various QC_GROUPs like: ADNI, CC, FAM 
	my(%race_assign_hash, %qc_assign_hash);																	
	my($id);																			# ID as read from $id_col in pedfile
	my($i, $j);
	my($stream);																		# read VCF
	my($sample_cnt);																	# valid samples, those with "1" in "KEEP" column
	my($qc_group_cnt);																	# Number of distinct QC groups
	my($race_group_cnt);																# Number of distinct RACE groups
	my($qc_group_cnt);
	my($size);																			# total IDs in ped/VCF
	my($key, $value);
	my($pattern) = qr/-/;																# regex pattern

	if($vcffile =~ m/gz$/){
		open($stream, "-|", "gzip -dc $vcffile") or die "Can not open infile $vcffile\n";
	}
	else{
		open($stream,"<$vcffile") or die "Can not open infile $vcffile\t$!\n";
	}
	while(chomp($line = <$stream>)){													# open VCF to find out which IDs are there
		if($line =~ m/^#/){
			if($line =~ m/^#CHROM/){
				$line =~ s/^(\S+\s+){9,9}//;
				@sampleids = split(/\s+/, $line);
				for($i = 0; $i < @sampleids; $i++){
					$vcf_ids1{$sampleids[$i]} = $i;										# enter each ID into hash
					if($sampleids[$i] =~ m/(?:.*?(-)){3}/){								# this is a mess b/c sometimes individual IDs and sample IDs were used which are different
						$temp = $&;														# catures everything BEFORE the 3rd "-"
						$temp =~ s/-$//;
						$parent_ids{$temp} = $i;
					}
					
				}
				last;
			}
		}
	}
	close($stream);
	fix_infile($infile);																# this function depends on the .fam file du jour
																						# the ever changing columns and column names of the .fam file made it necessary to have a function that permutes
																						# the columns so that the code below can process it. If the input.fam changes this function (fix_infile) needs
																						# to change as well.
	open(inf,"<temp_input.fam") or die "Can not open pedfile $infile.\n";
	$i = 0;
	$sample_cnt = 0;
	$line = <inf>;																		# header line
	while(chomp($line = <inf>)){
		$line =~ s/^\s+//;
		@lst = split(/\s+/, $line);
		$race_groups[$i] = -1;
		$qc_groups[$i] = -1;
		if($lst[7] eq "0"){																# do not use this sample/column in VCF, skip over it
			$keep[$i] = 0;
			$sex[$i] = -1;
			$real_sex[$i] = -1;
			$race_groups[$i] = -1;
			$affstat[$i] = -1;
			$usehwe[$i] = 0;
			$parent_index[$i][0] = -1;
			$parent_index[$i][1] = -1;
			$i++;
			next;
		}
		if($sampleids[$i] ne $lst[$id_col]){											# Need IDs in the same order in .fam and VCF
			print STDERR "Naming discrepancy between VCF ($sampleids[$i]) and pedfile ($lst[$id_col]). Input line $i\n";
			die;
		}
		if($lst[4] =~ m/\d/){
			$sex[$i] = $lst[4];
			$real_sex[$i] = $lst[4];
			$sex[$i] += 0;
			$real_sex[$i] += 0;
		}
		else{
			$sex[$i] = -1;
			$real_sex[$i] = -1;	
			if($isx){
				$keep[$i] = 0;
				$usehwe[$i] = 0;
			}
		}
		$affstat[$i] = $lst[8];
		if(!(($affstat[$i] == 0) || ($affstat[$i] == 1) || ($affstat[$i] == 2))){		# 0 = unknown; 1 = unaffected; 2 = affected
			print STDERR "Error: $lst[0] has affection status of $affstat[$i]\n";
		}
		$capture{$lst[5]} = $lst[10];
		$affstat[$i] += 0;
		$usehwe[$i] = $lst[9];
		$usehwe[$i] += 0;

		$temp = $race_gr_hash{$lst[$race_subset_col]};
		if($temp){
			$temp++;
			$race_gr_hash{$lst[$race_subset_col]} = $temp;
		}
		else{
			 $race_gr_hash{$lst[$race_subset_col]} = 1;
		}
		$temp = $qc_gr_hash{$lst[$qc_subset_col]};
		if($temp){
			$temp++;
			$qc_gr_hash{$lst[$qc_subset_col]} = $temp;
		}
		else{
			 $qc_gr_hash{$lst[$qc_subset_col]} = 1;
		}		
		$race_groups[$i] = $lst[$race_subset_col];
		$qc_groups[$i] = $lst[$qc_subset_col];
		$p1 = $vcf_ids1{$lst[2]};														# fatherID index
		$p2 = $vcf_ids1{$lst[3]};														# motherID index
		if($p1 eq ""){
			$p1 = -1;																	# NA
		}
		else{
			$p1 += 0;																	# convert to number
		}
		if($p2 eq ""){
			$p2 = -1;																	# NA
		}
		else{
			$p2 += 0;																	# convert to number
		}
		$parent_index[$i][0] = $p1;
		$parent_index[$i][1] = $p2;	
		$keep[$i] = 1;
		$i++;
		$sample_cnt++;
	}		
	close(inf);
	$size = $i;
	$i = 0;
	while(($key, $value) = each(%race_gr_hash)){
		$gr_lst[$i][0] = $key;
		$gr_lst[$i][1] = $value;			
		$i++;
	}
	@gr_lst = sort{$b->[1] <=> $a->[1]} @gr_lst;
	$j = $i;
	for($i = 0; $i < $j; $i++){
		$race_names[$i] = $gr_lst[$i][0];
	}
	$race_group_cnt  = $i;
	for($i = 0; $i < $race_group_cnt; $i++){
		$race_assign_hash{$gr_lst[$i][0]} = $i;
	}
	undef(@gr_lst);
	$i = 0;
	while(($key, $value) = each(%qc_gr_hash)){
		$gr_lst[$i][0] = $key;
		$gr_lst[$i][1] = $value;
		$i++;
	}
	@gr_lst = sort{$b->[1] <=> $a->[1]} @gr_lst;
	$qc_group_cnt  = $i;
	for($i = 0; $i < $qc_group_cnt; $i++){
		$qc_assign_hash{$gr_lst[$i][0]} = $i;
		$qc_names[$i] = $gr_lst[$i][0];
	}
	for($i = 0; $i < $qc_group_cnt; $i++){
		$male_cnt[$i] = 0;
	}
	for($i = 0; $i < $size; $i++){
		if($keep[$i]){
			$temp = $race_assign_hash{$race_groups[$i]};
			if($temp ne ""){
				$temp += 0;
				$race_groups[$i] = $temp;
			}
			else{
				$race_groups[$i] = -1;
			}
			$temp = $qc_assign_hash{$qc_groups[$i]};
			if($temp ne ""){
				$temp += 0;
				$qc_groups[$i] = $temp;
				if($sex[$i] == 0){
					$male_cnt[$temp]++;
				}
			}
			else{
				$qc_groups[$i] = -1;
			}
			$p1 = $parent_index[$i][0];															# check if parents are "problematic" and need to be excluded
			$p2 = $parent_index[$i][1];
			if(($p1 >= 0) && ($keep[$p1] == 0)){
				$parent_index[$i][0] = -1;
			}
			if(($p2 >= 0) && ($keep[$p2] == 0)){
				$parent_index[$i][1] = -1;
			}			
		}
		else{
			next;
		}
		if(($isx == 0) && ($sex[$i] <= 0)){														# chr 1-22, set all males and "NA" to female
			$sex[$i] = 1;
		}
	}
	undef(%parent_ids);
	undef(%race_gr_hash);
	undef(%qc_gr_hash);
	undef(@gr_lst);
	$ref = find_male_het_thresholds(\@male_cnt, $parameters{"error_threshold"}, $parameters{"error_rate"});
	my(@male_thresholds) = @$ref;
	print $dbg_handle "\tQC-subset\tmale_cnt\tmale_het_threshold(less_or_equal)\n";
	for($i = 0; $i < @qc_names; $i++){
		printf($dbg_handle "%80s\t%6i\t%3i\n",$qc_names[$i], $male_cnt[$i], $male_thresholds[$i]);
	}

#open(out,">qc_groups_sample.txt");
#for($i = 0; $i < @sampleids; $i++){
#	print out "$sampleids[$i]\t$qc_groups[$i]\tfa=$parent_index[$i][0]\tmo=$parent_index[$i][1]\n";
#}
#close(out);
	return(\@sex, \@parent_index, \@keep, \@race_groups, \@sampleids, \%vcf_ids1, \%race_assign_hash, \@qc_groups, \%qc_assign_hash, \@race_names, \@qc_names, \@affstat, $race_group_cnt, $qc_group_cnt, \@usehwe, \%capture, \@male_thresholds, \@real_sex);
}

##################################
=begin
N=c(16,27,39,72,99,190,259,276,302,902,1225,1979,2552) #respective sample sizes (known)
e=0.0001 #error rate
t=0.0001 #Prob of false positive
qbinom((1-t),N,e) #gives number c such that P(bin>c)<t. We reject a position if #male hets > c (strictly greater than)


#some examples
#e=0.001,t=0.001
# c=1 1 1 2 2 2 3 3 3 5 6 8 9
#e=0.001,t=0.0001 : lower t means higher c (i.e., throw out fewer positions)
# 2  2  2  2  3  3  4  4  4  6  7  9 10
#e=0.0001,t=0.001 : lower e means lower c (i.e., throw out more positions)
# 1 1 1 1 1 1 1 1 1 2 2 3 3
#e=0.0001,t=0.0001 : lower e means lower c (i.e., throw out more positions)
# 1 1 1 1 1 2 2 2 2 3 3 3 4

#Note that wth N=5000, e=0.0001, t=0.0001 we get c=5 which is what we used previously in the Discovery with ~5000 males

=cut
sub find_male_het_thresholds
{
	my($ref) = $_[0];
	my($error_threshold) = $_[1];
	my($error_rate) = $_[2];
	$error_threshold += 0;
	$error_rate += 0;
	my(@male_cnt) = @$ref;
	my($i, $j, $temp);
	my(@thresholds);
	for($i = 0; $i < @male_cnt; $i++){
		$j = 0;
#		print "### $male_cnt[$i] ###\n";
		while($error_threshold < cumP_lg_x($male_cnt[$i], $j, $error_rate)){
			$j++;
		}
		$thresholds[$i] = $j;
	}
	return(\@thresholds);
}

##################################
# some sort of binary search on regions

sub is_in_region
{
	my($chr) = $_[0];
	my($bp) = $_[1];
	my($ref) = $_[2];
	my($size) = $_[3];
	my($zeros) = "";
	my($l) = length("$bp");
	my($i);
	my($result) = 0;
	my($item, $middle);
	my(@list) = @$ref;
	my($lower) = 0;
	my($upper) = $size - 2;
	my($diff, $estimate);
	for($i = $l; $i < 9; $i++){
		$zeros .= "0";
	}
	$item = "$chr" . $zeros . "$bp";
	$item += 0;
	$middle = -1;

	while($lower <= $upper){
		if($middle == -1){
			$diff = $list[$upper][1] - $list[1][0];
			$estimate = ($item - $list[1][0]) / $diff;
			$middle = int($estimate * $size);
			if($middle < 1){
				$middle = 1;
			}
			elsif($middle >= $upper){
				$middle = $upper - 1;
			}
		}
		else{
			$middle = int(($lower + $upper) / 2);
		}
		if(($item >= $list[$middle][0]) && ($item < $list[$middle+1][0])){
			if($item <= $list[$middle][1]){
				$result = $middle;
			}
			last;
		}
		elsif($item < $list[$middle][0]){
			$upper = $middle - 1;
		}
		else{
			$lower = $middle + 1;
		}
	}
	return($result);
}

##################################

=begin
VFLAGS
0	Passed
1	Variant failed preliminary QC: GATK “FILTER” ≠ “PASS” or is in tranche ≥ 99.7%
2	Variant failed preliminary QC: All genotypes have DP<10 and/or GQ<20
3	Monomorphic
4	Call Rate ≤ 80%
5	Average mean depth > 500 reads
7	excessive male-hets
11	outside target region
			
=cut
# determines which VFLAGs are needed

sub get_vflag
{
	my($ref1, $depth, $ref2, $tranche, $altgt, $ontarget, $allelecnt, $malehets, $male_thres) = @_;
	my(@lst) = @$ref1;
	my(%param) = %$ref2;
	my($miss_thres) = 1 - $param{"min_call_rate"};
	my($depth_thres) = $param{"maxDP"};	
	my($maxtranche) = $param{"maxTranche"};
	my($vflag) = "";
	my($missing) = $lst[6];
	my($low_qual) = 0;									#$lst[3] + $lst[4] + $lst[5];
	my($good_qual) = 0;									#$lst[0] + $lst[1] + $lst[2];
	my($altgt) = 0;
	my($maf);
	my($sum);
	my($temp);
	my($passes) = 1;
	my($i);
	my($thres);
	my(%homo_hash);
	my($lone_gt) = -1;
	my($gtcnt) = $allelecnt * ($allelecnt+1) / 2;
	my($missing_gt_index) = 2 * $gtcnt;
	if($maxtranche ne "PASS"){
		$maxtranche += 0.0000001;
	}
	else{
		$maxtranche = -1;
	}
	
	for($i = 0; $i < $gtcnt; $i++){
		$good_qual += $lst[$i];
		if($lst[$i]){
			$altgt++;
			$lone_gt = $i;								# this ONLY matters is $altgt == 1. Then it records the only exisiting genotype
		}
	}
	for($i = $gtcnt; $i < (2*$gtcnt); $i++){
		$low_qual += $lst[$i];
	}
	for($i = 0; $i < $allelecnt; $i++){
		$temp = ($i * $allelecnt) + $i - (($i * ($i + 1)) / 2);
		$homo_hash{$temp} = 1;							# puts all homozygote GTs into a hash
	}
#print "$missing = ($lst[6] + $low_qual) / ($lst[6] + $low_qual + $good_qual)\n";
	$missing = ($lst[$missing_gt_index] + $low_qual) / ($lst[$missing_gt_index] + $low_qual + $good_qual);

	if($tranche ne "PASS"){
		$tranche =~ m/to/;
		$tranche = $';

		if($tranche > $maxtranche){
			$vflag .= "1,";
			$passes = 0;
		}
	}
	if($good_qual == 0){
		$vflag .= "2,";
		$passes = 0;		
	}
	if($altgt == 0){
		$vflag .= "3,";
		$passes = 0;	
	}
	elsif(($altgt == 1) && ($homo_hash{$lone_gt})){		# only one GT AND it is homozygote ==> monomorphic (VFLAG == 3)
		$vflag .= "3,";
		$passes = 0;
	}
	if($missing >= $miss_thres){
		$vflag .= "4,";
		$passes = 0;
	}
	if($depth > $depth_thres){
		$vflag .= "5,";
		$passes = 0;
	}
	if($malehets > $male_thres){
		$vflag .= "7,";
		$passes = 0;		
	}
#	if(($altgt == 1) && ($vflag !~ m/3,/)){				# only a single genotype, but NOT monomorphic
#		$vflag .= "9,";
#	}
	if(!$ontarget){
		$vflag .= "11,";
		$passes = 0;		
	}
	$vflag =~ s/,$//;
	if($vflag eq ""){
		$vflag = "0";
	}
	return($vflag, $passes);
}

##################################
# qc_group_cnt(\@keep, \@race_groups, \@qc_groups, \@race_names, \@qc_names, \@sex);
# so we will produces a 2-dim array [i][] = QC-groups, [][j] = races
# each entry will have the count of FEMALE race members per QC group. In the QC-group specific output only race with >= 1 member will be represented. 
# Those races with 0 (zero) FEMALE members will not be in the QC-group specific summary file

sub qc_group_cnt
{
	my($ref1, $ref2, $ref3, $ref4, $ref5, $ref6) = @_;
	my(@keep) = @$ref1;
	my(@race_groups) = @$ref2;
	my(@qc_groups) = @$ref3;
	my(@race_names) = @$ref4;
	my(@qc_names) = @$ref5;
	my(@sex) = @$ref6;												# -1 = unknown, 0 = male, 1 = female						
	my($i, $j);
	my(@qc_race_cnt);												# [i][] = QC-groups, [][j] = races
	my($race_id, $qc_id);
	for($i = 0; $i < @qc_names; $i ++){
		for($j = 0; $j < @race_names; $j++){
			$qc_race_cnt[$i][$j] = 0;
		}
	}
	for($i = 0; $i < @keep; $i ++){
		if($keep[$i]  && ($sex[$i] == 1)){
			$race_id = $race_groups[$i];
			$qc_id = $qc_groups[$i];
			$qc_race_cnt[$qc_id][$race_id]++;
		}
	}
	print "\n\t";
	for($j = 0; $j < @race_names; $j++){
		print $dbg_handle "\t$race_names[$j]";
	}
	print $dbg_handle "\nFemales only on X\n";
	for($i = 0; $i < @qc_names; $i++){
		print $dbg_handle "$qc_names[$i]";
		for($j = 0; $j < @race_names; $j++){
			print $dbg_handle "\t$qc_race_cnt[$i][$j]";
		}
		print $dbg_handle "\n";
	}	
	print $dbg_handle "\n";
	return(\@qc_race_cnt);
}

##################################
# snphwe.pl: A Perl implementation of the fast exact Hardy-Weinberg Equilibrium 
# test for SNPs as described in Wigginton, et al. (2005). 
#
# Copyright 2010 Joshua Randall
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Author:
#    This software was written Joshua C. Randall <jcrandall@alum.mit.edu>
#    and is a port to Perl of algorithms implemented by others in C and R.
#
# Attribution:
#    This software is based entirely on the C and R implementations of the 
#    algorithms for exact HWE tests as described in Wigginton, et al. (2005), 
#    which were originally written by Jan Wigginton and released into the 
#    public domain.  C, R, and Fortran implementations of these algorithms 
#    are available for download at:
#       http://www.sph.umich.edu/csg/abecasis/Exact/
#
# Citation: 
#    This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
#    Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
#    Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76(5): 887 - 893.
#    Please cite this work when using this code.
#
# Usage: 
#    This software is a Perl library, intended to be used within other programs.
#    To use this library directly from the command-line, you can run Perl with a 
#    one-line program such as:
#
#       perl -e 'require "snphwe.pl"; print(snphwe(@ARGV))' 57 14 50
#
#    Where the three numbers at the end are the observed counts of the three 
#    genotypes: first the heterozygote count, then one of the homozygote genotype 
#    counts, and finally the other homozygote genotype count, in that order.  
#
#    The example above, which would be for 57 Aa, 14 aa, and 50 AA, should print 
#    the resulting P-value, which in this case is 0.842279756570793, to the 
#    standard output.
#
# Note:
#    Code for the alternate P-value calculation based on p_hi/p_lo that was 
#    included in the Wigginton, et al. C and R implementations (but was 
#    disabled) has been included here, but has not been tested.  It is 
#    therefore commented out.  If you wish to make use of this code, please 
#    verify it functions as desired.
#
sub snphwe {
    my($obs_hets, $obs_hom1, $obs_hom2) = @_;

    if($obs_hom1 < 0 || $obs_hom2 < 0 || $obs_hets < 0) {
		return(-1);
    }

    # rare homozygotes
    my $obs_homr;

    # common homozygotes
    my $obs_homc;
    if($obs_hom1 < $obs_hom2) {
		$obs_homr = $obs_hom1;
		$obs_homc = $obs_hom2;
    } else {
		$obs_homr = $obs_hom2;
		$obs_homc = $obs_hom1;
    }

    # number of rare allele copies
    my $rare_copies = 2 * $obs_homr + $obs_hets;

    # total number of genotypes
    my $genotypes = $obs_homr + $obs_homc + $obs_hets;

    if($genotypes <= 0) {
		return(-1);
    }
    
    # Initialize probability array
    my @het_probs;
    for(my $i=0; $i<=$rare_copies; $i++) {
		$het_probs[$i] = 0.0;
    }

    # start at midpoint
    my $mid = int($rare_copies * (2 * $genotypes - $rare_copies) / (2 * $genotypes));

    # check to ensure that midpoint and rare alleles have same parity
    if(($rare_copies & 1) ^ ($mid & 1)) {
		$mid++;
    }
    
    my $curr_hets = $mid;
    my $curr_homr = ($rare_copies - $mid) / 2;
    my $curr_homc = $genotypes - $curr_hets - $curr_homr;

    $het_probs[$mid] = 1.0;
    my $sum = $het_probs[$mid];
    for($curr_hets = $mid; $curr_hets > 1; $curr_hets -= 2) {
		$het_probs[$curr_hets - 2] = $het_probs[$curr_hets] * $curr_hets * ($curr_hets - 1.0) / (4.0 * ($curr_homr + 1.0) * ($curr_homc + 1.0));
		$sum += $het_probs[$curr_hets - 2];

	# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
		$curr_homr++;
		$curr_homc++;
    }

    $curr_hets = $mid;
    $curr_homr = ($rare_copies - $mid) / 2;
    $curr_homc = $genotypes - $curr_hets - $curr_homr;
    for($curr_hets = $mid; $curr_hets <= $rare_copies - 2; $curr_hets += 2) {
		$het_probs[$curr_hets + 2] = $het_probs[$curr_hets] * 4.0 * $curr_homr * $curr_homc / (($curr_hets + 2.0) * ($curr_hets + 1.0));
		$sum += $het_probs[$curr_hets + 2];
	
	# add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
		$curr_homr--;
		$curr_homc--;
    }

    for(my $i=0; $i<=$rare_copies; $i++) {
		$het_probs[$i] /= $sum;
    }

    # Initialise P-value 
    my $p_hwe = 0.0;

    # P-value calculation for p_hwe
    for(my $i = 0; $i <= $rare_copies; $i++) {
		if($het_probs[$i] > $het_probs[$obs_hets]) {
			next;
		}
		$p_hwe += $het_probs[$i];
    }
    
    if($p_hwe > 1) {
		$p_hwe = 1.0;
    }

    return($p_hwe);
}

##################################
# read the capture_info file to figure out the names of the file that contain the capture regions

sub get_capture_relations
{
	my($ref1, $ref2, $infile) = @_;
	my(%capture) = %$ref1;
	my(%qc_assign) = %$ref2;
	my($i, $j, $line, $temp);
	my(@lst);
	my(@bed_lst);
	my(%name2int);
	my(@capture2);
	my($key, $val);
	open(inf,"<$infile") or die "Can not open $infile\n";
	$i = 1;
	$j = 0;
	while(chomp($line = <inf>)){
		@lst = split(/\s+/, $line);
		$name2int{$lst[0]} = $i;
		$bed_lst[$j] = $lst[1];
		while(($key, $val) = each(%capture)){
			if($lst[0] eq $val){
				$temp = $qc_assign{$key};
				if($temp ne ""){
					$temp += 0;
					$capture2[$temp] = $i;
				}
				else{
					print STDERR "get_capture_relations: unknown qc_group of $key\n";
				}
			}
		}
		$j++;
		$i *= 2;
	}
	close(inf);

	return(\@capture2, \@bed_lst);
}

##################################

sub read_capture_maps
{
	my($ref1, $ref2, $ref3, $ref4, $chr) = @_;
	my(@bp_lst) = @$ref1;
	my(@endpoints) = @$ref2;
	my(%bp_hash) = %$ref3;
	my(@map_lst) = @$ref4;
	my(@ontarget);
	my($i, $j, $temp);
	my($p) = 1;															# $p will take on various powers of 2 = {1,2,4,8,16,...}
	for($i = 0; $i < @map_lst; $i++){
	
		$ref1 = cap_map($map_lst[$i], \@bp_lst, \@endpoints, $chr);		# process each target map for the various QC-subsets
		@ontarget = @$ref1;
		for($j = 0; $j < @bp_lst; $j++){
			if($ontarget[$j]){											# the ontarget positions (=1) get translated tot he appropriate 2^$i
				$temp = $bp_hash{$bp_lst[$j]};
				$temp += $p;
				$bp_hash{$bp_lst[$j]} = $temp;
			}
		}
		$p *= 2;														# each target map will have its unique power of 2, namely 2^$i
	}
	return(\%bp_hash);
}

##################################
# return a list of 0 |1 for each bp if it is on the target region or not.
sub cap_map
{
	my($infile, $ref1, $ref2, $chr) = @_;
	my(@bp_lst) = @$ref1;
	my(@endpoints) = @$ref2;
	my(@ontarget);
	my($i) = 0;
	my($j);
	my(@lst);
	my($bp, $endpt);
	my($mapsize, $bpsize, $line);
	my(@map);
	open(inf,"<$infile") or die "Can not open $infile\n";
	while(chomp($line = <inf>)){
		@lst = split(/\s+/, $line);
		if($lst[0] ne $chr){
			next;
		}
		$lst[1] += 0;
		$lst[2] += 0;
		$map[$i][0] = $lst[1];
		$map[$i][1] = $lst[2];
		$i++;
	}
	$map[$i][0] = 999999999;											# set end point of map
	$map[$i][1] = 999999999;
	$mapsize = $i + 1;													# remember map size
	$bpsize = 1 + $#bp_lst;												# @bp_lst is sorted so we step through it and see what's on target and what's off target 
	$i = 0;																# index into @map
	$j = 0;																# index into @bp_lst and @endpoints
	while($j < $bpsize){
		$bp = $bp_lst[$j];
		$endpt = $endpoints[$j];
		while(($bp > $map[$i][1]) && ($i < ($mapsize-1))){
			$i++
		}
		if($bp == $endpt){
			if($bp < $map[$i][0]){
				$ontarget[$j] = 0;
				$j++;
				next;
			}
			if(($bp >= $map[$i][0]) && ($bp <= $map[$i][1])){
				$ontarget[$j] = 1;
				$j++;
				next;
			}
		}
		else{															# need to test segment instead of single point
			if($endpt >= $map[$i][0]){
				$ontarget[$j] = 1;
				$j++;
				next;
			}
			elsif($endpt < $map[$i][0]){
				$ontarget[$j] = 0;
				$j++;
				next;
			}			
			
		}
	}
	return(\@ontarget);
}

##################################
#

sub cumP_lg_x($male_cnt[$i], $j, $error_rate)
{
	my($N, $x, $p) = @_;
	my($ref) = ln_fact($N);
	my(@ln_fac) = @$ref;
	my($result) = 0;
	my($nCi, $temp);
	for($i = 0; $i <= $x; $i++){
		$nCi = $ln_fac[$N] - $ln_fac[$i] - $ln_fac[$N - $i];
		$temp = $i * log($p) + ($N - $i) * log(1-$p) + $nCi;
		$result += exp($temp);
	}
	return(1 - $result);	
}

##################################

sub ln_fact
{
	my($n) = $_[0];
	my(@ln_fac);
	my($i);
	$ln_fac[0] = 0;
	for($i = 1; $i <= $n; $i++){
		$ln_fac[$i] = log($i) + $ln_fac[$i-1];
	}
	return(\@ln_fac);
}

##################################

sub find_start_bp
{
	my($stream, $start_bp, $file_size) = @_;
	my($current) = tell($stream);
	my($found) = 0;
	my($start) = tell($stream);
	my($min) = 0;
	my($end) = $file_size;
	my($previous, $line, $temp, $diff, $bp);
	my($i) = 0;
	seek($stream, 0, 0);
	while($line = <$stream>){
		if($line =~ m/^#/){
			$previous = tell($stream);
			next;
		}
		else{
			$line =~ m/^\S+\s+\S+/;
			$temp = $&;
			$temp =~ m/\S+$/;
			$bp = $&;		
			seek($stream, $previous, 0);
			$min = $previous;
			$start = $min;
			if($bp >= $start_bp){
				$found = 1;
			}
			last;
		}
	}
	while(($found == 0) && ($i < 18)){
		chomp($line = <$stream>);
		$previous = $current;
		$current = tell($stream);
#print STDERR "\nstart= $start current= $current end= $end\n";
		$line =~ m/^\S+\s+\S+/;
		$temp = $&;
		$temp =~ m/\S+$/;
		$bp = $&;
#print STDERR "1st: bp= $bp looking for $start_bp current= $current i= $i\n";
		if($bp == $start_bp){
			seek($stream, $previous, 0);								# move back to previous position, done
			$found = 1;
			last;
		}
		if($bp > $start_bp){
			$end = tell($stream);
			$temp = int(($current - $start) / 2);
			$temp *= -1;
#print STDERR "seek(, $temp, 1)\n";
			if(($current + $temp) < $min){
				seek($stream, $min, 0);
				$found = 1;
				last;
			}
			seek($stream, $temp, 1);
			$line = <$stream>;
			
		}
		elsif($bp < $start_bp){
			$start = $current;
			$temp = int(($end - $current) / 2);
#print STDERR "seek(, $temp, 1)\n";
			seek($stream, $temp, 1);
			$line = <$stream>;
		}
		$i++;
	}
	if($found == 0){
		$previous = tell($stream);
		chomp($line = <$stream>);
		$current = tell($stream);
		$line =~ m/^\S+\s+\S+/;
		$temp = $&;
		$temp =~ m/\S+$/;
		$bp = $&;
		$diff = $current - $previous;
		while(($bp > $start_bp) && (!$found)){
			$temp = -1 * $diff * 10;
			$current = tell($stream);
			if(($current + $temp) < $min){
				seek($stream, $min, 0);
				$found = 1;
				last;
			}
			seek($stream, $temp, 1);
			chomp($line = <$stream>);
			$previous = tell($stream);
			chomp($line = <$stream>);
			$line =~ m/^\S+\s+\S+/;
			$temp = $&;
			$temp =~ m/\S+$/;
			$bp = $&;	
#print STDERR "2nd: bp= $bp looking for $start_bp\n";
			if($bp == $start_bp){
				seek($stream, $previous, 0);
				$found = 1;
				last;
			}				
		}
	}
	while($found == 0){
		$previous = tell($stream);
		chomp($line = <$stream>);
		$line =~ m/^\S+\s+\S+/;
		$temp = $&;
		$temp =~ m/\S+$/;
		$bp = $&;
#print STDERR "3rd: bp= $bp looking for $start_bp\n";
		if($bp >= $start_bp){
			seek($stream, $previous, 0);
			$found = 1;
		}
		else{
			seek($stream, $min, 0);
			$found = 1;
			last;		
		}
	}
#print STDERR "before return bp= $bp found= $found\n";
	return($stream);
}

##################################
