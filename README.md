# qc_vcf
To run QC: 
make_20K_multi_x_flagged_vcf.pl  name_of_control_file

The typical control file looks like this:
id_file:			                  name_of_fam_file_here
id_col:               		      0
qc_subset_col:                  5
race_subset_col:                6
indiv_summary:              	  indiv_summary.out
maxDP:                  	      500
maxTranche:                	    99.7
minDP:              		        10
minGQ:                    	    20
start_bp:			                  0
end_bp:				                  999999999
min_call_rate:                	0.8
read_capture:			              0
usegz:				                  0
snv_summary:                	  snv_summary.out
output_prefix:                  unique_prefix_here			      # some prefix so it doesn't overwrite previous results
vcf:                		        name_of_vcf_file_here			    # your VCF here
add2mapsegments:                0					                    # add this many bp to each map segment
capture_relations:              capture_info.txt
isX:                            0					                    # change to 1 if is chrX
PAR:                            10001 2781479 155701383 156030895			# depending on the genetic map used these numbers might change
chr:                            chrX							# must correspond to 1st column in VCF, chr1, chr2, ..., chrX


The .fam file has the following format:
SampleID	FID	FA	MO	SEX	AFF	QC-Subset	Ethnicity	useSample	useHWE	CaptureMap
A-ACT-AC000014-BL-NCR-15AD78694	0	0	0	0	1	one_subgroup	NonHispanicWhite	1	1	all.bed
A-ACT-AC000034-BL-NCR-16AD84906	0	0	0	0	1	one_subgroup	NonHispanicWhite	1	1	all.bed
A-ACT-AC000057-BL-NCR-15AD78356	0	0	0	0	2	one_subgroup	NonHispanicWhite	1	0	all.bed
...

SampleID = Exactly like it is displayed in the VCF
FID needed for mendincons 
FA needed for mendincons 
MO needed for mendincons 
SEX only needed for chrX
AFF needed for HWE. 0 = unknown; 1 = unaffected; 2 = affected
QC-Subset indicates distinct set that was typed at specific center and/or technology
Ethnicity for HWE and HetZ calculations
useSample a 0 or 1 if this sample is to be included
useHWE a 0 or 1 if this sample is to be included in HWE calculation
CaptureMap needed for targeted maps, otherwise I can provide a “catch all” map


sven@HIHG-Ubuntu:~/checkout/hihg_um/qc_vcf$ docker run -v /hihg:/hihg -v /Volumes/Synology:/Volumes/Synology -v /home/sven/checkout/hihg_um/qc_vcf:/app  -it hihg-um/sven/perl make_20K_multi_x_flagged_thread_vcf.pl cuadi_1.ctrl
Did not find printFLG: default = 1
Did not find minDP_female: default = 10
Did not find minDP_male: default = 10
Did not find start_bp: default = 0
Did not find start_bp: default = 999999999
Did not find minStat to calculate pop. stats: default = 6
Did not find usegz. Default= 0
Did not find error_rate: default = 0.0001
Did not find error_threshold: default = 0.0001
Did not find read_capture: default = 0
#      add2mapsegments =                              0
#    capture_relations =               capture_info.txt
#                  chr =                           chr1
#               end_bp =                      999999999
#           error_rate =                         0.0001
#      error_threshold =                         0.0001
#               id_col =                              0
#              id_file = /Volumes/Synology/shared/qc_vcf/CUADI/cuadi.fam
#        indiv_summary =              indiv_summary.out
#                  isX =                              0
#                maxDP =                            500
#           maxTranche =                           99.7
#         minDP_female =                             10
#           minDP_male =                             10
#                minGQ =                             20
#              minStat =                              6
#        min_call_rate =                            0.8
#           output_dir =     /Volumes/Synology/results/
#        output_prefix =                          chr1_
#                pHWEx =                              0
#             printFLG =                              1
#        qc_subset_col =                              5
#      race_subset_col =                              6
#         read_capture =                              0
#          scratch_dir =     /Volumes/Synology/scratch/
#          snv_summary =                snv_summary.out
#             start_bp =                              0
#              threads =                             16
#                usegz =                              0
#                  vcf = /hihg/studies/AD/cohort_vcf/EOAD/vcf/EOAD-wgs.chr1.snp.indel.recalibrated.vcf.gz
#              version = QC-version 10/25/2022 2:08pm by Mike Schmidt, mschmidt@med.miami.edu

reading /Volumes/Synology/shared/qc_vcf/CUADI/cuadi.fam
gzip: /hihg/studies/AD/cohort_vcf/EOAD/vcf/EOAD-wgs.chr1.snp.indel.recalibrated.vcf.gz: Permission denied
Naming discrepancy between VCF () and pedfile (A-MIA-UM010501-BL-MIA-201819558). Input line 0
Died at make_20K_multi_x_flagged_thread_vcf.pl line 2704, <inf> line 2.