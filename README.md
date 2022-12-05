# qc_vcf
To run QC: 
make_20K_multi_x_flagged_vcf.pl  name_of_control_file

### to run docker

docker run -v /data:/data -v /scratch:/scratch -v `pwd`:/app  -it libbio-db-hts-perl:latest make_20K_multi_x_flagged_thread_vcf.pl <ctrl_file>

### control file

id_file:			/data/cuadi.fam
id_col:               		0
qc_subset_col:                  5
race_subset_col:                6
indiv_summary:              	indiv_summary.out
maxDP:                  	500
maxTranche:                	99.7
minDP:              		10
minGQ:                    	20
min_call_rate:                	0.8
snv_summary:                	snv_summary.out
output_prefix:                  chr21_
scratch_dir:			/data/scratch
output_dir:			/data/results
vcf:                		/data/vcf/CUADI-wgs.chr21.snp.indel.recalibrated.vcf.gz
add2mapsegments:                0
capture_relations:              capture_info.txt
isX:                            0
PAR:                            10001 2781479 155701383 156030895
chr:                            chr21
threads:			4


The typical control file looks like this:
id_file:		name_of_fam_file_here
id_col:               	0
qc_subset_col:          5
race_subset_col:        6
indiv_summary:          indiv_summary.out
maxDP:                  500
maxTranche:             99.7
minDP:              	10
minGQ:                  20
start_bp:		0
end_bp:			999999999
min_call_rate:          0.8
read_capture:		0
usegz:			0
snv_summary:            snv_summary.out
output_prefix:          unique_prefix_here # some prefix so it doesn't overwrite previous results
vcf:                	name_of_vcf_file_here  # your VCF here
add2mapsegments:        0 # add this many bp to each map segment
capture_relations:      capture_info.txt
isX:                    0 # change to 1 if is chrX
PAR:                    10001 2781479 155701383 156030895 # depending on the genetic map used these numbers might change
chr:                    chrX # must correspond to 1st column in VCF, chr1, chr2, ..., chrX


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


Assumptions:
The following directories MUST exist:
results
scratch	(exclusive use by single instance of the script)
