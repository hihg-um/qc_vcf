# !/bin/bash

docker run -v /hihg:/hihg -v /Volumes/Synology:/Volumes/Synology -v /home/sven/checkout/hihg_um/qc_vcf:/app  -it hihg-um/sven/perl scripts/make_20K_multi_x_flagged_thread_vcf.pl $@

