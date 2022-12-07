# !/bin/bash

time docker run -v /hihg:/hihg -v /Volumes/Synology:/Volumes/Synology -v `pwd`:/app -it hihg-um/$USER/perl scripts/make_20K_multi_x_flagged_thread_vcf.pl $@

