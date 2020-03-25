#/bin/bash
set -euo pipefail 

head -n 10000 data/ccs/CD4_stim_DS75687.ccs.bam.zmws > data/test_data/test.zmws

time /net/eichler/vol26/projects/sequencing/pacbio/smrt-link/smrtcmds/bin/bamsieve --whitelist data/test_data/test.zmws data/ccs/CD4_stim_DS75687.ccs.bam data/test_data/ccs.bam

time /net/eichler/vol26/projects/sequencing/pacbio/smrt-link/smrtcmds/bin/bamsieve --whitelist data/test_data/test.zmws data/subreads/CD4_stim_DS75687.subreads.bam data/test_data/subreads.bam

