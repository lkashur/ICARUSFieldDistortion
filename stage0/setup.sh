#!/bin/bash 
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
source /exp/sbnd/app/users/vito/products/setup.sh
voms-proxy-destroy
kx509 
#voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/icarus/Role=Analysis -valid 168:00
setup icaruscode v09_82_01 -q e26:prof
#setup project_py v2_2_2
setup project_py v2_3_0
make


