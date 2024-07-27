#!/bin/bash                                                                                                                                                                                                
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
kx509                                                                                                                                                                                                     
voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/icarus/Role=Analysis -valid 168:00                                                                                                         
setup icaruscode v09_82_01 -q e26:prof
make
