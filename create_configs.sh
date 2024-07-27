#!/bin/bash
# FIRST LOOP TO SPLIT LARGE TEXT FILE
split --suffix-length=4 -d -C 500000 $1
for file in x*; do
    mv $file ${file}.txt
done

# SECOND LOOP TO GENERATE CFG/INI FILES
for file in x*; do
    RUN=$(echo ${file//[^0-9]/})

    cat > "lkashur_sceana_run$RUN.cfg" << EOF
[global]
# these are global parameters that can be used in any section
# and eventually overridden in stage sections
# users can add more parameters as needed to be used through the configuration
experiment = icarus
account = lkashur
project_name = sce_dx_study_run$RUN
run = $RUN
group = %(experiment)s
wrapper = file:///\${FIFE_UTILS_DIR}/libexec/fife_wrap
b_name = %(project_name)s
basename = override_me
outdir = override_me
stage_name = override_me

[env_pass]
# environment variable to pass to jobs
IFDH_DEBUG = 1
SAM_EXPERIMENT = %(experiment)s
# these are used by Project-py to identify output dataset
# XRootD env var to make it more resilient
XRD_CONNECTIONRETRY = 32
XRD_REQUESTTIMEOUT  = 14400
XRD_REDIRECTLIMIT   = 255
XRD_LOADBALANCERTTL = 7200
XRD_STREAMTIMEOUT   = 7200

[submit]                                       
# jobsub_submit options                                                 
G = %(group)s                                            
e = SAM_EXPERIMENT                                              
e_1 = IFDH_DEBUG                                     
e_2 = POMS4_CAMPAIGN_NAME                               
e_3 = POMS4_CAMPAIGN_STAGE_NAME                                 
generate-email-summary = True                               
email-to = %(account)s@fnal.gov
# default resource requirements, these can be overridden in stage sections                                        
expected-lifetime = 8h                                                   
disk = 10GB                                               
memory = 1500MB                                               
blacklist = RAL
singularity-image = /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest

[job_setup]
# options related to user code setup
ifdh_art = True
debug = True
find_setups = True
setup_local = True
source_1 = /cvmfs/%(experiment)s.opensciencegrid.org/products/%(experiment)s/setup_%(experiment)s.sh 
setup_1 = icaruscode v09_82_01 -q e26:prof

[sam_consumer]
# SAM related options
limit = 1
schema = root
appname = %(experiment)scode
appfamily = art

[job_output]
# options to handle ROOT output files
addoutput = offsets*
dest = %(outdir)s

[executable]
# executable options/arguments
# most of them are overridden in stage sections
name = true
arg_1 = overrideme
arg_2 = overrideme

# In stage sections we can override parameter from each section using the form:
# section.key = value
# most options are using placeholder values
# resources need to be assigned realistic values

[stage_cryoE]
global.stage_name = cryoE
job_setup.ifdh_art = False
job_setup.prescript = ls -lh \${CONDOR_DIR_INPUT}
job_setup.prescript_1 = chmod +x \${CONDOR_DIR_INPUT}/sce_dx_analyzer
global.outdir = /pnfs/%(experiment)s/scratch/users/%(account)s/%(project_name)s
executable.name = \\\${CONDOR_DIR_INPUT}/sce_dx_analyzer
executable.arg_1 = \\\${CONDOR_DIR_INPUT}/$file
executable.arg_2 = E
submit.f = dropbox:///exp/icarus/app/users/lkashur/sce_poms/step0/sce_dx_analyzer
submit.f_1 = dropbox://$PWD/$file
submit.N = 1
submit.disk = 2GB
submit.expected-lifetime = 4h
submit.memory = 2000MB

[stage_cryoW]
global.stage_name = cryoW
job_setup.ifdh_art = False
job_setup.prescript = ls -lh \${CONDOR_DIR_INPUT}
job_setup.prescript_1 = chmod +x \${CONDOR_DIR_INPUT}/sce_dx_analyzer
global.outdir = /pnfs/%(experiment)s/scratch/users/%(account)s/%(project_name)s
executable.name = \\\${CONDOR_DIR_INPUT}/sce_dx_analyzer
executable.arg_1 = \\\${CONDOR_DIR_INPUT}/$file
executable.arg_2 = W
submit.f = dropbox:///exp/icarus/app/users/lkashur/sce_poms/step0/sce_dx_analyzer
submit.f_1 = dropbox://$PWD/$file
submit.N = 1
submit.disk = 2GB
submit.expected-lifetime = 4h
submit.memory = 2000MB
EOF

    cat > "lkashur_sceana_run$RUN.ini" <<EOF
[campaign]
experiment = icarus
poms_role = analysis
name = lkashur_sce_dx_study_run$RUN
state = Active
campaign_stage_list = cryoE,cryoW

[campaign_defaults]
vo_role = Analysis
software_version = v09_82_01
dataset_or_split_data = from_parent
cs_split_type = None
completion_type = complete
completion_pct = 99
param_overrides = []
test_param_overrides = []
login_setup = lkashur_login
job_type = lkashur_sce_dx_study_run${RUN}_jobtype

[campaign_stage cryoE]
dataset_or_split_data = None
cs_split_type = None
completion_type = complete
completion_pct = 99
param_overrides = [["--stage ", "cryoE"]]
test_param_overrides = [["--stage ", "cryoE"]]
login_setup = lkashur_login
job_type = lkashur_sce_dx_study_run${RUN}_jobtype

[campaign_stage cryoW]
dataset_or_split_data = None
cs_split_type = None
completion_type = complete
completion_pct = 99
param_overrides = [["--stage ", "cryoW"]]
test_param_overrides = [["--stage ", "cryoW"]]
login_setup = lkashur_login 
job_type = lkashur_sce_dx_study_run${RUN}_jobtype

[login_setup lkashur_login_token]
host = pomsgpvm02.fnal.gov
account = poms_launcher
setup = source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh;setup fife_utils

[job_type lkashur_sce_dx_study_run${RUN}_jobtype]
launch_script = fife_launch
parameters = [["-c ", "$PWD/lkashur_sceana_run${RUN}.cfg"]]
output_file_patterns = %%.root

#[dependencies cryoW]
#campaign_stage1 == cryoE
#file_pattern_1 = %%
EOF

done
