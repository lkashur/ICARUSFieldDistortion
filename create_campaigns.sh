#!/bin/bash
# GET LIST OF RUN NOS.
for file in lkashur_sceana_run*.cfg; do
    echo ${file//[^0-9]/} >> runlist.txt
done

# LOOP OVER LIST OF RUN NOS.
while read run_num; do

    Project.py --create_campaign --ini_file lkashur_sceana_run$run_num.ini --cfg_file lkashur_sceana_run$run_num.cfg --replace_campaign

done < runlist.txt
