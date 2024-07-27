#!/bin/bash

while read run_num; do
    Project.py --submit --campaign lkashur_sce_dx_study_run$run_num --stage cryoE
    Project.py --submit --campaign lkashur_sce_dx_study_run$run_num --stage cryoW
done < runlist.txt
