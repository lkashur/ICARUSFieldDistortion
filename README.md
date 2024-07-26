# ICARUSFieldDistortion
Standalone C++ and python scripts to analyze space charge effects and other E-Field distortions at ICARUS.

### Prerequisite
Set up the proper container on an ICARUS GPVM:
```
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh
```

### Stage 0: Calculate Offsets
<ins>Input</ins>: Calibration Ntuples

<ins>Output</ins>: TTree(s) containing x,y,z and dx info for every hit


##### Setup
```bash
cd stage0
source setup.sh
```

##### Create working directory that corresponds to a dataset
```bash
cd ../
mkdir <dataset name>
```

##### Run analyzer
```bash
../stage0/sce_dx_analyzer filelist.txt W
```

The above command produces an output root file(s) named "offsets_W.root" containing the offset TTree.

### Stage 1: Aggregate Offsets

### Stage 2: Generate Space Charge Density Maps
