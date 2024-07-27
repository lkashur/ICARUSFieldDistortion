# ICARUSFieldDistortion
Standalone C++ and python scripts to analyze space charge effects and other E-Field distortions at ICARUS.

### Prerequisite
Set up the proper container on an ICARUS GPVM:
```
sh /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh
```

Create working directory that corresponds to a dataset
```bash
mkdir <dataset>
```

### Stage 0: Calculate Offsets
<ins>Input</ins>: Calibration Ntuples

<ins>Output</ins>: ROOT file(s) with TTree containing x,y,z and dx info for every hit


##### Setup
```bash
cd stage0
source setup.sh
```

##### Run analyzer
From dataset directory:
```bash
../stage0/sce_dx_analyzer filelist0.txt W
```

The above command produces an output ROOT file(s) named "offsets_W.root" containing the offset TTree.

### Stage 1: Aggregate Offsets
<ins>Input</ins>: ROOT files from Stage 0 (paths in a text file)

<ins>Output</ins>: 
1. TH3 with aggregated offset in each bin
2. Text file with list of offsets for each bin

##### Setup
```bash
cd stage1
source setup.sh
```

##### Run aggregator
From dataset directory:
```bash
../stage0/sce_dx_aggregator filelist1.txt WE
```

The above command produces an output ROOT file named "aggoffsets_WE.root" containing the aggregated TH3, as well as a text file named "alloffsets_WE.txt" containing all offsets per bin.

### Stage 2: Generate Space Charge Density Maps
