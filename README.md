# NicheCompass Reproducibility

## Installation
1) Clone the nichecompass-reproducibility repository and navigate into it: <br>
```git clone https://github.com/sebastianbirk/nichecompass-reproducibility.git``` <br>
```cd nichecompass-reproducibility```

2) Create the nichecompass-reproducibility conda environment: <br>
```conda env create -f envs/environment.yml```

3) Create the deeplinc conda environment (for deeplinc workloads which are
incompatible with rest): <br>
```conda env create -f envs/environment_deeplinc.yml```

4) Activate the nichecompass-reproducibility conda environment: <br>
```conda activate nichecompass-reproducibility```

5) Clone the nichecompass repository and navigate into it: <br>
```git clone https://github.com/sebastianbirk/nichecompass.git``` <br>
```cd nichecompass```

6) Install all NicheCompass Python dependencies via Poetry: <br>
```poetry install```



