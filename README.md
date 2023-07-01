# NicheCompass Reproducibility

## Installation

### Standard
1) Clone the nichecompass-reproducibility repository and navigate into it: <br>
```git clone https://github.com/sebastianbirk/nichecompass-reproducibility.git``` <br>
```cd nichecompass-reproducibility```

2) (Optional) Install the Libmamba solver to make the installation faster: <br>
```conda update -n base conda``` <br>
```conda install -n base conda-libmamba-solver``` <br>
```conda config --set solver libmamba```

3) Create the nichecompass-reproducibility conda environment: <br>
```conda env create -f envs/environment.yml```

4) Create the deeplinc conda environment (for benchmarking deeplinc method which relies on legacy packages and is
incompatible with dependencies of other methods): <br>
```conda env create -f envs/environment_deeplinc.yml```

5) Activate the nichecompass-reproducibility conda environment: <br>
```conda activate nichecompass-reproducibility```

6) Clone the nichecompass repository and navigate into it: <br>
```cd ..``` <br>
```git clone https://github.com/sebastianbirk/nichecompass.git``` <br>
```cd nichecompass```

7) Install all NicheCompass Python dependencies via Poetry: <br>
```poetry install```

### Docker / Charliecloud Container
1) Clone the nichecompass-reproducibility repository: <br>
```git clone https://github.com/sebastianbirk/nichecompass-reproducibility.git``` <br>

2) Clone the nichecompass repository: <br>
```git clone https://github.com/sebastianbirk/nichecompass.git``` <br>

3) From the root repository that contains both the nichecompass and nichecompass-reproducibility repositories, run: <br>
```docker buildx build --load --platform linux/amd64 --file nichecompass-reproducibility/envs/Dockerfile --tag nichecompass . --no-cache``` <br>



