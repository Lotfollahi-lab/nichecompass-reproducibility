# NicheCompass Reproducibility

This repository contains the code to reproduce the analyses and benchmarking experiments performed in the NicheCompass [manuscript](https://www.biorxiv.org/content/10.1101/2024.02.21.581428v3).
The NicheCompass source code can be found [here](https://github.com/Lotfollahi-lab/nichecompass).

## Installation

### Standard
1) Clone the nichecompass-reproducibility repository and navigate into it: <br>
```git clone https://github.com/Lotfollahi-lab/nichecompass-reproducibility.git``` <br>
```cd nichecompass-reproducibility```

2) (Optional) Install the Libmamba solver to make the installation faster: <br>
```conda update -n base conda``` <br>
```conda install -n base conda-libmamba-solver``` <br>
```conda config --set solver libmamba```

3) Create the nichecompass-reproducibility conda environment: <br>
```conda env create -f envs/environment.yaml```

4) Create the deeplinc conda environment (for benchmarking deeplinc method which relies on legacy packages and is
incompatible with dependencies of other methods): <br>
```conda env create -f envs/environment_deeplinc.yaml```

5) Activate the nichecompass-reproducibility conda environment: <br>
```conda activate nichecompass-reproducibility```

6) Install NicheCompass: <br>
Follow the instructions in the [documentation](https://nichecompass.readthedocs.io/en/latest/installation.html).
 
### Docker / Charliecloud Container
1) Clone the nichecompass-reproducibility repository: <br>
```git clone https://github.com/sebastianbirk/nichecompass-reproducibility.git``` <br>

2) Clone the nichecompass repository: <br>
```git clone https://github.com/sebastianbirk/nichecompass.git```

3) From the root repository that contains both the nichecompass and nichecompass-reproducibility repositories, run: <br>
```docker buildx build --load --platform linux/amd64 --file nichecompass-reproducibility/envs/Dockerfile --tag nichecompass . --no-cache``` (macOS)

4) (For Charliecloud) Export the Docker image to a tarball: <br>
```docker export $(docker create nichecompass) | gzip -c > ./nichecompass.tar.gz```

## Data & Models
All preprocessed data used in the manuscript and trained models are downloadable from [GDrive](https://drive.google.com/drive/folders/1sqoqCq1y5NMIbC1K7uq6v4PBWDPQQJgY).

## Reference
```
@article{Birk2024,
  title = {Quantitative characterization of cell niches in spatial atlases},
  author = {Birk, S. and others},
  journal = {bioRxiv},
  year = {2024},
  doi = {10.1101/2024.02.21.581428},
  url = {https://www.biorxiv.org/content/10.1101/2024.02.21.581428}
}
```


