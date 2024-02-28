# Analysis of the STARmap PLUS Mouse Central Nervous System spatial dataset with NicheCompass

## Training of the NicheCompass model

Model training was performed using NicheCompass (branch `59-add-clustering`) in a Singularity container, with the following commands and with the config described in `run-config.json`.

```
nichecompass build-gene-programs run-config.json
nichecompass build-dataset run-config.json
nichecompass train-reference run-config.json
nichecompass cluster run-config.json
```

## Alignment of the STARmap dataset with the Allen Mouse Brain Reference Atlas

Alignment of the STARmap dataset with the Allen Mouse Brain Reference Atlas was performed using `starmap_mouse_central_nervous_system_atlas_alignment.ipynb` with a python environment described in `requirements_2.txt`.

## Downstream interpretation of the NicheCompass model output

Downstream figure generation was performed using `starmap_mouse_central_nervous_system.ipynb` with a python environment described in `requirements_1.txt`.
