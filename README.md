# autotalker-reproducibility

## Benchmarking

### SCVI
```conda env create -f scvi/scvi_env.yml``` <br>
```conda activate scvi```

### expiMap
```conda env create -f expimap/expimap_env.yml``` <br>
```conda activate expimap```

### DeepLinc
```conda env create -f deeplinc/deeplinc_env.yml``` <br>
```conda activate deeplinc``` <br>
```python deeplinc.py -e ../dataset/squidpy_seqfish_mouse_organogenesis/counts.csv -a ../dataset/squidpy_seqfish_mouse_organogenesis/adj.csv -c ../dataset/squidpy_seqfish_mouse_organogenesis/coords.csv -r ../dataset/squidpy_seqfish_mouse_organogenesis/cell_types.csv -n run1 -i 40```

### STACI
```conda env create -f staci/staci_env.yml``` <br>
```conda activate staci```
