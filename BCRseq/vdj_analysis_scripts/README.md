# Scripts for BCR analysis
Kenneth B. Hoehn
kenneth.hoehn@yale.edu
3/8/2022

# Dependencies

OS:
Ubuntu 20.04 LTS, may work on other systems.

Docker images:
immcantation/suite:4.0.0

R 4.1.1 packages:
dowser v0.1.0
scoper v1.1.0
alakazam v1.2.0
shazam v1.1.0
ggtree v3.0.4
ape v5.5-3
dplyr v1.0.7
tidyr v1.1.4
ggpubr v0.4.0
ggplot2 v3.3.5

Other:
IgPhyML v1.1.3

Directory structure:
raw/Pam1_CITE_multi5P06/vdj_b/
raw/Pam1_CITE_multi5P07/vdj_b/
raw/Pam1_CITE_multi5P08/vdj_b/
intermediates
processed
results

# Processing VDJ data

To align sequences to the VDJ reference database, the the Immcantion docker image:

```bash
docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.0.0 bash
```

You may need to use sudo. Within the Docker image, run the processing pipeline:

```bash
bash processData.sh
```

Now, exit the Docker container.

# Downstream analysis

To combine information from the Seurat object and the VDJ sequences, run:

```bash
Rscript combineData.R
```

To perform clonal clustering, run:

```bash
Rscript cloneGermline.R
```

To perform all downstream analyses and generate all BCR figures, run:

```bash
Rscript analysis.R
```
