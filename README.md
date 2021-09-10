# wf_RNAvelo
Estimating RNA Velocity through distinct ambients

This workflow was developed by integrating three different tools: [Seurat](https://satijalab.org/seurat/index.html) (R toolkit for single-cell genomics), [velocyto](https://github.com/velocyto-team/velocyto.R) (estimates RNA velocities of single cells by distinguishing unspliced and spliced mRNAs in standard single-cell RNA sequencing protocols), and [scvelo](https://scvelo.readthedocs.io/) (RNA velocity generalized through dynamical model). 

We aimed to design a new pipeline in order to easily unify all the steps need to estimate RNA velocity through distinct tools and models from data coming from different scRNA-seq technologies (10x Chromium, SMART-seq2).

![english workflow](https://user-images.githubusercontent.com/56934471/132842558-ecaa4bf8-514a-4429-a2ae-61ea49776bde.jpg)


The pipeline presents the following main steps:
1. (Pseudo)alignment
2. Basic processing (QC, normalization & dimensional reduction)
3. Clustering
4. Cell subset (optional)
4. DEG
5. RNA Velocity (velocyto.R & scvelo)

See [TFM_AgustinGarcia.pdf](https://github.com/dtin13/wf_RNAvelo/files/7143225/TFM_AgustinGarcia.pdf) (Spanish) for more detailed information about the workflow applied to different tissues.
