# wf_RNAvelo
Estimating RNA Velocity through distinct ambients

This workflow was developed by integrating three different tools: https://satijalab.org/seurat/index.html (R toolkit for single-cell genomics), https://github.com/velocyto-team/velocyto.R (estimates RNA velocities of single cells by distinguishing unspliced and spliced mRNAs in standard single-cell RNA sequencing protocols), and https://scvelo.readthedocs.io/ (RNA velocity generalized through dynamical model). 

We aimed to design a new pipeline in order to easily unify all the steps need to estimate RNA velocity through distinct tools and models from data coming from different scRNA-seq technologies (10x Chromium, SMART-seq2).

![english workflow](https://user-images.githubusercontent.com/56934471/132841661-0c824acf-bdc1-45f4-a8c8-6b4e52c580ef.png)

http://github.com - automatic!
[GitHub](http://github.com)

The pipeline presents the following main steps:
1. (Pseudo)alignment
2. Basic processing (QC, normalization & dimensional reduction)
3. Clustering
4. Cell subset (optional)
4. DEG
5. RNA Velocity (velocyto.R & scvelo)

See [TFM_AgustinGarcia.pdf](https://github.com/dtin13/wf_RNAvelo/files/7143225/TFM_AgustinGarcia.pdf) (Spanish) for more detailed information about the workflow applied to different tissues.
