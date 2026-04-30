#!/bin/bash

outdir="../../data/ovarian-carcinoma/spaceranger"
mkdir -p $outdir

# Download files from 10x
wget https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
wget -qO- https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_spatial.tar.gz | tar -xz
wget -qO- https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma/CytAssist_11mm_FFPE_Human_Ovarian_Carcinoma_filtered_feature_bc_matrix.tar.gz | tar -xz

# Move to outdir
mv Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv $outdir/
mv spatial $outdir/
mv filtered_feature_bc_matrix $outdir/