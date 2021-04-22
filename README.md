# PMd_escape_vigor
This repository contains the analysis code for the manuscript "Dorsal premammillary hypothalamic projection to periaqueductal gray controls escape vigor from innate and conditioned threats."

Regarding a fMRI-specific analysis code, specifically:
bootpls_loadings.m
bootpls_yloadings.m
bootpls.m
extend_range.m
hyperalign_group_avg.m
hyperalign_pls.m
hyperalign.m
lighten.m
linspecer.m
model_brain_pathway_v1.m
model_hythal_PAG_pathways.m
plot_coords.m
r_test_paired.m
rgb.m
trajectory_plotter.m
vals2colors.m
varexp_diff.m
varexp.m


This code estimates a pathway from the hypothalamus to the PAG and
assess its functional association with different stimulus modalities and
behavioral ratings of avoidance - see model_hythal_PAG_pathways.m to run
the analyses reported in the manuscript.

Required MATLAB toolboxes
SPM (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)  
CANLab core tools (https://github.com/canlab/CanlabCore) 

Required data
The data file MPA_data_pHythal_PAG.mat (located at https://datadryad.org/stash/share/dYuSl2nnXsyi0nTDjCDeHR08gwW7paFL4Eo3TmF_aH4) contains the following variables:

CeM: a mask of the central amygdala from the SPM Anatomy toolbox - https://github.com/inm7/jubrain-anatomy-toolbox
PAG: a mask of the PAG from Kragel et al. 2019 - https://doi.org/10.1523/JNEUROSCI.2043-18.2019
hythal: a mask of the hypothalamus from the CIT168 Atlas - https://doi.org/10.1038/sdata.2018.63
mask: a combined mask of these regions
masked_dat: a CANLab data object for the single-trial estimates of fMRI activation within the above regions
S: a vector denoting subject number for each trial
COND: a vector denoting the condition (stimulus modality) for each trial
