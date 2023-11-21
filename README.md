# MicroFreqTool
This toolbox was developped to 

(1) investigate the microstructure-induced frequency shift in white matter (WM) with crossing fibers using hollow cylinder fiber model (HCFM)

(2) separate the microstructure-related frequency shift from the bulk susceptibility-induced frequency shift by model fitting the gradient-echo frequency evolution (i.e.,TE-dependent frequency fitting)

Refer to this paper:

Exploiting gradient-echo frequency evolution: probing white matter microstructure and extracting bulk susceptibility-induced frequency for quantitative susceptibility mapping, Magnetic Resonance in Medicine, Chen et al.

This toolbox has the following features:

(1) Create 3D white models with HCFM with two fiber orientations.
    
    See Build2FiberHCFM folder for generation of each model inside WM voxel
    Intra-voxel field perturbations and the complex signal were simulated in qsub_SimulateHCFMSignalVoxelParallel_batch.m

(2) Simulate the field perturbation and the corresponding complex       signal from 3D WM model.
     
     See SimulateSignal folder
     See ShowSignalEvolutionExample.m for an example

(3) Simulate full complex GRE signal inside brain with total frequency shift f_total=f_micro+f_macro and conduct TE-dependent frequency fitting.
    
    See BrainPhantom folder

(4) TE-dependent frequency fitting with regularization for in vivo data.

    See Invivo folder
    
