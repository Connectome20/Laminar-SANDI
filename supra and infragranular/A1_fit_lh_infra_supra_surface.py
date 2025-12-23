"""
Main script for generating laminar inferred cortical surfaces using
equidistance, equivolume, and curvature-based linear models

Reference:
Lee H, ..., Huang S, Communications Biology, 2026

Requirements:
- Python 3.11
- FreeSurfer installed and configured
- cortech Python package available (https://github.com/simnibs/cortech)
- nibabel, numpy, pandas
- Subject-specific FreeSurfer recon-all outputs
- fsaverage FreeSurfer template available
- Precomputed global and local layer fraction models (https://github.com/simnibs/cortech)
  (equidistance, equivolume, and linear curvature-based models)

Outputs:
- Inferred laminar surfaces for each subject using:
    * Equidistance (global and local)
    * Equivolume (global and local)
    * Curvature-based linear model
- Output surfaces written in FreeSurfer geometry format
  (with and without CRAS offset)

Author:

Oula Puonti (0000-0003-3186-244X)
Athinoula A. Martinos Center for Biomedical Imaging,
Massachusetts General Hospital

Hansol Lee (0000-0003-2112-1197)
Athinoula A. Martinos Center for Biomedical Imaging,
Massachusetts General Hospital / Harvard Medical School
"""

import sys

sys.path.append(r"/autofs/cluster/connectome2/Bay8_C2/bids/code/SANDI/granular/cortech/")
from pathlib import Path

from cortech.cortex import Hemisphere
from cortech.cortex import Surface
import nibabel as nib
import numpy as np
import numpy.typing as npt
import shutil

from pathlib import Path
from typing import Union
import pandas as pd


def _get_local_inf_prediction(surf, fractions, thickness, curv, placement_method):
    # Fit the surfaces
    n_steps = 40
    fracs = np.linspace(0.2, 0.8, n_steps)

    number_of_nodes = surf.white.n_vertices
    inf_preds = np.zeros((n_steps, number_of_nodes, 3))

    for i, f in enumerate(fracs):
        v = surf.estimate_layers(thickness, curv.H, f, method=placement_method)
        inf_preds[i, ...] = v

    inf_pred = inf_preds[
        np.round(fractions).astype(int), np.arange(0, number_of_nodes), :
    ]
    return inf_pred


# Paths to the various surface models, don't change these

model_path = Path(
    "/autofs/space/rauma_001/users/op035/data/exvivo/hires_surf/analysis/output_smooth20"
)
isovol_global_frac = np.loadtxt(
    model_path / "equivolume_model_smoothed" / "best_fraction_overall.csv"
)
isovol_global_frac = isovol_global_frac.item()
isodist_global_frac = np.loadtxt(
    model_path / "equidistance_model_smoothed" / "best_fraction_overall.csv"
)
isodist_global_frac = isodist_global_frac.item()

isovol_local_frac = (
    nib.load(model_path / "equivolume_model_smoothed" / "lh.best.local.isovol.frac.mgh")
    .get_fdata()
    .squeeze()
)
isodist_local_frac = (
    nib.load(
        model_path / "equidistance_model_smoothed" / "lh.best.local.isodist.frac.mgh"
    )
    .get_fdata()
    .squeeze()
)

# Linear model
beta0 = (
    nib.load(
        model_path / "linear_model_smoothed" / "lh.intercept.all.subjects.smoothed.mgh"
    )
    .get_fdata()
    .squeeze()
)
beta_k1 = (
    nib.load(
        model_path / "linear_model_smoothed" / "lh.beta_k1.all.subjects.smoothed.mgh"
    )
    .get_fdata()
    .squeeze()
)
beta_k2 = (
    nib.load(
        model_path / "linear_model_smoothed" / "lh.beta_k2.all.subjects.smoothed.mgh"
    )
    .get_fdata()
    .squeeze()
)
beta_k1k2 = (
    nib.load(
        model_path / "linear_model_smoothed" / "lh.beta_k1k2.all.subjects.smoothed.mgh"
    )
    .get_fdata()
    .squeeze()
)

betas = [beta0, beta_k1, beta_k2, beta_k1k2]

C2_sub = ["sub-001"
]

for sub in C2_sub:

    # Path to where the FS runs are. This one needs to be changed per data set.
    fs_path = Path(f"/autofs/cluster/connectome2/Bay8_C2/bids/{sub}/anat_process")
    print(fs_path)
    
    # Get fsaverage
    
    if not (fs_path / "fs" / "surf" / "lh.white").exists() or not (fs_path / "fs" / "surf" / "lh.pial").exists():
        continue
    
    fs_sub_path = fs_path / "fs"

    surf = Hemisphere.from_freesurfer_subject_dir(fs_sub_path, "lh")
    
    surf.white.taubin_smooth(n_iter=5, inplace=True)
    surf.pial.taubin_smooth(n_iter=5, inplace=True)
    fsavg = surf.from_freesurfer_subject_dir("fsaverage", "lh")
    
    surf.spherical_registration.compute_projection(fsavg.spherical_registration)
    fsavg.spherical_registration.compute_projection(surf.spherical_registration)
    
    metadata = surf.white.metadata
    
    
    # Calculate curv and thicnkess, then fit models, predict inf surf, compute thickness, curvature (again)
    # and read the intensity values from the scans
    # surf.compute_thickness()
    thickness = surf.compute_thickness()
    curv_wm = surf.white.compute_curvature(smooth_iter=20)
    curv_gm = surf.pial.compute_curvature(smooth_iter=20)
    curv = surf.compute_average_curvature(curv_kwargs={"smooth_iter": 20})
    
    
    # Okay fit isodistance models
    v = surf.estimate_layers(
        thickness, curv.H, isodist_global_frac, method="equidistance"
    )
    
    
    data_path = Path("/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/SANDI/granular")
    sub_path = data_path / sub
    
    sub_path.mkdir()
    
    # Create an output folder
    outpath = sub_path / "isodistance_smooth20"
    
    if outpath.exists():
        shutil.rmtree(outpath)
    
    
    outpath.mkdir()
    
    # Write it out
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.isodistance.global"),
        v + metadata['cras'],
        surf.white.faces,
    )
    
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.isodistance.global.nocras"),
        v,
        surf.white.faces,
    )
    
    # Map isodist fractions to sub space
    # This is actually the frac index
    fracs_on_sub = fsavg.spherical_registration.resample(isodist_local_frac)
    v = _get_local_inf_prediction(
        surf, fracs_on_sub, thickness, curv, "equidistance"
    )
    
    # Write it out
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.isodistance.local"),
        v + metadata['cras'],
        surf.white.faces,
    )
    
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.isodistance.local.nocras"),
        v,
        surf.white.faces,
    )
    
    # Create an output folder
    outpath = sub_path / "isovolume_smooth20"
    
    if outpath.exists():   
        shutil.rmtree(outpath)
    
    
    outpath.mkdir()
    
    # Okay fit isodistance models
    v = surf.estimate_layers(
        thickness, curv.H, isovol_global_frac, method="equivolume"
    )
    
    # Write it out
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.isovolume.global"),
        v + metadata['cras'],
        surf.white.faces,
    )
    
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.isovolume.global.nocras"),
        v,
        surf.white.faces,
    )
    
    # Map isodist fractions to sub space
    # This is actually the frac index
    fracs_on_sub = fsavg.spherical_registration.resample(isovol_local_frac)
    v = _get_local_inf_prediction(
        surf, fracs_on_sub, thickness, curv, "equivolume"
    )
    
    # Write it out
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.volume.local"),
        v + metadata['cras'],
        surf.white.faces,
    )
    
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.volume.local.nocras"),
        v,
        surf.white.faces,
    )
    
    
    # Create an output folder
    outpath = sub_path / "linear_model_smooth20"
    
    if outpath.exists():
        shutil.rmtree(outpath)
    
    
    outpath.mkdir()
    betas_on_sub = []
    
    for b in betas:
        beta_tmp = fsavg.spherical_registration.resample(b)
        betas_on_sub.append(beta_tmp)
    
    betas_on_sub = np.array(betas_on_sub).T
    
    k1 = curv.k1
    k2 = curv.k2
    clip_range = (0.01, 99.9)
    
    pk1 = np.percentile(k1, clip_range)
    pk2 = np.percentile(k2, clip_range)
    k1 = np.clip(k1, pk1[0], pk1[1])
    k2 = np.clip(k2, pk2[0], pk2[1])
    
    k1 = k1 - k1.mean()
    k2 = k2 - k2.mean()
    k1 = k1.T
    k2 = k2.T
    k1k2 = k1*k2
    dummy = np.ones_like(k1)
    predictors = np.stack([dummy, k1, k2, k1k2], axis=1)
    frac = np.einsum("ij, ij -> i", predictors, betas_on_sub)
    # We know the frac can't be <0 or >1
    frac = np.clip(frac, 0, 1)
    v = surf._layer_from_distance_fraction(frac)
    
    # Write it out
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.linear.model"),
        v + metadata['cras'],
        surf.white.faces,
    )
    
    nib.freesurfer.write_geometry(
        outpath / Path("lh.inf.linear.model.nocras"),
        v,
        surf.white.faces,
    )
