#!/usr/bin/env python3
# 02_trajectory_analysis.py
# Structural Trajectory Kinetics and Protein-Ligand Interaction Fingerprints (ProLIF)

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import prolif as plf
import plotly.graph_objects as go
import numpy as np

# 1. Load Molecular Dynamics Trajectory
u = mda.Universe("topol.top", "md_production.xtc")
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname UNL")

# 2. Compute 1D Backbone RMSD Calculation
R = rms.RMSD(u, select="backbone")
R.run()

# Plot RMSD Kinetics
fig = go.Figure()
fig.add_trace(go.Scatter(x=R.results.rmsd[:, 1]/1000, y=R.results.rmsd[:, 2], mode='lines', name='Backbone'))
fig.update_layout(title="Rv0784 Backbone RMSD Stability Curve", xaxis_title="Time (ns)", yaxis_title="RMSD (nm)")
fig.write_html("03_md_simulation/plots/rmsd_stability.html")

# 3. Protein-Ligand Interaction Fingerprints via ProLIF
ligand_plf = u.atoms.select_atoms("resname UNL")
protein_plf = u.atoms.select_atoms("protein")

fp = plf.Fingerprint()
fp.run(u.trajectory[::10], ligand_plf, protein_plf)
df = fp.to_dataframe()

# Filter and save high occupancy contacts (>90%)
occ = fp.to_dataframe().mean().to_dict()
stable_contacts = {k: v for k, v in occ.items() if v >= 0.90}
print("Conserved Binding Interactions Matrix:", stable_contacts)
