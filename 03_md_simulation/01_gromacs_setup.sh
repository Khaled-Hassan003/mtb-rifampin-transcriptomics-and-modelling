#!/bin/bash
# 01_gromacs_setup.sh
# GROMACS MD Simulation Pipeline for Rv0784 Complex

# 1. Generate Receptor Topology
gmx pdb2gmx -f rv0784_raw.pdb -o rv0784_processed.gro -water spce

# 2. Define Simulation Box and Solvate System
gmx editconf -f rv0784_processed.gro -o rv0784_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp rv0784_newbox.gro -cs spc216.gro -o rv0784_solv.gro -p topol.top

# 3. Add Ions to Neutralize Charges
gmx grompp -f mdp_files/minimization.mdp -c rv0784_solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o rv0784_ions.gro -p topol.top -pname NA -nname CL -neutral

# 4. Energy Minimization Run
gmx grompp -f mdp_files/minimization.mdp -c rv0784_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# 5. Dual Phase Equilibration (NVT & NPT)
gmx grompp -f mdp_files/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

gmx grompp -f mdp_files/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# 6. Production MD Simulation (65 ns)
gmx grompp -f mdp_files/production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_production.tpr
gmx mdrun -v -deffnm md_production
