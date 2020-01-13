import numpy as np
import MDAnalysis as mda
import os
from numpy.linalg import norm
from scipy.optimize import curve_fit
from MDAnalysis.analysis.rms import rmsd
import matplotlib.pyplot as plt
from scipy import stats

os.chdir('/export/home/im/prof/mariana.ionita/Final/GO/ds_o1/')
u = mda.Universe('em.tpr', 'em.trr')

sel_DNA_not_backbone = u.select_atoms('nucleic and not nucleicbackbone and not resname GGG')
sel_DNA_backbone = u.select_atoms('nucleicbackbone and not resname GGG')
sel_DNA = u.select_atoms('nucleic')
sel_SOL = u.select_atoms('name MG or name NA or name CL or resname SOL')
sel_P1A = u.select_atoms('resname P1A')

sel_stuff = u.select_atoms('(not resname GGG or resname SOL or name NA) or resname E1A or resname H1A or resname C1A')
sel_stuff = {}
for residue in u.residues:
    if ( (str(residue.resname) == 'C1A') or (str(residue.resname) == 'H1A') or (str(residue.resname) == 'E1A') ):
        sel_stuff[str(residue.resname) + "_" + str(residue.resid)] = u.select_atoms('resname ' + str(residue.resname) + ' and resid ' + str(residue.resid))

sel_indep = {}
for residue in sel_DNA_not_backbone.residues:
    for resid in sel_DNA_not_backbone.resids:
        if ( (str(residue)[9:-7] + '_' + str(resid)) not in sel_indep):
            sel_indep[str(residue)[9:-7] + '_' + str(resid)] = sel_DNA_not_backbone.select_atoms('resname ' + str(residue)[9:-7] + ' and resid ' + str(resid))

sel_indep_copy = sel_indep.copy()
for key in sel_indep_copy.keys():
        if (len(sel_indep_copy[key]) == 0):
            sel_indep.pop(key)

backbone_indep = {}
for residue in sel_DNA_backbone.residues:
    for resid in sel_DNA_not_backbone.resids:
        if ((str(residue)[9:-7] + '_' + str(resid)) not in backbone_indep):
            backbone_indep[str(residue)[9:-7] + '_' + str(resid)] = sel_DNA_backbone.select_atoms('resname ' + str(residue)[9:-7] + ' and resid ' + str(resid))

backbone_indep_copy = backbone_indep.copy()
for key in backbone_indep_copy.keys():
    if (len(backbone_indep_copy[key]) == 0):
        backbone_indep.pop(key)

nucleic_residue_indep = {}
for residue in sel_DNA.residues:
    for resid in sel_DNA.resids:
        if ( (str(residue)[9:-7] + '_' + str(resid)) not in nucleic_residue_indep):
            nucleic_residue_indep[str(residue)[9:-7] + '_' + str(resid)] = sel_DNA.select_atoms('resname ' + str(residue)[9:-7] + ' and resid ' + str(resid))

nucleic_residue_indep_copy = nucleic_residue_indep.copy()
for key in nucleic_residue_indep_copy.keys():
    if (len(nucleic_residue_indep_copy[key]) == 0):
        nucleic_residue_indep.pop(key)

sel_CX = u.select_atoms('name CX')
sel_couple = u.select_atoms('name CX or resname E1A or resname C1A or resname H1A or resname P1A')

sel_dict = {}
for key in sel_indep:
    sel_dict["sel_indep_" + str(key)] = sel_indep[key]
for key in backbone_indep:
    sel_dict["backbone_indep_" + str(key)] = backbone_indep[key]
for key in nucleic_residue_indep:
    sel_dict["nucleic_residue_indep_" + str(key)] = nucleic_residue_indep[key]
sel_dict['CX'] = sel_CX
sel_dict['ssdna'] = sel_DNA
sel_dict['backbone'] = sel_DNA_backbone
sel_dict['SOL_and_ions'] = sel_SOL
sel_dict['P1A'] = sel_P1A
for key in sel_stuff:
    sel_dict["sel_stuff_" + str(key)] = sel_stuff[key]
sel_dict['sel_couple'] = sel_couple

filename = "index.ndx"

index_file = mda.selections.gromacs.SelectionWriter(filename, mode='w', numterms=None, preamble=None )

for key in sel_dict.keys():
    index_file.write(sel_dict[key], number=None, name=key, frame=None, mode=None)

index_file.close()


