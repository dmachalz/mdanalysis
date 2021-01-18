# This script uses MDAnalysis to count the number of C-alpha atoms at the interface of two proteins A and B. 
# author @dmachalz
import MDAnalysis
import os
import numpy as np

# Variables
# The working directory containing the PDB files. Proteins A and B should be stored in one PDB file.
foldername='REPLACEME'
# Outputfile in CSV format
outcsv='REPLACEME.csv'
# Project required renaming
replacer='1jma_prep_corr_orient_'

outlist=[]
# Change to working directory
os.chdir(foldername)
pdblist = [p for p in os.listdir(foldername) if p.endswith('.pdb')]
for pdbfile in pdblist:
    # Getting the System props - highly project specific
    sysname = pdbfile.replace('.pdb','')
    sysnumber = sysname.replace(replacer, '')
    # Create MDAnalysis universe
    u = MDAnalysis.Universe(pdbfile)
    # Get atom selection
    asel = u.select_atoms("((name CA and segid A) and around 4.5 segid B) or ((name CA and segid B) and around 4.5 segid A)")
    # Get number of Calpha atoms (acount) in selection
    acount = len(asel)
    # Write pdbname, and acount to csv file
    #print(sysname, acount)
    outlist.append([int(sysnumber), acount])

# Converting output list to numpy array for easier writing
outlist = np.array(outlist)
# Actually writing out the list using 'integer' formatting and delimiter ',' to OUTCSV!
np.savetxt(outcsv, outlist, fmt='%i', delimiter=',')
