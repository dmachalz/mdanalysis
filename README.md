# mdanalysis

Collection of scripts to be used with MDAnalysis. MDAnalysis allows to analyse Molecular Dynamics simulations in python. (https://github.com/MDAnalysis)

## perimeter.py
Monitor the perimeter between a selection of atoms for a MD simulation

## mdanalysis/main.py & mdanalysis/config/config-template.cfg
This provides config template that can be used to execute standard analysis on MD simulation trajectories via main.py, which again applies MDAnalysis.  
For atom selection syntax of MDAnalysis see https://docs.mdanalysis.org/stable/documentation_pages/selections.html

### Monitoring of atom distances
To monitor more than one pair of atoms state the selection command of atom pairs separated by semicolon:  
selectionlist = \<atomselection1>, \<atomselection2>; \<atomselection3>, \<atomselection4>  
labellist = atom1-atom2, atom3-atom4

### Example
python main.py config.cfg
