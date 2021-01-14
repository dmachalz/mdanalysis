"""This scripts reads in the conf file and executes the appropriate MDAnalysis functionalities"""
import os
import sys
import configparser
import warnings

sys.path.extend(['/mdspace/davidm/tools/mdanalysis'])
from lib.functions import *
from MDAnalysis.analysis import distances
warnings.filterwarnings("ignore")  # This is bad programming: ignoring warnings

#configpath = '/home/davidm/Tools/MDanalysis/run-mda-test-single.cfg'
configpath = None

### Reading in the config file
# Reading in config file
if len(sys.argv) == 2 and os.path.exists(sys.argv[1]):
    configpath = sys.argv[1]  # Config file path provided via command line
    print('')
    print('Using config file {}'.format(sys.argv[1]))
elif configpath is not None:
    pass  # this is used for the case, that the configfile is defined within the script
else:
    config = None
    print('No (suitable) config file provided... Exiting')
    exit()

config = configparser.ConfigParser()
config.sections()
config.read(configpath)


### Checking input files ###
if config is not None:
    if config.has_section('Reading'):
        if config.has_option('Reading', 'topology') and config.has_option('Reading', 'trajectory'):
            topology = config.get('Reading', 'topology')
            trajectory = config.get('Reading', 'trajectory')
            if os.path.isfile(topology) is False or os.path.isfile(trajectory) is False:
                print('Please provide valid input MD file paths in the conf file... Exiting')
                exit()
                topology = trajectory = None
        else:
            print('Please provide valid input paths in the conf file... Exiting')
            exit()
            topology = trajectory = None
    else:
        print('Please provide valid input path section in the conf file... Exiting')
        exit()
        topology = trajectory = None

    ### Looking for functionalities to execute ###
    if config.has_section('RMSD'):
        if config.has_option('RMSD', 'proteinbb') and config.get('RMSD', 'proteinbb') == '1':
            runprotrmsd = True
            print('Protein backbone RMSD calculation active')
        else:
            runprotrmsd = False
        if config.has_option('RMSD', 'ligand') and config.get('RMSD', 'ligand') == '1':
            ligandname = config.get('RMSD', 'ligandname')
            if ligandname is not None and ligandname != '':
                runligrmsd = True
                print('Ligand heavy atoms RMSD calculation active')
            else:
                runligrmsd = False
                print('Ligand selection command is empty... Exiting')
                exit()
        else:
            runligrmsd = False
    else:
        runprotrmsd = runligrmsd = False

    if config.has_section('RMSF'):
        if config.has_option('RMSF', 'proteinbb') and config.get('RMSF', 'proteinbb') == '1':
            runprotrmsf = True
            print('Protein backbone RMSF calculation active')
        else:
            runprotrmsf = False
    else:
        runprotrmsf = False

    if config.has_section('Distance monitoring'):
        if config.has_option('Distance monitoring', 'selectionlist') and \
                config.get('Distance monitoring', 'selectionlist') is not None:
                    selectionlist = config.get('Distance monitoring', 'selectionlist')
                    selectionlist = [i.split(', ') for i in selectionlist.split('; ')]
        else:
            selectionlist = None
        if config.has_option('Distance monitoring', 'labellist') and \
            config.get('Distance monitoring', 'labellist') is not None:
            labellist = config.get('Distance monitoring', 'labellist').split(', ')
        else:
            labellist = None
    else:
        selectionlist = labellist = None

    if config.has_section('Plotting'):
        if config.has_option('Plotting', 'plot-rmsd') and config.get('Plotting', 'plot-rmsd') == '1':
            plotprotrmsd = True
        else:
            plotprotrmsd = False
        if config.has_option('Plotting', 'plot-rmsf') and config.get('Plotting', 'plot-rmsf') == '1':
            plotprotrmsf = True
        else:
            plotprotrmsf = False
        if config.has_option('Plotting', 'plot-dist') and config.get('Plotting', 'plot-dist') == '1':
            plotdist = True
        else:
            plotdist = False
    else:
        plotprotrmsd = plotprotrmsf = plotdist = False
    # TODO: Include plotting of ligand rmsd

    if config.has_section('Writing'):
        if config.has_option('Writing', 'outdir') and config.get('Writing', 'outdir') is not '':
            writeoutdir = config.get('Writing', 'outdir')
            writeoutdir = writeoutdir.rstrip("/")

            if os.path.isdir(writeoutdir) is False:
                print('Creating output directory under')
                os.makedirs(writeoutdir)
        else:
            writeoutdir = None
    else:
        writeoutdir = None
else:
    print('No config provided... Exiting')
    exit()
    runligrmsd = runprotrmsd = runprotrmsf = False
    topology = trajectory = writeoutdir = labellist = selectionlist = None
    plotprotrmsd = plotprotrmsf = plotdist = False

# Defining output file names

if writeoutdir is not None:
    rmsfout = writeoutdir + '/mda-rmsf-allframes.dat'
    rmsdout = writeoutdir + '/mda-rmsd-allframes.dat'
    rmsdligout = writeoutdir + '/mda-rmsd-lig-allframes.dat'
else:
    rmsfout = rmsdout = rmsdligout = None

if selectionlist is not None and labellist is not None:
    if len(selectionlist) == len(labellist):
        rundistmonitor = True
    else:
        print('Distance monitoring selection list and data label list are not same length! Exiting')
        exit()
        labellist = selectionlist = None
        rundistmonitor = False
else:
    rundistmonitor = False
print('')

### Running MDAnalysis ###
print('Constructing MDAnalysis universe')
u = MDAnalysis.Universe(topology, trajectory)  # Creating the MDAnalysis universe

protein = u.select_atoms("protein")  # Creating atom selection of protein atoms

resnames = protein.atoms.residues.resnames  # numpy array of residue names
resids = protein.atoms.residues.resids  # numpy array of residue IDs
reslist = list(zip(resnames, resids))

# Executing Functionalities
if runprotrmsd is True:
    print('\
Calculating protein backbone RMSD')
    rmsdcalc(u, rmsdout)  # calculating protein backbone rmsd

if runprotrmsf is True:
    print('\
Calculating protein backbone RMSF')
    rmsfcalc(u, rmsfout)  # calculationg the backbone rmsf

if runligrmsd is True:
    print('\
Calculating ligand heavy atom RMSD')
    rmsdcalclig(u, rmsdligout=rmsdligout, ligandsel=ligandname)  # calculating ligand heavy atom rmsd

if rundistmonitor is True:
    for distsel, datalabel in zip(selectionlist, labellist):
        sel1 = distsel[0]
        sel2 = distsel[1]
        distmonitorout = '{0}/mda-{1}-dist.dat'.format(writeoutdir, str(datalabel))
        print('Getting {} - distance'.format(datalabel))
        atomdist(u, sel1, sel2, outputpath=distmonitorout)

print('Done with analysis. :)')

if runprotrmsd is True and plotprotrmsd is True:
    print('Plotting...')
    rmsdplot()
if plotprotrmsf is True:
    print('Plotting...')
if plotdist is True:
    print('Plotting...')