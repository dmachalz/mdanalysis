### This file contains all functions used for common MDAnalysis ###

import MDAnalysis
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from unittest.mock import MagicMock
import os


def getmatchingfiles(pattern1, pattern2, filelist):
    """
    From a list of files get a new sorted list of files matching both patterns 1 and 2
    :param pattern1:
    :param pattern2:
    :param filelist:
    :return:
    """
    return [file for file in sorted(filelist) if pattern1 in file and pattern2 in file]


def stderr(data, N):
    """
    Standard Error (SE) calculation
    :param data: Data to perform the standard deviation on
    :param N:    Number of samples
    :return:     Standard error value
    """
    SE = np.std(data) / np.sqrt(N)
    return SE


def cm2inch(value):
    """
    Convert centimeters to inches for figure sizes
    :param value:   Value in centimeters
    :return:        Value in inches
    """
    return value/2.54


def trajload(top, trj):
    """
    Loading trajectories. This is not really used.
    :param top: Topology file
    :param trj: Trajectory file
    :return: Returns a MDAnalysis universe
    """
    u = MDAnalysis.Universe(top, trj)  # creating the universe

    protein = u.select_atoms("protein")  # creating atom selection of protein atoms

    resnames = protein.atoms.residues.resnames  # numpy array of residue names
    resids = protein.atoms.residues.resids  # numpy array of residue IDs
    reslist = list(zip(resnames, resids))
    return u


def rmsdcalc(u, rmsdout):
    """
    Calculation of Protein RMSD
    :param u: MDAnalysis universe
    :param rmsdout: File path to output file
    :return: None
    """
    prot = u.select_atoms("protein")
    u.trajectory[0]
    R = RMSD(prot, select="backbone", filename=rmsdout)
    R.run()
    np.savetxt(rmsdout, R.rmsd, delimiter=',')
    print('Finished Protein RMSD calculation successfully!')


def rmsdcalclig(u, rmsdligout=None, ref=None, ligandsel=None):
    """
    Calculaiton of ligand RMSD
    :param u:           MDAnalysis universe
    :param rmsdligout: File path to output file
    :param ref:         Reference in the format of a universe with identical atom selection!
    :param ligandsel:   Selection command for ligand
    :return:            None
    """
    if ligandsel is None:
        ligand = u.select_atoms("not protein and not ((resname T3P or resname C*) or (resname N* or resname HEM))")
    else:
        ligand = u.select_atoms("{}".format(ligandsel))
    ligandheavy = ligand.select_atoms("not name H*")

    u.trajectory[0]
    if ref is not None:
        Rlig = RMSD(ligandheavy, reference=ref, select='all')  # output: frame, time (ps), RMSD (A)
    elif ref is None:
        Rlig = RMSD(ligandheavy, select='all')  # output: frame, time (ps), RMSD (A)
    Rlig.run()
    raw = Rlig.run().rmsd
    if rmsdligout is not None:
        np.savetxt(rmsdligout, raw, delimiter=',')
    print('Finished Ligand RMSD calculation successfully!')

    return raw


def rmsfcalc(u, rmsfout, modus="calphas"):
    """
    Calculation of Protein RMSF
    :param u:       MDAnalysis universe
    :param rmsfout: File path to output file
    :param modus:   Selection to align and calculation ("calphas" or "allatoms")
    :return:    None
    """
    calphas = u.select_atoms("protein and name CA")  # defining alpha atom selection modus
    allatoms = u.select_atoms("protein")            # defining all all heavy atom selection modus
    if modus is "calphas":
        select = calphas
    elif modus is "allatoms":
        select = allatoms

    rmsfer = RMSF(select, verbose=True, start=0, stop=50).run()  # defining and executing the rmsf calculation

    rawrmsf = rmsfer.rmsf  # saving the output for writing

    # writing out the obtained data to file

    with open(rmsfout, 'w+') as out:
        for line in zip(calphas.residues, rawrmsf):  # labeling the RMSF value with residue name and id
            newline = str(line).replace('<Residue', '').replace('>', '').replace('(', '')\
                .replace(')', '')  # making the string nice
            out.write(str(newline) + '\n')  # writing the string
        out.close()
    print('Finished Protein RMSF calculation successfully!')


def atomdist(u, sel1, sel2, outputpath=None):
    """
    Atom distance monitoring
    :param u: MDAnalysis atom universe
    :param sel1: Atom selection statement 1
    :param sel2: Atom selection statement 2
    :param outputpath: If given, the computed numpy array will be written to this location
    :return: numpy distance array
    """

    a = u.select_atoms("{}".format(str(sel1)))
    b = u.select_atoms("{}".format(str(sel2)))
    distlist = []
    step = -1
    for ts in u.trajectory:
        step += 1
        # distance array cropped to the distance between selections rounded to two decimal digits
        dist = round(MDAnalysis.analysis.distances.dist(a, b)[2, 0], 2)
        distlist.append((step, u.trajectory.time, dist))
    distarray = np.array(distlist)
    if outputpath is not None:
        np.savetxt(outputpath, distarray, delimiter=',')  # writing to file
    return distarray


def readdata(inputlist, datatype="distance", concatenate=True, separator=','):
    """
    Reading already written out data
    :param inputlist: List of datapaths to read in.
    :param datatype: The kind of data that is read in. Determines the column names
    :param concatenate: Concatenate provided data
    :param separator: Field separator
    :return: Returns dictionary of dataframes or dataframe
    """
    if datatype == "distance":
        columnnames = ['frame', 'dist']
    elif datatype == "rmsd":
        columnnames = ['frame', 'RMSD']
    elif datatype == "rmsf":
        # TODO: Complete reading in RMSF data.
        columnames = []
    else:
        columnnames = ['frame', 'data']

    if len(inputlist) == 1:
        df = pd.read_csv(inputlist[0],
                         sep='{}'.format(separator),
                         header=None,
                         usecols=[0, 2],
                         names=columnnames)
        return df
    elif len(inputlist) > 1:
        if concatenate is False:
            dfdict = {}
            for c in range(len(inputlist)):
                # TODO: datatype/ column names in not yet implemented..
                dfdict['df{}'.format(c)] = pd.read_csv(inputlist[c],
                                                       sep='{}'.format(separator),
                                                       header=None,
                                                       usecols=[0, 2],
                                                       names=['frame{}'.format(c), 'RMSD{}'.format(c)])
            return dfdict
        else:
            dfs = []  # list of dataframes
            for i in range(len(inputlist)):
                dfs.append(pd.read_csv(inputlist[i],
                                       sep='{}'.format(separator),
                                       header=None,
                                       usecols=[0, 2],
                                       names=['frame', 'dist']))
            concdf = pd.concat(dfs, ignore_index=True)
            return concdf
    else:
        print('Something went wrong..')
        return None


def rmsdplot(inputlist, simtime, labellist=None, plotstep=None, timeticks=50, figuresize=(32, 12), hidefig=None,
             style='default', outfile=''):
    """
    This function takes the following input and plots the results to screen (if active) and/or writes to file
    :param inputlist:   List of rmsd data file paths from MDAnalysis
    :param simtime:     Total or longest simulation in ns
    :param labellist:   List of data labels to be used for the legend
    :param plotstep:    Step size for data selection
    :param timeticks:   Ticks on the time axis every x ns
    :param figuresize:  Define figure size in cm
    :param hidefig:     If turned on, the Figure will not be displayed
    :param style:       Defining the style, default or seaborn
    :param outfile:     Path for the output file
    :return:            None
    """
    dfdict = {}
    figuresize = int(cm2inch(figuresize[0])), int(cm2inch(figuresize[1]))
    # Insane!!! A dictionary of pandas dataframes created by reading in the raw data!
    for c in range(len(inputlist)):
        dfdict['df{}'.format(c)] = pd.read_csv(inputlist[c],
                                               sep=',',
                                               header=None,
                                               usecols=[0, 2],
                                               names=['Time{}'.format(c), 'RMSD{}'.format(c)])

    df = pd.concat(dfdict.values(), axis=1)  # concatenating the imported dataframes
    df = df.T.drop_duplicates()  # getting rid of duplicated residue identifier, note that the dataframe is transposed for the work
    df = df.T  # transposing the dataframe back to its original orientation

    # preparing data for plotting ###
    # looking for the last entry of each row titled *Time* and getting the name of the row with the maximum value back
    timerow = max([(df[x].iloc[-1], x) for x in df if 'Time' in x])[1]  # list comprehension and if statement combined
    time = pd.Series(data=df[timerow], dtype='float64')
    f = simtime / round(max(time), -2)  # getting scaling factor for time, note that max is rounded to get rid of additional frames

    time = time.apply(lambda x: x * f)  # scaling time series
    if plotstep != None:
        time = time[::plotstep]
    ### plotting ###
    plt.style.use(style)
    plt.figure(figsize=figuresize)
    plt.ylabel('Root Mean Square Deviation [Å]')
    plt.xlabel('Time [ns]')
    #plt.yticks(np.arange(0, max(), 0.5))
    plt.xticks(np.arange(0, max(time) + 1, timeticks))
    # This command prints all RMSD values at once
    #[plt.plot(time[::plotstep], df[x][::plotstep], linewidth=1) for x in df if 'RMSD' in x]
    rmsdlist = [x for x in df if 'RMSD' in x]
    for x in rmsdlist:

        c = int(x.replace('RMSD', ''))  # getting the index number in order to get the labellist entry
        rmsvalues = pd.Series(data=df[x], dtype='float64')
        if plotstep is not None:
            rmsvalues = rmsvalues[::plotstep]
        plt.plot(time, rmsvalues, linewidth=1, label=labellist[c])
    plt.legend(loc='best')
    plt.tight_layout()

    # if a outfile path was provided, the figure will be saved
    if outfile != '':
        plt.savefig(fname=outfile, dpi=400, format='png', orientation='landscape')
    if hidefig is None:
        plt.show()


def rmsdplot_diffstepsize(inputlist, labellist=None, plotstep=None, timeticks=50, figuresize=(32, 12), hidefig=None,
                          style='default', outfile=''):
    """
    Special function for printing RMSD plots for simulations of different step size
    :param inputlist: List of lists with [rmsd data files from MDAnalysis, simulation time] format
    :param labellist: List of data labels to be used for the legend
    :param plotstep:  Step size for data selection
    :param timeticks: Ticks on the time axis every x ns
    :param figuresize: Define figure size in cm
    :param hidefig:   If turned on the Figure will not be displayed
    :param style:     Defining the style, default or seaborn
    :param outfile:   Path to output file
    :return: None
    """
    simtimelist= [x[1] for x in inputlist]
    inputlist=[x[0] for x in inputlist]
    dfdict = {}
    figuresize = int(cm2inch(figuresize[0])), int(cm2inch(figuresize[1]))
    for c in range(len(inputlist)):
        dfdict['df{}'.format(c)] = pd.read_csv(inputlist[c], sep=',', header=None, usecols=[0, 2],
                                               names=['Time{}'.format(c), 'RMSD{}'.format(c)])

    df = pd.concat(dfdict.values(), axis=1)  # concatenating the imported dataframes


    ### preparing data for plotting ###
    # looking for the last entry of each row titled *Time* and getting the name of the row with the maximum value back
    timerow = max([(df[x].iloc[-1], x) for x in df if 'Time' in x])[1]  # list comprehension and if statement combined
    time = pd.Series(data=df[timerow], dtype='float64')

    # Preparing plotting
    plt.style.use(style)
    plt.figure(figsize=figuresize)
    plt.ylabel('Root Mean Square Deviation [Å]')
    plt.xlabel('Time [ns]')
    plt.xticks(np.arange(0, int(time.max()), timeticks))
    plt.xticks(np.arange(0, max(time) + 1, timeticks))


    # Time rescaling
    for c in range(len(inputlist)):
        simtime = simtimelist[c]  # getting simulation time
        f = simtime / round(df['Time{}'.format(c)].max(skipna=True), -2)  # getting scaling factor
        time = pd.Series(data=df['Time{}'.format(c)], dtype='float64')  # extracting time into series for plotting
        time = time.apply(lambda x: x * f)  # actual time rescaling
        rmsvalues = pd.Series(data=df['RMSD{}'.format(c)], dtype='float64')
        if plotstep is not None:
            rmsvalues = rmsvalues[::plotstep]
        plt.plot(time, rmsvalues, label=labellist[c])

    plt.legend(loc='best')
    plt.tight_layout()

    # if a outfile path was provided, the figure will be saved
    if outfile != '':
        plt.savefig(fname=outfile, dpi=400, format='png', orientation='landscape')
    if hidefig is None:
        plt.show()


def rmsfplot(inputlist, labellist=None, residueticks=50, figuresize=(32, 12), ymax=None, outfile=''):
    """
    Plotting results of RMSF calculation
    :param inputlist:       List of rmsf data files from MDAnalysis
    :param labellist:       List of data labels to be used for the legend
    :param residueticks:    Ticks on the x axis
    :param figuresize:      Figure size in cm
    :param ymax:            # maximum y value in Å. Add + 0.5 to desired value
    :param outfile:         Path to output file PNG
    :return:                None
    """


    dfdict = {}
    N = len(inputlist)  # needed for statistical calculations
    figuresize = int(cm2inch(figuresize[0])), int(cm2inch(figuresize[1]))
    # Insane!!! A dictionary of pandas dataframes created by reading in the raw data!
    for c in range(len(inputlist)):
        dfdict['df{}'.format(c)] = pd.read_csv(inputlist[c], sep=',', header=None, names=['res', 'resid', 'rmsf{}'.format(c)])

    df = pd.concat(dfdict.values(), axis=1)  # concatenating the imported dataframes
    df = df.T.drop_duplicates()  # getting rid of duplicated residue identifier, note that the dataframe is transposed for the work
    df = df.T  # transposing the dataframe back to its original orientation
    df['resname'] = df['res'].map(str) + df['resid'].map(str)  # concatenating residue name and saving it in the 'resname' column

    # Commencing statistical sorcery aka mean and standard error calculation

    ### plotting ###



    # plotting settings
    plt.figure(figsize=figuresize)
    plt.ylabel('Root Mean Square Fluctuation [Å]')
    plt.xlabel('Residue number')

    for c in range(len(inputlist)):
        x = pd.Series(data=df['resid'], dtype='float64')
        rms = pd.Series(data=df['rmsf{}'.format(c)], dtype='float64')
        if labellist is None:
            plt.plot(x, rms, linewidth=1.5)
        else:
            plt.plot(x, rms, linewidth=1.5, label=labellist[c])


    if ymax is not None:
        plt.yticks(np.arange(0, ymax, 0.5))
    plt.xticks(np.arange(0, max(x) + 1, residueticks))  # rotation='vertical')
    plt.legend(loc='best')
    plt.tight_layout()
    if outfile != '':
        plt.savefig(fname=outfile, dpi=400, format='png', orientation='landscape')
    plt.show()


def rmsfplotSE(inputlist, residueticks=50, figuresize=(32, 12), ymax=None, outfile=''):
    """
    Plotting results of RMSF calculation using mean +- standard error
    :param inputlist:       List of rmsf data files from MDAnalysis
    :param residueticks:    Ticks on the x axis
    :param figuresize:      Figure size in cm
    :param ymax:            # maximum y value in Å. Add + 0.5 to desired value
    :param outfile:         Path to output file PNG
    :return:                None
    """
    dfdict = {}
    N = len(inputlist)  # needed for statistical calculations
    figuresize = int(cm2inch(figuresize[0])), int(cm2inch(figuresize[1]))
    # Insane!!! A dictionary of pandas dataframes created by reading in the raw data!
    for c in range(len(inputlist)):
        dfdict['df{}'.format(c)] = pd.read_csv(inputlist[c], sep=',', header=None,
                                               names=['res', 'resid', 'rmsf{}'.format(c)])
    df = pd.concat(dfdict.values(), axis=1)  # concatenating the imported dataframes
    df = df.T.drop_duplicates()  # getting rid of duplicated residue identifier, note that the dataframe is transposed for the work
    df = df.T  # transposing the dataframe back to its original orientation
    df['resname'] = df['res'].map(str) + df['resid'].map(
        str)  # concatenating residue name and saving it in the 'resname' column
    # Commencing statistical sorcery aka mean and standard error calculation
    df['m'] = df[[i for i in df.columns if 'rmsf' in i]].mean(axis=1)  # calculating the mean rmsf in column 'm'
    df['SEu'] = df['m'].map(float) + stderr(df['m'].map(float), N)  # calculating Standard error upper boundary
    print(df)
    df['SEl'] = df['m'].map(float) - stderr(df['m'].map(float), N)  # calculating Standard error lower boundary
    ### plotting ###
    x = pd.Series(data=df['resid'], dtype='float64')
    m = pd.Series(data=df['m'], dtype='float64')
    SEu = pd.Series(data=df['SEu'], dtype='float64')
    SEl = pd.Series(data=df['SEl'], dtype='float64')
    # plotting settings
    plt.figure(figsize=figuresize)
    plt.plot(x, m, color='black', linewidth=1.5, label='Mean of all MD simulation replicas')
    plt.plot(x, SEu, color='grey', linewidth=1, label='Standard error of the mean')
    plt.plot(x, SEl, color='grey', linewidth=1, label='_SEl')
    plt.fill_between(x, SEl, SEu, color='grey', alpha='0.5')
    plt.ylabel('Root Mean Square Fluctuation [Å]')
    plt.xlabel('Residue number')
    if ymax is None:
        ymax = max(SEu)
    plt.yticks(np.arange(0, ymax, 0.5))
    plt.xticks(np.arange(0, max(x) + 1, residueticks))  # rotation='vertical')
    plt.legend( loc='best')
    plt.tight_layout()
    if outfile != '':
        plt.savefig(fname=outfile, dpi=400, format='png', orientation='landscape')
    plt.show()


def disthistplot(inputlist, style='seaborn', labellist=None, outfile=None, ylabel=None, xlabel=None):
    """
    Plotting distance distribution between two atoms for multiple trajectories
    :param inputlist:
    :param style:
    :param labellist:
    :param outfile:
    :param ylabel:
    :param xlabel:
    :return:
    """
    sns.set_context("paper")
    # creating a dictionary of dataframes
    d = {}
    for i in range(len(inputlist)):
        d['dist{}'.format(i)] = pd.read_csv(inputlist[i], sep='\t', header=None, names=['frame', 'dist{}'.format(i)])

    plt.style.use(style)  # switiching the representation style
    if labellist is None:
        for key, i in zip(d.keys(), range(len(d.keys()))):
            sns.distplot(d[key][key], hist=False, label='replica{}'.format(str(int(i)+1)))  # plotting each dataframe
    else:
        for key, i in zip(d.keys(), range(len(d.keys()))):
            sns.distplot(d[key][key], hist=False, label=labellist[i])  # plotting each dataframe

    # plotting the distribution in one diagram
    if ylabel is None:
        plt.ylabel('Kernel density function')  # default yaxis label
    else:
        plt.ylabel(ylabel)  # using provided yaxis label
    if xlabel is None:
        plt.xlabel('Distance [Å]')  # default xaxis label
    else:
        plt.xlabel(xlabel)  # using provided xaxis label
    plt.tight_layout()
    if outfile is not None:
        plt.savefig(fname=outfile, dpi=400, format='png', orientation='landscape')  # saving the figure, if path given
    plt.show()


def rmsdligsubplot(inputlist, labellist=None, plotstep=None, timeticks=50, figuresize=(32, 12), hidefig=None,
    style='seaborn', outfile=None):
    """
    This functionality accepts ligand RMSD data obtained with MDAnalysis and plots the data over time and as histogram (2 subplots)
    :param inputlist:   List of lists of input in the format of [['file1', simulationtime1], ['file2', simulationtime2]]
    :param labellist:   Labels to be used for data
    :param plotstep:    Step size for data selection
    :param timeticks:   Ticks on the time axis every x ns
    :param figuresize:  Define figure size in cm
    :param hidefig:     If turned on the Figure will not be displayed
    :param style:       Defining the style, default or seaborn
    :param outfile:     Path for the output file
    :return:
    """
    simtimelist= [x[1] for x in inputlist]
    inputlist=[x[0] for x in inputlist]
    dfdict = {}
    figuresize = int(cm2inch(figuresize[0])), int(cm2inch(figuresize[1]))

    # Creating a dictionary to read rmsd values as individual dataframes
    for c in range(len(inputlist)):
        dfdict['df{}'.format(c)] = pd.read_csv(inputlist[c], sep=',', header=None, usecols=[0, 2],
                                               names=['Time{}'.format(c), 'RMSD{}'.format(c)])
    df = pd.concat(dfdict.values(), axis=1)  # concatenating the imported dataframes

    ### preparing data for plotting ###
    # looking for the last entry of each row titled *Time* and getting the name of the row with the maximum value back
    timerow = max([(df[x].iloc[-1], x) for x in df if 'Time' in x])[1]  # list comprehension and if statement combined
    time = pd.Series(data=df[timerow], dtype='float64')

    # Preparing plotting
    plt.style.use(style)  # setting style for plotting
    plt.figure(figsize=figuresize)
    plt.ylabel('Root Mean Square Deviation [Å]')
    plt.xticks(np.arange(0, int(time.max()), timeticks))
    plt.xticks(np.arange(0, max(time) + 1, timeticks))
    plt.subplot2grid((1, 3), (0, 0), colspan=2)  # Initializing subplot with two horizontal subfigures

    # Time rescaling
    for c in range(len(inputlist)):
        simtime = simtimelist[c]  # getting simulation time
        f = simtime / round(df['Time{}'.format(c)].max(skipna=True), -2)  # getting scaling factor
        time = pd.Series(data=df['Time{}'.format(c)], dtype='float64')  # extracting time into series for plotting
        time = time.apply(lambda x: x * f)  # actual time rescaling
        rmsvalues = pd.Series(data=df['RMSD{}'.format(c)], dtype='float64')
        if plotstep is not None:
            rmsvalues = rmsvalues[::plotstep]
        if labellist is not None:
            plt.plot(time, rmsvalues, label=labellist[c])
        else:
            plt.plot(time, rmsvalues, label='Replica  {}'.format(c+1))
        plt.xlabel('Time [ns]')

    plt.legend(loc='best')
    plt.ylabel('Root Mean Square Deviation [Å]')
    plt.tight_layout()


    # reading data for kernel dist /histogram in the form of dictionary of dataframes
    d = {}
    for i in range(len(inputlist)):
        d['RMSD{}'.format(i)] = pd.read_csv(inputlist[i], sep=',', header=None, usecols=[0, 2], names=['frame', 'RMSD{}'.format(i)])

    # preparting plotting histogram
    plt.style.use(style)
    plt.subplot2grid((1, 3), (0, 2), rowspan=1)
    for key, i in zip(d.keys(), range(len(d.keys()))):
        sns.distplot(d[key][key], hist=False, vertical=True)
        plt.xlabel('Kernel density function')

    # plotting the distribution in one diagram
    plt.ylabel('')
    plt.xlabel('Kernel density function')
    plt.tight_layout()

    # if a outfile path was provided, the figure will be saved
    if outfile != '':
        plt.savefig(fname=outfile, dpi=400, format='png', orientation='landscape')

    if hidefig is None:
        plt.show()


def histplot(inputlist, style='seaborn', context=None, labellist=None, outfile=None, ylabel=None, xlabel=None, concatenate=False):
    """

    :param inputlist:   List of files containing distances obtained by MDAnalysis
    :param style:       Defining the style, default or seaborn
    :param labellist:   Labels to be used for data
    :param outfile:     Path for the output file
    :param ylabel:      Label to be used for the y axis; default is 'Kernel density function'
    :param xlabel:      Label to be used for the x axis; default is 'Distance [Å]'
    :param concatenate: Concatenate input list and plot resulting overall distribution.
    :return:            None
    """
    plt.style.use(style)  # switching the representation style
    if context is not None:
        sns.set_context("{}".format(context))

    if concatenate is False:
        # creating a dictionary of dataframes
        d = {}
        for i in range(len(inputlist)):
            d['dist{}'.format(i)] = pd.read_csv(inputlist[i], sep=',', header=None, usecols=[0, 2],
                                                names=['frame', 'dist{}'.format(i)])
        if labellist is None:
            for key, i in zip(d.keys(), range(len(d.keys()))):
                sns.distplot(d[key][key], hist=False, label='replica{}'.format(str(int(i) + 1)))  # plotting each dataframe
        else:
            for key, i in zip(d.keys(), range(len(d.keys()))):
                sns.distplot(d[key][key], hist=False, label=labellist[i])  # plotting each dataframe
    # Plotting concatenated data
    else:
        dfs = []  # list of dataframes
        for i in range(len(inputlist)):
            dfs.append(pd.read_csv(inputlist[i], sep=',', header=None, usecols=[0, 2], names=['frame', 'dist']))
        concdf = pd.concat(dfs, ignore_index=True)
        sns.distplot(concdf['dist'], hist=False, label=labellist)

        # plotting the distribution in one diagram
    if ylabel is None:
        plt.ylabel('Kernel density function')  # default yaxis label
    else:
        plt.ylabel(ylabel)  # using provided yaxis label
    if xlabel is None:
        plt.xlabel('Distance [Å]')  # default xaxis label
    else:
        plt.xlabel(xlabel)  # using provided xaxis label
    plt.tight_layout()
    if outfile is not None:
        plt.savefig(fname=outfile, dpi=400, format='png', orientation='landscape')  # saving the figure, if path given
    plt.show()


### Plotting distance between two atoms over time for multiple trajectories ###

### Hbonding analysis

### run mdanalysis ###
