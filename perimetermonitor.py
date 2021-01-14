'''Perimeter monitoring for MD simulations using MDAnalysis'''
import glob
import MDAnalysis
from MDAnalysis.analysis import distances
import numpy as np

# Please change the workdir here..
workdir = 'REPLACEME'

top = glob.glob("{0}/**/{1}".format(workdir.rstrip("/"), '*-in-noh2o.pdb'), recursive=True)[0]
traj = glob.glob("{0}/**/{1}".format(workdir.rstrip("/"), "*.dcd"), recursive=True)[0]
outputpath = "{0}/{1}".format(workdir.rstrip("/"), "perimeter.csv")

# For manual mode..
#top = '/mdspace/mstahnke/2007_2YDV_Na_MD/2YDV_Na_MD-in-noh2o.pdb'
#traj = '/mdspace/mstahnke/2007_2YDV_Na_MD/2YDV_Na_MD_trj/allframes-noh2o.dcd'
#outputpath = '/mdspace/mstahnke/2007_2YDV_Na_MD/perimeter.csv'


u = MDAnalysis.Universe(top, traj)  # Creating MDAnalysis universe aka loading in MD data
print('Created MDUniverse successfully.')

# creating atom selections for later use
s1 = u.select_atoms("resid 246 and resname TRP and name CA")
s2 = u.select_atoms("resid 91 and resname SER and name CA")
s3 = u.select_atoms("resid 52 and resname ASP and name CA")
s4 = u.select_atoms("resid 280 and resname ASN and name CA")
lig = u.select_atoms("resid 400")
sodium = u.select_atoms("resid 2402")

# getting distances
distdict = {}
distpairs = [(s1, s2), (s2, s3), (s3, s4), (s4, s1)]
step = -1
perilist = []
print('Starting perimeter calculations.')
for ts in u.trajectory:
    step += 1
    peri = 0.0
    for n, distpair in enumerate(distpairs):
        a = distpair[0]
        b = distpair[1]
    # distance array cropped to the distance between selections rounded to two decimal digits
        dist = round(MDAnalysis.analysis.distances.dist(a, b)[2, 0], 2)
        peri += dist
    perilist.append([step, u.trajectory.time, peri])
perilist = np.array(perilist)

print('Finished perimeter calculations and writing out file to \"{}\"'.format(outputpath))
# adding distances
# writing out
if outputpath is not None:
    np.savetxt(outputpath, perilist, delimiter=',')  # writing to file
