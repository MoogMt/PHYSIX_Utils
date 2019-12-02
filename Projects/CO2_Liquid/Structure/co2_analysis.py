# Loading all important libraries

import numpy as np
import scipy as scp
import dscribe as ds
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import time

# folder stuff

home_folder=str("/home/moogmt");
co2_folder=home_folder+str("/CO2/CO2_AIMD/");
co2_folder="/media/moogmt/Stock/Mathieu/CO2/AIMD/Liquid/PBE-MT/"

# Temperature and volumes

V=8.82;
T=3000;

# final folder

data_folder=co2_folder+str(V)+"/"+str(T)+"K/";

# final .xyz

filepath=data_folder+str("TRAJEC_wrapped.xyz");

# Useful functions

def getNbSteps( filepath ):
    # Values to return
    step_number=0
    nb_atoms=0
    # Reading file
    fp=open( filepath, "r" )
    read=fp.readline()
    if read != "" :
        nb_atoms=int(((read.rstrip("\n")).split())[0])
    while ( read != "" ) :
        step_number += 1
        read=fp.readline()
    fp.close()
    return nb_atoms, int(step_number/(nb_atoms+2))

def getTypes( filepath, nb_atoms ):
    # Initialize file
    types_names = ["" for x in range(nb_atoms)]
    fp=open( filepath, "r" )
    fp.readline()
    fp.readline()
    for atom in range(nb_atoms):
        types_names[atom]=((fp.readline().rstrip("\n")).split())[0]
    fp.close()
    types=[ types_names[1] ]
    types_vector=[0]
    for atom in range(nb_atoms):
        found = False
        for i in range(len(types)):
            if types_names[atom] == types[i] :
                types_vector=np.append(types_vector,i)
                found = True
                break
        if not(found):
            types=np.append(types,types_names[atom])
            types_vector=np.append(types_vector,i)
    return types, types_vector

def readXYZ( filepath ):
    nb_atoms, nb_steps = getNbSteps( filepath )
    types, types_vector = getTypes( filepath, nb_atoms )
    positions=np.zeros(( nb_steps, nb_atoms, 3 ))
    fp=open( filepath, "r" )
    for step in range( nb_steps ):
        # Reading two comments lines
        fp.readline()
        fp.readline()
        for atom in range( nb_atoms ):
            line=(fp.readline()).rstrip("\n").split()
            # Reading two comments lines
            for i in range(3):
                positions[step,atom,i] = float(line[i+1])
    fp.close()
    return positions, types, types_vector

def cellOrthoDistance( x1, x2, a):
    dx=x1-x2
    if dx > a*0.5:
        dx -= a
    if dx < -a*0.5:
        dx += a
    return dx*dx

def distance( position1, position2, cell_param ):
    dist=0
    for i in range(3):
        dist += cellOrthoDistance( position1[i], position2[i], cell_param[i] )
    return np.sqrt(dist)

def distanceMatrix( positions, nb_atoms, cell_lengths ):
    distance_matrix=np.zeros((nb_atoms,nb_atoms))
    for i in range(nb_atoms):
        for j in range(i+1,nb_atoms):
            distance_matrix[i,j]=distance( positions[i,:], positions[j,:], cell_lengths )
            distance_matrix[j,i]=distance_matrix[i,j]
    return distance_matrix

def writeData( filepath, data ):
    fp=open( filepath, "w" )
    for i in range(np.shape(data)[0]):
        for j in range(np.shape(data)[1]):
            fp.write(str(data[i,j])+" ")
        fp.write("\n")
    fp.close()
    return

def writeCubeFile( filepath, vol_data, cell_param ):
    fp=open( filepath, "w" )
    fp.write("Data - Volume 1 \n")
    fp.write("Data 2 - X Y Z \n")
    # Only one atom at origin (useless, virtual, we care about the vol data here)
    fp.write("1 0.0 0.0 0.0 \n")
    # Voxel and cell
    voxels=np.shape(vol_data)
    for i in range(3):
        fp.write(str(voxels[i])+" ")
        for j in range(i):
            fp.write("0.0 ")
        fp.write(str(cell_param[i])+" ")
        for j in range(i+1,3):
            fp.write("0.0 ")
        fp.write("\n")
    # Fictitious atoms
    fp.write("1 0.0 0.0 0.0\n")
    # Volumetric data
    for i in range(voxels[0]):
        for j in range(voxels[1]):
            for k in range(voxels[2]):
                fp.write(str(vol_data[i,j,k])+" ")
                if k % 6 == 5 :
                    fp.write("\n")
            fp.write("\n")
    # Closing file
    fp.close()
    return

def computeDensity( data, nb_points ):
    density=np.zeros(( nb_points ))
    mins=np.zeros((len(nb_points)))
    maxs=np.zeros((len(nb_points)))
    deltas=np.zeros((len(nb_points)))
    for i in range(len(nb_points)):
        mins[i] = np.min(data[:,i])
        maxs[i] = np.max(data[:,i])
        deltas[i] = (maxs[i]-mins[i])/nb_points[i]
    nb_data=len(data[:,0])
    for i in range(nb_data):
        indexs=np.zeros((3))
        for j in range(3):
            indexs[j]=(data[i,j]-mins[j])/deltas[j]-1
        density[int(indexs[0]),int(indexs[1]),int(indexs[2])] += 1
    density /= np.max(density)
    return density, mins, maxs, deltas

    # Reading file
    positions,types,types_vector=readXYZ(filepath)
    nb_steps,nb_atoms,dim = np.shape(positions)

    # Computing Data

max_neigh=4
nbC=32
nb_steps_max=5000
nb_data=nb_steps_max*nbC
data=np.zeros((nb_data,max_neigh))

cell=[V,V,V]

count=0
for step in range(nb_steps_max):
    distance_matrix=distanceMatrix(positions[step,:,:],nb_atoms,cell)
    for atom in range(nbC):
        indexs=np.argsort(distance_matrix[atom,:])
        distance_matrix[atom,:]=distance_matrix[atom,indexs]
        data[count,:]=distance_matrix[atom,1:max_neigh+1]
        count += 1

        # Writting data to file
file_out=data_folder+str("data.dat");
writeData(file_out,data)

nb_boxes=[100,100,100]
# Compute the point density
density,mins,maxs,deltas=computeDensity(data[:,1:4],nb_boxes)

        # Writing to cubefile
file_out=data_folder+str("data.cube");
writeCubeFile( file_out, density, cell )

file_read="/home/moogmt/Downloads/positions_coupledCluster.npy"

data=np.load(file_read)
np.shape(data)

nb_steps, nb_atoms, dim = np.shape(data)
for i in range(nb_steps):
    fp=[]
    if i < 25000:
        fp=open("/home/moogmt/H2O_positions/sec_1/"+str(i)+"_pos.dat","w")
    else:
        fp=open("/home/moogmt/H2O_positions/sec_2/"+str(i)+"_pos.dat","w")
    for j in range(nb_atoms):
        if j < 2 :
            fp.write("O ")
        else:
            fp.write("H ")
        for k in range(dim):
            fp.write(str(data[i,j,k])+" ")
        fp.write("\n")
    fp.close()
