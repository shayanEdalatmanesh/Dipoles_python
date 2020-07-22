#!/usr/bin/python3

ppafm = '/home/matyas/Scripts/Python/pypath/ProbeParticleModel'

import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import os

sys.path.append(ppafm)
import pyProbeParticle.GridUtils as GU


class Dipole():

    def __init__(self, momentum, possition):

        self.momentum = momentum
        self.possition = possition

        


def ReadInput():
    #Reads Variables from input files and assigns them Scarry how easy it is to do this in python                                      

    fil = open(InputFile, "r")
    doc = fil.read()
    fil.close()

    #spilt by ; or \n then ignore coments (#), erases coments and split into variable name and value
    doc = doc.replace("\n",";").replace(" ", "")
    doc = [x.split('#')[0].split("=") for x in doc.split(';')]

    #erase empty
    doc = [x for x in doc if x != ['']]

    #Assign
    for opt in doc:
        if opt[0] in VarOptions:

            try:
                exec("global %s; %s = %s" % (opt[0], opt[0], opt[1]))
            except:
                exec('global %s; %s ="%s"' % (opt[0], opt[0], opt[1]))
        elif opt[0] in ListOptions:
            try:
                exec("global %s; %s = [%s]" % (opt[0], opt[0], opt[1]))
            except:
                val = "'" + opt[1].replace(',',"','") + "'"
                exec("global %s; %s = [%s]" % (opt[0], opt[0], val))
        else:
            print("Error. Option %s in %s is unknown" % (opt[0], InputFile))
            quit()


def CheckExistence(variables):
    #takes in varible name and cheks if it exists                                                                                      
    glob = [x for x in globals()]#globals is a dictonary an glob is a list. in does not seem to work for dictonaries
    
    for var in variables:
        if var not in glob:
            print("Error. %s not defined" % var)
            quit()


def goThroughVar():
    #changes format and existance of some variables
    CheckExistence(['geo_input', 'cutoff','point_grid','eval_height', 'long_range'])#Check existance on nonoptional variables

    
    
def ReadGeo():

    fil = open(geo_input, 'r')
    doc = [x.split('#')[0] for x in fil.read().splitlines()]
    fil.close()

    dipole_inp = [x.split(',') for x in doc[1:] if x != '']

    Lattice = np.array([[float(x) for x in vec.split()] for vec in doc[0].split(',')])

    Dipoles = []

    for dip in dipole_inp:
        pos = np.array([float(x) for x in dip[0].split()])
        moment =  np.array([float(x) for x in dip[1].split()])
        Dipoles.append(Dipole(moment, pos))
    
    return Dipoles, Lattice
    


def sample_vector(vector, N):
    #creates N vectors in direcection of the initial vector with length uniformly distributed
    # from origin(non 0) to length of initial vector
    
    coord = [np.linspace(0, x, N) for i, x in enumerate(vector)]
    vector_grid = [np.array([x, coord[1][i], coord[2][i]]) for i, x in enumerate(coord[0])]
    
    return vector_grid
    

def CreateGrid(Lattice):
    #creates a grid o points over space where potential will be evaluated.
    #Those will be in a cell base of the unit cell and user defined height
   
    shape_origin = np.array([0,0,eval_height[0]])
    eval_shape =  [x for x in Lattice]
    eval_shape.append(np.array([0,0,eval_height[1]-eval_height[0]]))
   
    agrid = sample_vector(eval_shape[0], point_grid[0])
    bgrid = sample_vector(eval_shape[1], point_grid[1])
    cgrid = sample_vector(eval_shape[2], point_grid[2])
    
    PointGrid = np.zeros((point_grid[2], point_grid[1], point_grid[0], 3))
    
    
    for z, c in enumerate(cgrid):
        for y, b in enumerate(bgrid):
            for x, a in enumerate(agrid):
                
                PointGrid[z, y, x] = shape_origin + a + b + c

    return PointGrid


def FindMaxLatice(Lattice):
    maxLat = [int(cutoff/np.sqrt(Lattice[i].dot(Lattice[i]))) + 3 for i in range(2)]

    return maxLat

def getELongRange(Dipoles, Lattice):

    if long_range:
        
        cross = np.cross(Lattice[0], Lattice[1])
        cell_area = np.sqrt(cross.dot(cross))
        Net_moment = 0
        
        for Dip in Dipoles:
            Net_moment += Dip.momentum 

        
        E_long_range = 4*np.pi*Net_moment/(3*cutoff*cell_area)
            
    else:
        E_long_range = np.zeros([1,3])
    
    return E_long_range


def CalculateGrid(Dipoles, PointGrid, Lattice):
    #Calculates the potential at each point. NOT TESTED I have no input.
    
    maxLat = FindMaxLatice(Lattice)#mamximum number of unit cell in one direction from origin being evaluated [+-x,=-y]

    cutoff_cube = abs(cutoff**3)
    
    Potential = np.zeros((point_grid[2], point_grid[1], point_grid[0]))
    
    E_long_range = getELongRange(Dipoles, Lattice)


    for z, Z in enumerate(PointGrid):
        for y, Y in enumerate(Z):
            for x, pos in enumerate(Y):

                Potential[z, y, x] += np.dot(E_long_range, pos)
                #I will send this equation on paper on Tuesday/Wednesday

                for Rx in range(-maxLat[0], maxLat[0]+1):
                    for Ry in range(-maxLat[1], maxLat[1]+1):
                        lat_dist = pos - Rx*Lattice[0] - Ry*Lattice[1]
                        
                        for Dip in Dipoles:

                            rel_pos = lat_dist - Dip.possition
                            
                            rel_pos_cube = (np.sqrt(rel_pos.dot(rel_pos)))**3
                            
                            if rel_pos_cube < cutoff_cube:
                                Potential[z, y, x] += np.dot(Dip.momentum, rel_pos)/rel_pos_cube

                Potential[z, y, x]*=K


        # I am hungry

    return Potential

def CalculateGrid2(Dipoles, PointGrid):

    Potential = np.zeros((point_grid[2], point_grid[1], point_grid[0]))
    for z, Z in enumerate(PointGrid):
        for y, Y in enumerate(Z):
            for x, pos in enumerate(Y):
                
                for Dip in Dipoles:
                    
                    rel_pos = pos - Dip.possition
                    rel_pos_cube = (np.sqrt(rel_pos.dot(rel_pos)))**3
                    
                    Potential[z, y, x] += np.dot(Dip.momentum, rel_pos)/rel_pos_cube
                    
                Potential[z, y, x]*=K
    return Potential

def WriteOutput(Potential, Lattice):

    lat_mat = np.array([[0,0, eval_height[0]],
                       Lattice[0],
                       Lattice[1],
                       [0,0, eval_height[1] - eval_height[0]]])
    
    GU.saveXSF('dipole_potential.xsf', Potential, lat_mat)

    
def DrawDipoles(Dipoles, Lattice, delta=0.5, axis='on'):

    size = max(Lattice[0][0],Lattice[1][1])
    plt.xlim(0, size)
    plt.ylim(0, size)

    plt.axis('off')
    
    for Dip in Dipoles:
        pos = Dip.possition
        mom = Dip.momentum*delta
        plt.arrow(pos[0],pos[1],mom[0],mom[1], color='black', head_width=.5)
        
    plt.savefig('DipolesNoAxes.png')


### FLOW STARTS HERE ###

if __name__ == '__main__':

    try:
        InputFile = sys.argv[1]
    except:
        print('Please select Input File')
        quit()
    
    
    ListOptions = ['point_grid', 'eval_height']
    VarOptions = ['geo_input', 'cutoff', 'long_range', 'cluster', 'plotDip']
    K = 14.3996

    ReadInput()
    goThroughVar()

    Dipoles, Lattice = ReadGeo()
    PointGrid = CreateGrid(Lattice)

    try:
        if cluster:
            Potential = CalculateGrid2(Dipoles, PointGrid, Lattice)
        else:
            Potential = CalculateGrid(Dipoles, PointGrid, Lattice)
    except:
        Potential = CalculateGrid(Dipoles, PointGrid, Lattice)

    WriteOutput(Potential, Lattice)

    
    if plotDip:
        DrawDipoles(Dipoles, Lattice)
        os.sytem('convert Dipoles.png -transparent white Dipoles.png')

