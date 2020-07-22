#!/usr/bin/python3

ppafm = '/home/matyas/Scripts/Python/pypath/ProbeParticleModel'

import numpy as np
import re
import sys

sys.path.append(ppafm)
import pyProbeParticle.GridUtils as GU

arg = sys.argv
art = False

for ar in arg:
    if ar == '-art':
        art = True
        arg.remove('-art')

try:
    InputFile = sys.argv[1]
except:
    print('Please select Input File')
    quit()
    
    
ListOptions = ['point_grid', 'eval_height']
VarOptions = ['geo_input', 'cutoff', 'long_range', 'cluster']
K = 14.3996


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
    CheckExistence(['geo_input', 'point_grid','eval_height'])#Check existance on nonoptional variables

    clus = False
    longr = False
    
    if 'cluster' in [x for x in globals()]:
        if cluster == 1:
            clus = True
        else:
            try:
                clus = True if cluster else False
            except:
                pass
                
    if 'long_range' in [x for x in globals()]:
        if long_range == 1:
            longr = True
        else:
            try:
                longr = True if long_range else False
            except:
                pass
                
    if not clus:
        CheckExistence(['cutoff'])
    else:
        global cutoff
        cutoff = 1000000 
        
    return clus, longr


def ReadGeo():

    fil = open(geo_input, 'r')
    doc = [x.split('#')[0] for x in fil.read().splitlines()]
    fil.close()

    Lattice = np.array([[float(x) for x in vec.split()] for vec in doc[0].split(',')])   
    dipole_inp = [x.split(',') for x in doc[1:] if x != '']
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
    print (1)
    return PointGrid


def ArtSample(vector, origin, N):
       #creates N vectors in direcection of the initial vector with length uniformly distributed                                     
    # from origin(non 0) to length of initial vector                                                                              

    coord = [np.linspace(origin[i], x, N) for i, x in enumerate(vector)]
    vector_grid = [np.array([x, coord[1][i], coord[2][i]]) for i, x in enumerate(coord[0])]
    return vector_grid


def Art(Lattice):
    #creates a grid o points over space where potential will be evaluated.                                                        
    #Those will be in a cell base of the unit cell and user defined height                                                        

    shape_origin = np.array([0,0,eval_height[0]])
    eval_shape = [x + shape_origin for x in Lattice]
    eval_shape.append(np.array([0,0,eval_height[1]]))

    agrid = ArtSample(eval_shape[0], shape_origin, point_grid[0])
    bgrid = ArtSample(eval_shape[1], shape_origin, point_grid[1])
    cgrid = ArtSample(eval_shape[2], shape_origin, point_grid[2])

    PointGrid = np.zeros((point_grid[2], point_grid[1], point_grid[0], 3))


    for z, c in enumerate(cgrid):
        for y, b in enumerate(bgrid):
            for x, a in enumerate(agrid):
                PointGrid[z, y, x] = a + b + c
    
    return PointGrid


def FindMaxLatice(Lattice):
    maxLat = [int(cutoff/np.sqrt(Lattice[i].dot(Lattice[i]))) + 3 for i in range(2)]

    return maxLat

def getELongRange(Dipoles, Lattice):

    if longr:
        
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
    
    if clus:
        print(3)
        maxLat = [0,0]
        Lattice = np.zeros([2,3])
        E_long_range = np.zeros([3])

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


def WriteOutput(Potential, Lattice):

    lat_mat = np.array([[0,0, eval_height[0]],
                       Lattice[0],
                       Lattice[1],
                       [0,0, eval_height[1] - eval_height[0]]])
    
    GU.saveXSF('dipole_potential.xsf', Potential, lat_mat)



### FLOW STARTS HERE ###


ReadInput()
clus, longr = goThroughVar()
Dipoles, Lattice = ReadGeo()

if art:
   PointGrid = Art(Lattice)
else:
    PointGrid = CreateGrid(Lattice)
    
Potential = CalculateGrid(Dipoles, PointGrid, Lattice)
WriteOutput(Potential, Lattice)

