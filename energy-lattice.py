# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 11:28:55 2020

@author: Maywell2019
"""

import re
import sys
import os
import shutil
import numpy as np

path = os.getcwd()
f_list = os.listdir(path)
for i in f_list:
    if os.path.splitext(i)[1] == '.castep':
        seedname = os.path.splitext(i)[0]

amplitude = 0.1
number = 11

# regular expression to match the whole of the final cell from a .castep file
dotcastep_latt_RE = re.compile("""\sL?BFGS\s*:\sFinal\sConfiguration:\s*\n
                                =+\s*\n\s*\n\s+\-+\s*\n\s+Unit\sCell\s*\n\s+\-+\s*\n
                                \s+Real\sLattice\(A\)\s+Reciprocal\sLattice\(1/A\)\s*\n
                                \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
                                \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
                                \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n""",re.VERBOSE)
dotcell_lattice_start_RE = re.compile("^\s*%BLOCK\s+LATTICE_(?:CART|ABC)",re.VERBOSE)
dotcell_lattice_end_RE = re.compile("^\s*%ENDBLOCK\s+LATTICE_(?:CART|ABC)",re.VERBOSE)
volume_RE = re.compile("^\s*Current\s+cell\s+volume\s+=\s+([0-9\.]+)\s*A\*\*3")
density_RE = re.compile("^\s*=\s+([0-9\.]+)\s*g\/cm\^3+\n")
energy_RE = re.compile("^\s*BFGS:\s+Final\s+Enthalpy\s+=\s+([\-?0-9\.]+E[\+\-0-9]+)\s*eV")

inputfile = open(seedname+".castep","r")
result = dotcastep_latt_RE.findall(inputfile.read())[-1]
lattice = []
lattice.append([float(result[0]),float(result[1]),float(result[2])])
lattice.append([float(result[6]),float(result[7]),float(result[8])])
lattice.append([float(result[12]),float(result[13]),float(result[14])])
lattice = np.array(lattice)
in_lattice = False
inputfile.close()


########################### obtain lattice change ########################
change_list = []
n = 0
while n < number:
    change = amplitude * 2 / (number - 1) * amplitude * n - amplitude + 1
    change_list.append(change)
    n = n + 1

for change in change_list:
    new_lattice = np.dot(lattice,change)
    inputfile_0 = open(seedname+".cell", "r")
    outputfile = open(seedname+"_"+change+".cell","w")
    for line in inputfile_0:
######### write lattice parameters ##############################        
        if (re.search(dotcell_lattice_end_RE,line) and in_lattice):
            in_lattice = False
        elif (re.search(dotcell_lattice_start_RE,line) and not in_lattice):
            outputfile.write("%BLOCK LATTICE_CART\n")
            outputfile.write("       " + '%0.15f'%new_lattice[0][0] + "       " + '%0.15f'%new_lattice[0][1] + "       " + '%0.15f'%new_lattice[0][2] + "\n")
            outputfile.write("       " + '%0.15f'%new_lattice[1][0] + "       " + '%0.15f'%new_lattice[1][1] + "       " + '%0.15f'%new_lattice[1][2] + "\n")
            outputfile.write("       " + '%0.15f'%new_lattice[2][0] + "       " + '%0.15f'%new_lattice[2][1] + "       " + '%0.15f'%new_lattice[2][2] + "\n")
            outputfile.write("%ENDBLOCK LATTICE_CART\n")
            outputfile.write("FIX_ALL_CELL TRUE\n")
            in_lattice = True
        elif (not in_lattice):
            outputfile.write(line)
    inputfile_0.close()
    outputfile.close()
    shutil.copy(seedname+".param",seedname+"_"+change+".param")
    os.system('scastep -np 24 -walltime 0:00:00 -mem 64000 -seedname '+seedname+'_'+change+' wm'+seedname+'_'+change)
    inputfile_1 = open(seedname+'_'+change+'.castep',"r")
    
    volume = []
    density = []
    inputfile_1 = open("W_OTFGultra.castep", "r")
    outputfile_1 = open(seedname+".summary","a")
    for line in inputfile_1:
        #print(line)
        if re.search(volume_RE,line):
            volume.append(re.search(volume_RE,line).group(1))
        elif re.search(density_RE,line):
            density.append(re.search(density_RE,line).group(1))
        elif re.search(energy_RE,line):
            energy = re.search(energy_RE,line).group(1)
    
    outputfile_1.write(volume[-1]+' '+density[-1]+' '+energy+'\n')
    inputfile_1.close()
    outputfile_1.close()
