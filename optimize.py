# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 10:43:19 2020

@author: Maywell2019
"""

import re
import shutil
import os
import numpy as np
from itertools import combinations

################################## regular rules ##########################
dotcell_lattice_RE = re.compile("""^\s*%BLOCK\s+LATTICE_(?:CART|ABC)*\n
                                \s+([\+\-]?\d+.\d+)\s+\s+([\+\-]?\d+.\d+)\s+\s+([\+\-]?\d+.\d+)\s*\n
                                \s+([\+\-]?\d+.\d+)\s+\s+([\+\-]?\d+.\d+)\s+\s+([\+\-]?\d+.\d+)\s*\n
                                \s+([\+\-]?\d+.\d+)\s+\s+([\+\-]?\d+.\d+)\s+\s+([\+\-]?\d+.\d+)\s*\n""",re.VERBOSE)
dotcell_lattice_start_RE = re.compile("^\s*%BLOCK\s+LATTICE_(?:CART|ABC)",re.VERBOSE)
dotcell_lattice_end_RE = re.compile("^\s*%ENDBLOCK\s+LATTICE_(?:CART|ABC)",re.VERBOSE)
dotcell_atoms_start_RE = re.compile("^\s*%BLOCK\s+POSITIONS_(?:FRAC|ABS)", re.VERBOSE)
dotcell_atoms_end_RE = re.compile("^\s*%ENDBLOCK\s+POSITIONS_(?:FRAC|ABS)", re.VERBOSE)
MASS_start_RE=re.compile("^\s*%BLOCK\s+SPECIES_(?:MASS|ABS)", re.VERBOSE)
MASS_end_RE = re.compile("^\s*%ENDBLOCK\s+SPECIES_(?:MASS|ABS)", re.VERBOSE)
POT_start_RE=re.compile("^\s*%BLOCK\s+SPECIES_(?:POT|ABS)", re.VERBOSE)
POT_end_RE = re.compile("^\s*%ENDBLOCK\s+SPECIES_(?:POT|ABS)", re.VERBOSE)
LCAO_start_RE=re.compile("^\s*%BLOCK\s+SPECIES_LCAO_STATES", re.VERBOSE)
LCAO_end_RE = re.compile("^\s*%ENDBLOCK\s+SPECIES_LCAO_STATES", re.VERBOSE)
########################## atom database #######################################
atom_Al={'name':'Al','mass':'26.9820003510','POT':'3|2|3.5|5|7|30UU:31UU:32UU[]','LCAO':'2'}
atom_Co={'name':'Co','mass':'58.9329986572','POT':'3|2.2|2|1|9|12|14|40UU:32UU:41UU(qc=5.5)[]','LCAO':'3'}
atom_Cr={'name':'Cr','mass':'51.9959983826','POT':'3|2|2|0.7|11|12|14|30U:40UU:31UU:32UU(qc=6)[]','LCAO':'5'}
atom_Hf={'name':'Hf','mass':'178.4900054932','POT':'3|2.1|12|14|16|50U:52UU:60UU:51UU:43UU(qc=6)[]','LCAO':'3'}
atom_Mg={'name':'Mg','mass':'24.3050003052','POT':'3|1.8|10|14|18|20U:30UU:21UU:32UU[]','LCAO':'4'}
atom_Mo={'name':'Mo','mass':'95.9400024414','POT':'3|1.6|11|14|15.5|40U:50UU:41UU:42UU(qc=6)[]','LCAO':'5'}
atom_Nb={'name':'Nb','mass':'92.9059982300','POT':'3|1.6|10|12.9|14.7|40U:50UU:41UU:42UU(qc=6)[]','LCAO':'5'}
atom_Ta={'name':'Ta','mass':'180.9479980469','POT':'3|2.4|9|11.4|13.2|50U:52UU:60UU:51UU:43UU(qc=6)[]','LCAO':'3'}
atom_Ti={'name':'Ti','mass':'47.9000015259','POT':'3|1.8|10|12|14|30U:40UU:31UU:32UU(qc=5.5)[]','LCAO':'5'}
atom_W={'name':'W','mass':'183.8500061035','POT':'3|2.4|9.5|11|13|50U:60UU:51UU:52UU:43UU(qc=6)[]','LCAO':'3'}
atom_Zr={'name':'Zr','mass':'91.2200012207','POT':'3|2.1|7|8.5|10|40U:50UU:41UU:42UU[]','LCAO':'5'}

atom_list=[atom_Cr,atom_Hf,atom_Mo,atom_Nb,atom_Ta,atom_Ti,atom_W,atom_Zr]
############################### list all alloys ###############################
quaternary_combination=list(combinations(atom_list,4))
quinary_combination=list(combinations(atom_list,5))
hex_combination=list(combinations(atom_list,6))
##############################atom position in body-centered cubic####################
position_1=['0.0000000000000000','   ','0.0000000000000000','   ','0.0000000000000000']
position_2=['0.5000000000000000','   ','0.5000000000000000','   ','0.5000000000000000']

################## obtain lattice parameters from old.cell file #################
inputfile_0 = open("W_OTFGultra.cell", "r")
result = dotcell_lattice_RE.findall(inputfile_0.read())[0]
lattice = []
lattice.append([float(result[0]),float(result[1]),float(result[2])])
lattice.append([float(result[3]),float(result[4]),float(result[5])])
lattice.append([float(result[6]),float(result[7]),float(result[8])])
lattice = np.array(lattice)
new_lattice = np.dot(lattice,1.0)
inputfile_0.close()

###################### make new .cell file ##########################
in_lattice = False
in_atoms = False
in_mass = False
in_pot = False
in_LCAO = False
for alloy in quinary_combination:
    alloy_name = alloy[0]['name']+alloy[1]['name']+alloy[2]['name']+alloy[3]['name']+alloy[4]['name']
    inputfile = open("W_OTFGultra.cell", "r")
    outputfile = open(alloy_name+".cell","w")
    present_dir = os.getcwd()
    os.mkdir(alloy_name)
    for line in inputfile:
######### write lattice parameters ##############################        
        if (re.search(dotcell_lattice_end_RE,line) and in_lattice):
            in_lattice = False
        elif (re.search(dotcell_lattice_start_RE,line) and not in_lattice):
            outputfile.write("%BLOCK LATTICE_CART\n")
            outputfile.write("       " + '%0.15f'%new_lattice[0][0] + "       " + '%0.15f'%new_lattice[0][1] + "       " + '%0.15f'%new_lattice[0][2] + "\n")
            outputfile.write("       " + '%0.15f'%new_lattice[1][0] + "       " + '%0.15f'%new_lattice[1][1] + "       " + '%0.15f'%new_lattice[1][2] + "\n")
            outputfile.write("       " + '%0.15f'%new_lattice[2][0] + "       " + '%0.15f'%new_lattice[2][1] + "       " + '%0.15f'%new_lattice[2][2] + "\n")
            outputfile.write("%ENDBLOCK LATTICE_CART\n")
            in_lattice = True
########## write atoms position ################################        
        elif (re.search(dotcell_atoms_end_RE,line) and in_atoms):
            in_atoms = False
        elif (re.search(dotcell_atoms_start_RE,line) and not in_atoms):
            outputfile.write("%BLOCK POSITIONS_FRAC\n")
            for element in alloy:
                outputfile.write("  "+element['name']+"   "+position_1[0]+position_1[1]+position_1[2]+position_1[3]+position_1[4]+" "+"MIXTURE:(   1  0.200000)"+ "\n")
                outputfile.write("  "+element['name']+"   "+position_2[0]+position_2[1]+position_2[2]+position_2[3]+position_2[4]+" "+"MIXTURE:(   2  0.200000)"+ "\n")
            outputfile.write("%ENDBLOCK POSITIONS_FRAC\n")
            in_atoms = True
########## wrire atoms mass ######################################        
        elif (re.search(MASS_end_RE,line) and in_mass):
            in_mass = False
        elif (re.search(MASS_start_RE,line) and not in_mass):
            outputfile.write("%BLOCK SPECIES_MASS\n")
            for element in alloy:
                outputfile.write("       "+element['name']+"    "+element['mass']+ "\n")
            outputfile.write("%ENDBLOCK SPECIES_MASS\n")
            in_mass = True
########## write atoms potential #################################
        elif (re.search(POT_end_RE,line) and in_pot):
            in_pot = False
        elif (re.search(POT_start_RE,line) and not in_pot):
            outputfile.write("%BLOCK SPECIES_POT\n")
            for element in alloy:
                outputfile.write("       "+element['name']+"  "+element['POT']+ "\n")
            outputfile.write("%ENDBLOCK SPECIES_POT\n")        
            in_pot = True
########### write atoms LCAO_STATES ###############################
        elif (re.search(LCAO_end_RE,line) and in_LCAO):
            in_LCAO = False
        elif (re.search(LCAO_start_RE,line) and not in_LCAO):
            outputfile.write("%BLOCK SPECIES_LCAO_STATES\n")
            for element in alloy:
                outputfile.write("       "+element['name']+"         "+element['LCAO']+ "\n")
            outputfile.write("%ENDBLOCK SPECIES_LCAO_STATES\n")
            in_LCAO = True
            
        elif (not (in_lattice or in_atoms or in_mass or in_pot or in_LCAO)):
            outputfile.write(line)
    inputfile.close()
    outputfile.close()
    shutil.move(alloy_name+".cell",present_dir+"/"+alloy_name+"/"+alloy_name+".cell")
    shutil.copy("W_OTFGultra.param",present_dir+"/"+alloy_name+"/"+alloy_name+".param")
#    os.chdir(present_dir+"/"+alloy_name+"/")
#    os.system('scastep -np 24 -walltime 0:00:00 -mem 64000 -seedname '+alloy_name+' wm'+alloy_name)
#    os.system('cp ../energy-lattice.py .')
#    os.system('python energy-lattice.py')
#    os.system('cd ../')