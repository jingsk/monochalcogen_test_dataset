#Name: Tina Mihm
#Date: Nov 11, 2024
#Description: Pulls the atoms that are in between the layers and shifts them by a set amount in the x and y direction to create the double well structures 

import os
import numpy as np
import pandas as pd
import math as ma
from ase import Atoms
from ase.geometry.analysis import Analysis
from ase.io import read, write
import matplotlib.pyplot as plt
from ase.geometry import wrap_positions
from ase.io.vasp import write_vasp

#---------------------------------------------------------
## Reads in the POSCAR and pulls out required info
#---------------------------------------------------------
## Read in POSCAR
atoms = read("CONTCAR-Original.vasp")
atoms2 = read("CONTCAR-Original.vasp")

## pull out atoms and cell 
pos = atoms.get_positions()
Cell = atoms.get_cell()

## get lattice parameters to scale coordinates 
print("lattice:")
print(Cell)
A = Cell[0][0]
B = Cell[1][1]
C = Cell[2][2]
print(A, B, C)

## separates the upper and lower layer via indicies 
print("atomic positions:")
print(pos)
print("atomic positions for Ge:")
pos_Ge = pos[:4]
pos_Se = pos[4:]
print(pos[:4])
pos_list = pos.tolist()
# print(pos_list)
Ind_l = []
Ind_u = []

for a in pos: 
    if a[2] < 5.0: 
        # print("lower layer:")
        # print(a)
        # ind = np.where(pos == a)
        # print(ind)
        ind = pos_list.index(a.tolist())
        # print(ind[0][2])
        # print(ind)
        # Ind_l +=[ind[0][2]]
        Ind_l +=[ind]
    if a[2] > 5.0: 
        # print("upper layer:")
        # print(a)
        # ind = np.where(pos == a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        # Ind_u +=[ind[0][2]]
        Ind_u +=[ind]
print("Lower layer index:")
print(Ind_l)
print("Upper layer index:")
print(Ind_u)

Ind_l_Ge = []
Ind_u_Ge = []

# for a in pos_Ge: 
#     if a[2] < 5.0: 
#         # print("lower layer:")
#         print(a)
#         # ind = np.where(pos == a)
#         # print(ind)
#         ind = pos_list.index(a.tolist())
#         # print(ind[0][2])
#         print(ind)
#         # Ind_l +=[ind[0][2]]
#         Ind_l_Ge +=[ind]
#     if a[2] > 5.0: 
#         # print("upper layer:")
#         # print(a)
#         # ind = np.where(pos == a)
#         ind = pos_list.index(a.tolist())
#         # print(ind[0][0])
#         # Ind_u +=[ind[0][2]]
#         Ind_u_Ge +=[ind]
# print("Lower layer Ge index:")
# print(Ind_l_Ge)
# print("Upper layer Ge index:")
# print(Ind_u_Ge)

for a in pos_Ge: 
    if a[2] < 3.0: 
        # print("lower layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_l_b_Ge =ind
    if a[2] < 5.0 and a[2] > 3.0 : 
        # print("lower layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_l_t_Ge =ind
    if a[2] > 5.0 and a[2] < 8.0: 
        # print("upper layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_u_b_Ge =ind
    if a[2] > 8.0: 
        # print("upper layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_u_t_Ge =ind
print("Lower layer Ge atoms index:")
print(Ind_l_b_Ge, Ind_l_t_Ge)
print("Upper layer Ge atoms index:")
print(Ind_u_b_Ge, Ind_u_t_Ge)

for a in pos_Se: 
    if a[2] < 3.0: 
        # print("lower layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_l_b_Se =ind
    if a[2] < 5.0 and a[2] > 3.0 : 
        # print("lower layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_l_t_Se =ind
    if a[2] > 5.0 and a[2] < 8.0: 
        # print("upper layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_u_b_Se =ind
    if a[2] > 8.0: 
        # print("upper layer:")
        # print(a)
        ind = pos_list.index(a.tolist())
        # print(ind[0][0])
        Ind_u_t_Se =ind
print("Lower layer Se atoms index:")
print(Ind_l_b_Se, Ind_l_t_Se)
print("Upper layer Se atoms index:")
print(Ind_u_b_Se, Ind_u_t_Se)

## set up the difference for a and b coordinates based on an angle:

#Max difference: 
Diff = pos[4] - pos[3] 
# print(Diff)
bond = np.sqrt(Diff[0]**2 + Diff[1]**2 + Diff[2]**2 ) ##hypotenus of triangle = bond length 
print("this is the bond length (i.e hypotenus):", bond)
B_Delta_b = abs(pos[4][1] - pos[3][1]) ## oposite side of tringle from theta
print("this is the difference in b-coordinates (i.e opposite side from theta):", B_Delta_b)
B_theta = ma.asin(B_Delta_b/bond)
print("this is the theta for the B point in the double well:", round(B_theta, 2))

theta_pos= [0.3, 0.26, 0.22, 0.21, 0.2, 0.18, 0.16, 0.14, 0.12, 0.1, 0.08, 0.06, 0.04, 0.02, 0.0]
theta_neg = [t *-1 for t in theta_pos[:-1]]
# print(theta_neg)

print("This is the change in b for the positive theta:")
Delta_b_pos = []
for t in theta_pos: 
    # print(t)
    db = bond * ma.sin(t)
    print(db)
    Delta_b_pos +=[db]
# print(Delta_b_pos)

print("This is the change in b for the negative theta:")
Delta_b_neg = []
db_0 = Delta_b_pos[-1]
for t in reversed(theta_neg): 
    # print(t)
    db = bond * ma.sin(t)
    # db = db + B_Delta_b
    print(db)
    Delta_b_neg +=[db]

Delta_b = Delta_b_pos + Delta_b_neg
print(Delta_b)
theta = theta_pos 
theta.extend(reversed(theta_neg))
print(theta)

#-------------------------------------------------------------------
## Calculates the coordinates for the new atoms and generates POSCARs
#-------------------------------------------------------------------

bpos_l_l = pos[Ind_l_t_Se][1]
bpos_l_r = pos[Ind_l_b_Se][1]
bpos_u_l = pos[Ind_u_t_Se][1]
bpos_u_r = pos[Ind_u_b_Se][1]
print(bpos_l_l)
print(bpos_l_r)
print(bpos_u_l)
print(bpos_u_r)
c = 1
for d in Delta_b:   
    # print("lower level:")
    # print(pos[i][1])
    # print("old pos")
    # print(pos[Ind_l_l_b])
    # print(b)
    newb = bpos_l_l + d
    # print(newb)
    pos_u_new = [pos[Ind_l_b_Ge][0], newb, pos[Ind_l_b_Ge][2]]
    # print("new pos")
    # print(pos_u_new)
    atoms[Ind_l_b_Ge].position = pos_u_new
    #########
    # print("old pos")
    # print(pos[Ind_l_r_t])
    # print(b)
    newb2 = bpos_l_r + d
    # print(newb)
    pos_u_new = [pos[Ind_l_t_Ge][0], newb2, pos[Ind_l_t_Ge][2]]
    # print("new pos")
    # print(pos_u_new)
    atoms[Ind_l_t_Ge].position = pos_u_new
    ##################
    # print("upper level:")
    # print(pos[i][1])
    # print("old pos")
    # print(pos[Ind_u_l_b])
    # print(b)
    newb = bpos_u_l - d
    # print(newb)
    pos_u_new = [pos[Ind_u_b_Ge][0], newb, pos[Ind_u_b_Ge][2]]
    # print("new pos")
    # print(pos_u_new)
    atoms[Ind_u_b_Ge].position = pos_u_new
    #########
    print("old pos")
    print(pos[Ind_u_t_Ge])
    # print(b)
    newb = bpos_u_r - d
    # print(newb)
    pos_u_new = [pos[Ind_u_t_Ge][0], newb, pos[Ind_u_t_Ge][2]]
    print("new pos")
    print(pos_u_new)
    atoms[Ind_u_t_Ge].position = pos_u_new
    # #---------------------------------------------------------
    # ##  wraps new postions to stay within the unitcell
    # ##  printing out POSCARs
    # #---------------------------------------------------------
    pos2 = atoms.get_positions()
    print(pos2)
    atoms_new = wrap_positions(pos2, Cell)
    atoms2.set_positions(atoms_new, apply_constraint=False)
    write_vasp("POSCARs/POSCAR_Theta_"+str(theta[Delta_b.index(d)]), atoms2)
    ##
    #------------------------------------------------------
    ## Check-point for moving along the b axix (zig-zag)
    #------------------------------------------------------
    pos3 = atoms2.get_positions()
    print(pos3)
    x_u = [pos3[Ind_u[1]][1], pos3[Ind_u[2]][1], pos3[Ind_u[0]][1], pos3[Ind_u[3]][1]]
    y_u = [pos3[Ind_u[1]][2], pos3[Ind_u[2]][2], pos3[Ind_u[0]][2], pos3[Ind_u[3]][2]]
    df = {"x": x_u, "y": y_u}
    df = pd.DataFrame(df)
    df = df.sort_values("x")
    x_l = [pos3[Ind_l[3]][1], pos3[Ind_l[0]][1], pos3[Ind_l[2]][1], pos3[Ind_l[1]][1]]
    y_l = [pos3[Ind_l[3]][2], pos3[Ind_l[0]][2], pos3[Ind_l[2]][2], pos3[Ind_l[1]][2]]
    df2 = {"x": x_l, "y": y_l}
    df2 = pd.DataFrame(df2)
    df2 = df2.sort_values("x")
    # ##
    # # x_u = [pos2[Ind_u[1]][1], pos2[Ind_u[2]][1], pos2[Ind_u[0]][1], pos2[Ind_u[3]][1]]
    # # y_u = [pos2[Ind_u[1]][2], pos2[Ind_u[2]][2], pos2[Ind_u[0]][2], pos2[Ind_u[3]][2]]
    # # x_l = [pos2[Ind_l[3]][1], pos2[Ind_l[0]][1], pos2[Ind_l[2]][1], pos2[Ind_l[1]][1]]
    # # y_l = [pos2[Ind_l[3]][2], pos2[Ind_l[0]][2], pos2[Ind_l[2]][2], pos2[Ind_l[1]][2]]
    # ## Check-point for moving along the b axix (zig-zag)
    plt.figure(1)
    fig, ax1 = plt.subplots()
        ##
    pl1 = ax1.plot("x", "y", data = df, linestyle = "--", marker = "o",label = "B', Original", color = "#0000FF") ## blue
    pl1 = ax1.plot("x", "y", data = df2, linestyle = "--", marker = "o", color = "#FF007F") ## pink  
    ax1.set_xlim(0, B)
    ax1.set_ylim(0, C)
    ax1.set_ylabel(r"C coordinate")
    ax1.set_xlabel(r"B coordinate")
    plt.savefig("Figs/"+str(c)+"-Theta_"+str(theta[Delta_b.index(d)])+".png", bbox_inches='tight')
        ##
    plt.show()
    c = c+1

print("There should be ", len(theta), " POSCARs in the POSCARs/ directory")