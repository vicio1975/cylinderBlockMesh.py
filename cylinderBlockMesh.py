# -*- coding: utf-8 -*-
"""
Author: Vincenzo Sammartano
"""
#### Libraries
import numpy as np
import shutil
#import matplotlib.pyplot as pl
####


#### Function definitions
def headerLines(cl,loc,obj,con="convertToMeters 1;"):
    h = [
    "/*--------------------------------*- C++ -*----------------------------------*\\",
    "| =========                 |                                                 |",
    "| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |",
    "|  \\\    /   O peration     | Version:  v4.x                                  |",
    "|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |",
    "|    \\\/     M anipulation  |                                                 |",
    "\*---------------------------------------------------------------------------*/",
    "FoamFile",
    "{",
    "    version     2.0;",
    "    format      ascii;",
    "    class       {};",
    '    location    "{}";',
    "    object      {};",
    "}",
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",
    "",
    "{}",
    ""]
   
    h[11] = h[11].format(cl)
    h[12] = h[12].format(loc)
    h[13] = h[13].format(obj)
    h[17] = h[17].format(con)    
    return h
        
    
def Geometry():
    ############################
    ## vertices and their coord
    px = []
    py = []
    pz = []

    px.append(C[0]-R*np.cos(alpha))  #0
    px.append(C[0]+R*np.cos(alpha))  #1
    px.append(C[0]+a/2)              #2
    px.append(C[0]-a/2)              #3

    px.append(C[0]-a/2)              #4
    px.append(C[0]+a/2)              #5
    px.append(C[0]+R*np.cos(alpha))  #6
    px.append(C[0]-R*np.cos(alpha))  #7

    
    py.append(C[1]-R*np.sin(alpha))     #0
    py.append(C[1]-R*np.sin(alpha))     #1
    py.append(C[1]-b/2)                 #2
    py.append(C[1]-b/2)                 #3

    py.append(C[1]+b/2)                 #4
    py.append(C[1]+b/2)                 #5
    py.append(C[1]+R*np.sin(alpha))     #6
    py.append(C[1]+R*np.sin(alpha))     #7

    pz.append(C[2])
    for i in range(7):
        pz.append(pz[0])
    pz.append(pz[-1] + H)
    for i in range(7):
        pz.append(pz[-1])
        
    px = px + px
    py = py + py
   
    Ixyz = [index for index,i in enumerate(px)]
    
    ##########################################
    ### Write file block #####################
    fid = open("blockMeshDict","w")
    
    #header
    header = headerLines(cl,loc,obj)
    for l in range(len(header)):
        fid.write("{}\n".format(header[l]))
    ####
    
    #vertices
    fid.write("vertices\n")    
    fid.write("(\n")
    for i in range(Nnod):
        fid.write("\t({} {} {}) \t\t //{}\n".format(px[i],py[i],pz[i],i))
    fid.write(");\n")
    fid.write("")
    ####
    
    #blocks
    fid.write("\nblocks\n")
    fid.write("(\n")
    block0 = [3,2,5,4] # Nx & Ny 
    block0 = block0 + [i+int(Nnod/2) for i in block0]

    block1 = list(range(4)) # Nx = Nx & Ny = Nd
    block1 = block1 + [i+int(Nnod/2) for i in block1]

    block2 = [2,1,6,5] # Nx = Nd & Ny = Ny
    block2 = block2 + [i+int(Nnod/2) for i in block2]

    block3 = [4,5,6,7] # Nx = Nx & Ny = Nd 
    block3 = block3 + [i+int(Nnod/2) for i in block3]

    block4 = [0,3,4,7] # Nx = Nd & Ny = Ny
    block4 = block4 + [i+int(Nnod/2) for i in block4]

    #block list and number of elements for each block
    block = [block0,block1,block2,block3,block4]  
    NN = [[Nx0,Ny0,Nz0],[Nx0,Nd0,Nz0],[Nd0,Ny0,Nz0],[Nx0,Nd0,Nz0],[Nd0,Ny0,Nz0]]
    nn = []

    for ind,ii in enumerate(block):
        nn = NN[ind]
        fid.write("\t hex (") 
        for i in ii:
            fid.write(" {} ".format(Ixyz[i]))
        fid.write(")\t")
        ##number of subdivisions
        fid.write("({} {} {})\t".format(nn[0],nn[1],nn[2]))
        ##simpleGrading
        fid.write("simpleGrading ({} {} {})\n".format(sx1,sy1,sz1))
    
    fid.write(");\n")
    fid.write("")
    ###
    
    #### edges
    fid.write("\nedges\n")
    fid.write("(\n")
    ee = [[0,1],[1,6],[6,7],[7,0]]
    eeTemp = []
    for eee in ee:
        eeTemp.append(list(map(lambda x: x+8,eee)))
    ee = ee + eeTemp

    for j in range(8):
        E = ee[j]
        fid.write("\t arc {} {}\t({} {} {})\n".format(E[0], E[1], inpX[j],inpY[j],inpZ[j]))
    fid.write(");\n")   
    fid.write("")
   
    #### boundary
    #bases
    ss0 = [0,1,2,3]
    ss1 = [3,2,5,4]
    ss2 = [1,6,5,2]  
    ss3 = [6,7,4,5] 
    ss4 = [0,3,4,7] 
    SS = [ss0,ss1,ss2,ss3,ss4]
    morb=[]
    for S in SS:
        morb.append([i+8 for i in S])
    SS = SS + morb
       
    fid.write("\nboundary\n")
    fid.write("(\n")
    for j in range(2):
        fid.write("\t{}\n".format(bcs[j]))
        fid.write("\t{\n")
        fid.write("\t\ttype {};\n".format(tipos[bcs[j]]))
        fid.write("\t\tfaces\n")
        fid.write("\t\t(\n")  
        for i in list(range(j*5,5*(j+1))):
            sS = SS[i]
            st = "\t\t\t({a} {b} {c} {d})\n".format(a=str(sS[0]),b=str(sS[1]),c=str(sS[2]),d=str(sS[3]))
            fid.write(st)
        fid.write("\t\t);\n")
        fid.write("\t}\n")
    #walls
    ss = [0,1,1,6,6,7,7,0]
    fid.write("\t{}\n".format(bcs[2]))
    fid.write("\t{\n")
    fid.write("\t\ttype {};\n".format(tipos[bcs[2]]))
    fid.write("\t\tfaces\n")
    fid.write("\t\t(\n")
    for i in list(range(0,8,2)):
        fid.write("\t\t\t({a} {b} {c} {d})\n".format(a=str(ss[i]),b=str(ss[i+1]),c=str(ss[i]+8),d=str(ss[i+1]+8)))
    fid.write("\t\t);\n")
    fid.write("\t}\n")
           
    fid.write(");\n")
    fid.write("")
    
    #mergPatchPairs
    fid.write("\nmergPatchPairs\n")
    fid.write("(\n")
    fid.write(");\n")
    fid.write("\n")
    ####
    fid.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")
    print("\n--> Moving the blockMeshDict file in system ...")
    shutil.move("blockMeshDict","system/blockMeshDict" )
############################# END Functions ##############################################################


################ MAIN Program
cl = "dictionary"
loc = "system"
obj = "blockMeshDict"
tipos = {}
bcs = "inlet outlet walls".split(" ")

print("\n--> Selection of BCs")
for bc in bcs:
    q = "    - {}\t: ".format(bc)
    tipos[bc] = input(q).strip()

########################
#Geometry specifications
print("\n--> Geometry definition")
##a = float((input("   - Core side       a: ").strip()))
#Radius of cylinder
R = float((input("   - Cylinder Radius R: ").strip()))
H = float((input("   - Cylinder height H: ").strip()))
#center position
print("\n--> Axis position:")
c1 = float((input("   - x0 = ").strip()))
c2 = float((input("   - y0 = ").strip()))
c3 = float((input("   - z0 = ").strip()))
C = [c1,c2,c3]
#core side
encore = 10*H/R 
a = R * 1/(1+1/encore)
b = a
#diagonal geometry
alpha = np.arctan(b/a)
diag = 0.5 * a * np.sqrt(2)
#Number of nodes
Nnod = 16

#cells for block 0 
lx = float((input("\n--> Length of cells dx = dy = ").strip()))
Nx0 = int(np.ceil(a/lx))
Ny0 = Nx0 
Nz0 = int(np.floor(H/lx))
llat = R - diag
Nd0 = int(np.floor(llat/lx))

if Nd0 == 0:
    Nd0 = 2

#simpleGrading
sx1=1
sy1=1
sz1=1

#interpolation points
bis =  [ 1.5*np.pi, 2*np.pi, 0.5*np.pi, 1.0*np.pi ]
inpX = [ R*np.cos(t)+C[0] for t in bis ] 
inpY = [ R*np.sin(t)+C[1] for t in bis ]
inpZ =  []
for i in range(4):
    inpZ.append(C[2])
inpX = inpX + inpX
inpY = inpY + inpY
inpZ = inpZ + [i+H for i in inpZ]

####################### 

Geometry()

####################### End Code
