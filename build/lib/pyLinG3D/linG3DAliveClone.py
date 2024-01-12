
"""
This is a companion code for the paper  
"LinG3D: Visualizing the Spatio-Temporal Dynamics of Clonal Evolution" 
A. Hu, A.M.E. Ojwang', K.D. Olumoyin, and K.A. Rejniak

This code generates the 3D lineage tree of one clone of number
specified in 'cloneNum' taking into account only the cells that
survived to the end of simulation. It uses data from directory
 'pathdata'. 

The following parameters needs to be specified:  
  pathdata -- directory with input data  
  CloneNum  -- clone number to be drawn  
  
It requires the following data in the pathdata/data/ directory:  
  cell_history.txt -- file with info about each cell  
  cellID_##.txt    -- cell IDs in a file with index number ##  
  cellXY_##.txt    -- cell coordinates in a file with index ##  
  drug.txt         -- concentration of a drug for background  
  
The following parameters are project-dependent and should be  
specified for a given project:  
  xmin,xmax,ymin,ymax -- dimensions of the spacial domain  
  tmin, tmax          -- dimensions of the temporal domain  
  fileStep            -- frequency of the sampled data

for the examples discussed in the paper use:  
 example 1: pathdata='exampleB05';  cloneNum between 0 and 9
 example 2: pathdata='exampleB005'; cloneNum between 0 and 147 

October 31, 2022
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from os import path

from .utils import patch_b
from .utils import patch_c
from .utils import patch_y
from .utils import patch_r
from .utils import DrawBackground
from .utils import DefineColorPalette
from .utils import find_b

def linG3DAliveClone(pathdata,cloneNum,IsGradient,xmin,xmax,ymin,ymax,tmin,tmax,fileStep,toPrint):
    dataDirectory = '/data/' # directory with cell and drug data
    timeStep=(tmax-tmin)/(2.5*(xmax-xmin))

    # plt.grid()
    # ax.set_xlim([xmin, xmax])
    # ax.set_ylim([tmin, tmax // timeStep])
    # ax.set_zlim([ymin, ymax])
    
    col=DefineColorPalette()
    Ncol=col[:,0].shape[0]

    if toPrint==1:
        pathFigs=pathdata+'/fig_clonesAlive';
        if path.exists(pathFigs):
            pass
        else:
            os.mkdir(pathFigs)
    
    # draw background with drug gradient
    if IsGradient == 1:
        drug = np.loadtxt(pathdata + dataDirectory + 'drug.txt')
        DrawBackground(drug, tmax, timeStep, xmin, xmax, ymin, ymax, ax)
        
    # load cell history file
    hist = np.loadtxt(pathdata + dataDirectory + 'cell_history.txt')
    # [cell ID, clone ID, mother ID, birth iter, div / death iter]
    
    # load indices of all survived cells
    cellID = np.loadtxt(pathdata + dataDirectory + 'cellID_' + str(tmax) + '.txt')
    # print(cellID.shape)  
    
    print('clone=' + str(cloneNum))
    
    # identify all survived cells from a given clone
    indLast = []
    Nlast = 0
    
    for ii in range(len(cellID)):  # all cells in the last file
        if hist[int(cellID[ii])-1,1] == cloneNum:
            indLast.append(cellID[ii])
            Nlast += 1
    Numlast = Nlast
    for ii in range(Numlast):
            
        moNum = hist[int(indLast[ii])-1,2]
        while moNum > 0 and (hist[int(moNum)-1,1] == cloneNum):  # all predecesor cells
            indLast.append(moNum)
            Nlast += 1
            moNum = hist[int(moNum)-1,2]
    indLast = np.unique(indLast)  # remove repeated cells
    # print(indLast)
        
    # define matrix of line segments (3D branches) to draw
    matrix_to_draw = np.zeros((1, 6))  # [x1, t1, y1, x2, t2, y2] X - time - Y axes
    Nmatrix = 0
        
    for ii in range(len(indLast)):  # for every cell with index in indLast
        if ii % 100 == 0:
            print('... calculating')
                
        cellNum = hist[int(indLast[ii])-1,0]  # cell ID
        # print(cellNum)
        mothNum = hist[int(indLast[ii])-1,2]  # mother ID
        strtNum = hist[int(indLast[ii])-1,3]  # cell birth
        endNum = hist[int(indLast[ii])-1,4]
        endNum = max(tmax, min(endNum, tmax)) # cell div / death / tmax
            
        # find all appearances of the cellNum
        kkStart = fileStep * np.floor(strtNum / fileStep)  # initial file number
        kkEnd = fileStep * np.floor(endNum / fileStep)  # final file number
        kkStart = int(kkStart)
        kkEnd = int(kkEnd)
        # print(kkStart)
            
        for kk in range(kkEnd, kkStart, - fileStep):  # inspect all files
            # cell ID and cell XY from the first file
            fileMeID = np.loadtxt(pathdata + dataDirectory + 'cellID_' + str(kk) + '.txt')  
            fileMeXY = np.loadtxt(pathdata + dataDirectory + 'cellXY_' + str(kk) + '.txt')
                
            indMe = find_b(fileMeID == cellNum) # find current indices of cellID
            # print(indMe)
                
            fileMe2ID = np.loadtxt(pathdata + dataDirectory + 'cellID_' + str(kk - fileStep) + '.txt')
            fileMe2XY = np.loadtxt(pathdata + dataDirectory + 'cellXY_' + str(kk - fileStep) + '.txt')
                
            indMe2 = find_b(fileMe2ID == cellNum)  # find current indices of cellID
            # print(indMe2)
                
            if indMe.size == 0:
                pass
            elif indMe2.size == 0:
                # print(indMe)                  
                while kkStart < hist[int(mothNum)-1,3]:  # find file with the grand-mother cell
                    mothNum = hist[int(mothNum)-1,2]
                fileMe2ID = np.loadtxt(pathdata + dataDirectory + 'cellID_' + str(kkStart) + '.txt')
                fileMe2XY = np.loadtxt(pathdata + dataDirectory + 'cellXY_' + str(kkStart) + '.txt')
                    
                indMe2 = find_b(fileMe2ID == mothNum)  # find current indices of mother cellID
                # print(indMe2)
                if indMe2.size == 0:
                    pass
                else:
                    matrix_to_draw = np.row_stack((matrix_to_draw, np.zeros((1, 6))))
                    matrix_to_draw[Nmatrix, :] = [fileMeXY[indMe][0][0], kkStart + fileStep,
                                                fileMeXY[indMe][0][1], fileMe2XY[indMe2][0][0],
                                                kkStart, fileMe2XY[indMe2][0][1]]
                    Nmatrix += 1  # save branch to draw [x1,t1,y1,x2,t2,y2]
            else:
                matrix_to_draw = np.row_stack((matrix_to_draw, np.zeros((1, 6))))
                matrix_to_draw[Nmatrix, :] = [fileMe2XY[indMe2][0][0], kk - fileStep,
                                                  fileMe2XY[indMe2][0][1], fileMeXY[indMe][0][0],
                                                  kk, fileMeXY[indMe][0][1]]
                Nmatrix += 1  # save branch to draw [x1,t1,y1,x2,t2,y2]


    
    # drawing clones
    NumCol = cloneNum % Ncol  # clone color
    # plt.title('survived cells from clone=' + str(cloneNum))
    plt.ylabel('iterations/time x ' + str(timeStep),fontsize=8)

    for ii in range(Nmatrix):
        x = [matrix_to_draw[ii,0], matrix_to_draw[ii,3]]
        y = [matrix_to_draw[ii,1], matrix_to_draw[ii,4]]
        y = [yi / timeStep for yi in y]
        z = [matrix_to_draw[ii,2], matrix_to_draw[ii,5]]
        ax.plot(x, y, z, c=col[NumCol,:3],linewidth=1.3)
    ax.view_init(37.5,-130)
    ax.set_box_aspect([1,2.5,1])
    ax.set_xticks([-100,-50,0,50,100])
    ax.set_yticks([500,450,400,350,300,250,200,150,100,50,0])
    ax.set_zticks([-100,-50,0,50,100])
    ax.set_ylim(0,500)
        
    if toPrint==1:
        plt.savefig(pathFigs + "/tree_alive_clone_"+ str(cloneNum)+'.jpg', dpi=300)
    
    plt.show()
    
    
fig = plt.figure()
ax = plt.axes(projection='3d')

           