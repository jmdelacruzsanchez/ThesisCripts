# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 15:25:30 2021

@author: jdelacruz
"""

import pickle
import os
from read_roi import read_roi_zip
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import varname


def CalculateStimulationEficacy(InputDic):
    PulsesNumber = len(InputDic["Pulses"])
    OutputDic = {}
    
    for roi, FluoPulses in InputDic.items():
        if roi == "Pulses":
            continue
        OutputDic[roi.split("(")[1].split(")")[0]] = len(FluoPulses[0:-1])/PulsesNumber
        
    return OutputDic


def ROI_DicToArray(dic): 
    scale = 0.3283 #pixel/micron
    a = np.zeros(shape=(26,26)) 
    
    for roi, value in dic.items():
        (x,y) = (int(rois[roi]["left"]/(scale*30)), int(rois[roi]["top"]/(scale*30)))
        a[x][y] = value
        
    return a


def PlotVTA(Name, Array, SavePath):
    
    cmap=plt.get_cmap("hot")
    
    fig, ax = plt.subplots()
    im = ax.imshow(Array, cmap=cmap, vmin=0, vmax=1)
    
    locs, labels = plt.xticks()
    labels = [int(item)*30 for item in locs]
    plt.xticks(locs, labels)
    locs, labels = plt.yticks()
    labels = [int(item)*30 for item in locs]
    plt.yticks(locs, labels)
    
    ax.set_xlim(-0.5, 25.5)
    ax.set_ylim(25.5, -0.5)
    ax.set_xlabel("$\mu$m")
    ax.set_ylabel("$\mu$m")
    ax.set_title(Name)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cax.set_ylabel("Stimulation Eficacy")
    plt.colorbar(im, cax=cax)
    
    plt.tight_layout()
    plt.savefig(SavePath + "\\Figures\\" + Name + ".svg")
    plt.savefig(SavePath + "\\Figures\\" + Name + ".png", dpi = 300)
    


#%%



Path = 'SLCM\\5th Culture\\HEX08\DIV16\\VTA\\Current'


My_Dicts = []

for path, subdirs, files in os.walk(Path):
    for file in files:   
        if "FluoPeaks.p" in file:
            if "png" not in file:
                My_Dicts.append(file)

Dicts_of_Dicts = {}

for dic in My_Dicts:
    Dicts_of_Dicts[dic] = pickle.load(open(Path+"\\"+dic, "rb"))
    
    
ROI_Name = "\\RoiSet_BackUp.zip"

rois = read_roi_zip(Path+ROI_Name)



# Grid 780x780/30




#%% Calculate stimulation eficacy arrays

ArraysToPlot = {}

Current10uA_320us_01 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[0]]))
Current10uA_320us_02 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[1]]))
Current10uA_320us = np.mean( np.array([ Current10uA_320us_01, Current10uA_320us_02 ]), axis=0 )
ArraysToPlot["Current 10uA 320 us"] = Current10uA_320us

Current10uA_520us_01 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[2]]))
Current10uA_520us_02 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[3]]))
Current10uA_520us = np.mean( np.array([ Current10uA_520us_01, Current10uA_520us_02 ]), axis=0 )
ArraysToPlot["Current 10uA 520 us"] = Current10uA_520us

Current5uA_120us_01 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[4]]))
Current5uA_120us_02 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[5]]))
Current5uA_120us = np.mean( np.array([ Current5uA_120us_01, Current5uA_120us_02 ]), axis=0 )
ArraysToPlot["Current 5uA 120 us"] = Current5uA_120us

Current5uA_320us = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[6]]))
ArraysToPlot["Current 5uA 320 us"] = Current5uA_320us

Current5uA_520us_01 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[7]]))
Current5uA_520us_02 = ROI_DicToArray(CalculateStimulationEficacy(Dicts_of_Dicts[My_Dicts[8]]))
Current5uA_520us = np.mean( np.array([ Current5uA_520us_01, Current5uA_520us_02 ]), axis=0 )
ArraysToPlot["Current 5uA 520 us"] = Current5uA_520us


#%% Plot chicha

for Name, VTA_Array in ArraysToPlot.items():
    PlotVTA(Name, VTA_Array, Path)

