# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 12:13:02 2018

@author: jdelacruz
"""

import os
import pandas as pnd
from datetime import datetime
import neo
import quantities as pq
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalProcess as RPro
import matplotlib.colors as mpcolors
import PhyREC.SignalAnalysis as Ran
from matplotlib import cm
from PhyREC.NeoInterface import NeoSegment, NeoSignal
import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
from scipy import signal
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats




def LoadAnalisys(FileName):
    
    Seg = NeoSegment()    
    
    date_format = "%d/%m/%Y" 
    
    SegAnnotations = {}    
    Data = pnd.read_csv(FileName, header=1, sep="\t", keep_default_na=False, encoding = "ISO-8859-1")
    
    SegAnnotations['AnimalNumber'] = Data.loc[:, 'Value'][9]
    SegAnnotations['Age'] = str(datetime.strptime("15/11/2018", date_format) - datetime.strptime(Data.loc[:, 'Value'][8].split(" ")[0], date_format)).split(" ")[0]
    if FileName.split("\\")[-1].split(".")[0].startswith("scoto"):
        ExperimentName = FileName.split("\\")[-1].split(".")[0][0:-1]
    else:
        ExperimentName = FileName.split("\\")[-1].split(".")[0]
    SegAnnotations["Experiment"] = ExperimentName
    
    Seg.Seg.annotate(**SegAnnotations)

    for header in Data.columns:
        if header.startswith("Chan 1") or header.startswith("Chan 2"):
            Sig = neo.AnalogSignal(Data.loc[:, header][1:-1].astype(float)*1e-3,
                           name=header,
                           sampling_rate=1e3*pq.Hz,
                           t_start=-0.05*pq.s,
                           units='uV')
            if header.endswith(".1"):
                Step = 2
            elif header.endswith(".2"):
                Step = 3
            else:
                Step = 1
            SigAnnotations = {}    
            SigAnnotations["Stimulus"] = Data.loc[:,'cd.s/m\xb2'][Step]
            SigAnnotations.update(SegAnnotations)
            Sig.annotate(**SigAnnotations)
            Seg.AddSignal(Sig)
    return Seg

def InterpolatePSD(PSD, Ch, npoints, ax, color = 'k', linestyle = '-'):
        log_ff = np.log10(PSD[Ch]['ff'][1:].squeeze().magnitude)
        log_psd = np.log10(PSD[Ch]['psd'][1:].squeeze().magnitude)
        interpolation = interp1d(log_ff, log_psd)
        ff1 = np.linspace(min(log_ff), max(log_ff), num=npoints, endpoint=True)
        psd1 = interpolation(ff1)
        ff2 = 10**ff1
        psd2 = 10**psd1
        
        print(10**ff1)
        ax.loglog(ff2, psd2, color = color, LineWidth = 0.5, linestyle = linestyle)
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Power [$\mu$$V^2$/Hz]')

            
        return ff2, psd2
    
def InterpolateSlot(Wave, npoints, ax, color = 'k', linestyle = '-'):
        Time = np.linspace(Wave.Signal.t_start.magnitude, Wave.Signal.t_stop.magnitude, num = len(Wave.Signal))
        interpolation = interp1d(Time, Wave.Signal.flatten())
        Time1 = np.linspace(min(Time), max(Time), num=npoints, endpoint=True)
        SignalInt = interpolation(Time1)
        
        ax.plot(Time1, SignalInt, color = color, LineWidth = 0.5, linestyle = linestyle)
#        ax.set_xlabel('Time [T]')
#        ax.set_ylabel('Signal [V]')
            
        return SignalInt    

def AvgSignals(Sigs):
    acc = np.array([])
    for sig in Sigs:       
        si = sig.GetSignal(None).magnitude
        acc = np.hstack((acc, si)) if acc.size else si   
  
    return NeoSignal(np.mean(acc, axis=1),
                     name='Mean',
                     sampling_rate=sig.sampling_rate,
                     t_start=sig.t_start,
                     units='uV')
    
def StdSignals(Sigs):
    acc = np.array([])
    for sig in Sigs:       
        si = sig.GetSignal(None).magnitude
        acc = np.hstack((acc, si)) if acc.size else si   
  
    return NeoSignal(np.std(acc, axis=1),
                     name='Std',
                     sampling_rate=sig.sampling_rate,
                     t_start=sig.t_start,
                     units='uV')
    
def SearchData(Data, Filter):
    Seg = []
    for Rat, Dats in Data.items():
        for dat in Dats:
            for sig in dat.Signals():
                shared_items = {k: sig.annotations[k] for k in sig.annotations if k in Filter and sig.annotations[k] == Filter[k]}
                if len(shared_items) == len(Filter):
                    Seg.append(sig)
                    print(sig.annotations["AnimalNumber"] + " added!")
    SigsGrph = [] # derecha
    SigsGold = []  #izquirda
    for sig in Seg:
        if sig.name.startswith('Chan 1'):
            SigsGrph.append(sig)
        elif sig.name.startswith('Chan 2'):
            SigsGold.append(sig)
    return SigsGrph, SigsGold

def CreateFilter(Age, Experiment, Stimulus):
    return     {"Age":Age,
              "Experiment":Experiment,
              "Stimulus": Stimulus}

def DrawBarScaleDV2(Ax, Location='Bottom Left',
                 xsize=None, ysize=None, xoff=0.1, yoff=0.1,
                 xlabelpad=-0.04, ylabelpad=0.04,
                 xunit='s', yunit='$\mu$V', LineWidth=5, Color='k',
                 FontSize=None):

    # calculate length of the bars
    xmin, xmax, ymin, ymax = Ax.axis()
    AxTrans = Ax.transAxes
    if xsize is None:
        xsize = (xmax - xmin)/5
        xsize = int(np.round(xsize, 0))
    if ysize is None:
        ysize = (ymax - ymin)/5
        ysize = int(np.round(ysize, 0))
    xlen = 1/((xmax - xmin)/xsize)  # length in axes coord
    ylen = 1/((ymax - ymin)/ysize)

    # calculate locations
    if Location == 'Bottom Rigth':
        xoff = 1 - xoff
        ylabelpad = - ylabelpad
        xlen = - xlen
    elif Location == 'Top Left':
        yoff = 1 - yoff
        ylen = - ylen
        xlabelpad = -xlabelpad
    elif Location == 'Top Rigth':
        xoff = 1 - xoff
        ylabelpad = - ylabelpad
        xlen = - xlen
        yoff = 1 - yoff
        ylen = - ylen
        xlabelpad = -xlabelpad
    xdraw = xoff + xlen
    ydraw = yoff + ylen

    # Draw lines
    Ax.hlines(yoff, xoff, xdraw,
              Color,
              linewidth=LineWidth,
              transform=AxTrans,
              clip_on=False)

    Ax.text(xoff + xlen/2,
            yoff + xlabelpad,
            str(float(xsize)) + ' ' + xunit,###•here
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=FontSize,
            transform=AxTrans)

    Ax.vlines(xoff+xlen, yoff, ydraw,
              Color,
              linewidth=LineWidth,
              transform=AxTrans,
              clip_on=False)

    Ax.text(xoff + xlen + ylabelpad,
            yoff + ylen/2,
            str(int(ysize)) + ' ' + yunit,
            horizontalalignment='center',
            verticalalignment='center',
            rotation='vertical',
            fontsize=FontSize,
            transform=AxTrans) 
    

#%%   
if __name__ == '__main__': 

    date_format = "%d/%m/%Y"   
    root = "Macro"
    
    MyData = {}
    Rats = [ ]
    for path, subdirs, files in os.walk(root):
        if len(path.split('\\')) < 3:
            continue
        Rats.append(path.split('\\')[-1])    
    for r in set(Rats):
        MyData[r] = []
            
    
    for path, subdirs, files in os.walk(root):
        for File in files:
            if File.endswith('.pdf'):
                continue
            Rat = path.split('\\')[-1]
            FilePath = path + "\\" + File
            print(FilePath + " imported!")
            data = LoadAnalisys(FilePath)
                
            MyData[Rat].append(data)


#%%  
    PeakFilter = [(148, 152), (248, 252), (298, 302), (348, 352)]

    Fband = (100,150)
    F1 = [{'function': RPro.Filter, 'args': {'Type':'bandpass',
                                                    'Order':2,
                                                    'Freqs':(Fband[0], Fband[1])}},
#         {'function': RPro.RemoveDC, 'args': {}}
         ]
    
#    for item in PeakFilter:
#        F1.append({'function': RPro.Filter, 'args': {'Type':'bandstop',
#                                                  'Order':2,
#                                                  'Freqs':(item[0], item[1])}})

    
#%%
    
    SigsGrph = []
    SigsGold = []
    SlotsGraphene = []
    SlotsGold = []
    
    colors = cm.rainbow(np.linspace(0, 1, 5))

    monocolors =["#8000ff",
                "#9932ff",
                "#b266ff",
                "#bf7fff",
                "#d8b2ff"]
    
    LineWidth = 1

#    Grph, Gold = SearchData(MyData, CreateFilter("27", "scoto", "10"))
#    SigsGrph.append(Grph), SigsGold.append(Gold)
#
#    Grph, Gold = SearchData(MyData, CreateFilter("61", "scoto", "10"))    
#    SigsGrph.append(Grph), SigsGold.append(Gold)
#    
#    Grph, Gold = SearchData(MyData, CreateFilter("90", "scoto", "10"))    
#    SigsGrph.append(Grph), SigsGold.append(Gold)
#    
#    Grph, Gold = SearchData(MyData, CreateFilter("165", "scoto", "10"))
#    SigsGrph.append(Grph), SigsGold.append(Gold)
#    
#    Grph, Gold = SearchData(MyData, CreateFilter("221", "scoto", "10"))
#    SigsGrph.append(Grph), SigsGold.append(Gold)
    
    Grph, Gold = SearchData(MyData, CreateFilter("27", "scoto", "10"))
    SigsGrph.append(Grph), SigsGold.append(Gold)

    # Grph, Gold = SearchData(MyData, CreateFilter("27", "scoto", "0.3"))    
    # SigsGrph.append(Grph), SigsGold.append(Gold)
    
    # Grph, Gold = SearchData(MyData, CreateFilter("27", "scoto", "0.03"))    
    # SigsGrph.append(Grph), SigsGold.append(Gold)
    
    # Grph, Gold = SearchData(MyData, CreateFilter("27", "scoto", "0.0001"))
    # SigsGrph.append(Grph), SigsGold.append(Gold)
    
    # Grph, Gold = SearchData(MyData, CreateFilter("27", "scoto", "3E-6"))
    # SigsGrph.append(Grph), SigsGold.append(Gold)
    
    
#%%
#    plt.close('all')
    # plt.rcParams.update({'font.size': 11})   
    # Fig, Ax = plt.subplots(1,2, sharex= True, sharey = True, figsize=(9,4.5))
    # Fig.subplots_adjust(wspace=0, left=0, right=0.85, bottom = 0.03, top = 0.903)
    # Ax[1].set_xlim((-0.075,0.075))
    # Ax[0].set_xlim((-0.05,0.1))

    
    StatGold = []
    StatGraphene = []
    
    SlotsGraphene = []
    SlotsGold = []
    
    # Offset = np.array([0,
    # -85,
    # -135,
    # -220,
    # -280
    # ])*pq.uV
    
    for sl, item in enumerate(SigsGrph):
        for sig in item:
            sig.ProcessChain = F1
    #        SlotsGraphene.append(Rplt.WaveSlot(sig,
    #                                           Position=0,
    #                                           Alpha=0.5,
    #                                           ))
            name = sig.annotations["Age"] + " days"
            SlotsGraphene.append(Rplt.WaveSlot(sig,
                                               Position=0,
                                               Alpha=1,
                                               LineWidth = LineWidth,
                                               DispName= name,
                                               Ax = 0,
                                               Ylim= (-20,20),
                                               # Color=colors[sl]
                                               ))
            # Ax[0].plot(AvgSignals(item).times, AvgSignals(item)+Offset[sl], color = monocolors[sl], label = None)
            StatGraphene.append(sig)
    
    for sl, item in enumerate(SigsGold):
        for sig in item:
            sig.ProcessChain = F1
    #        SlotsGraphene.append(Rplt.WaveSlot(sig,
    #                                           Position=0,
    #                                           Alpha=0.5,
#                                                ))
            name=sig.annotations["Age"] #+ " days"                                       
            SlotsGold.append(Rplt.WaveSlot(sig,
                                               Position=0,
                                               Alpha=1,
                                               LineWidth = LineWidth,
                                               DispName= name,
                                               Ax = 0,
                                               Ylim= (-20,20),
                                               # Color=colors[sl]
                                               ))
            # Ax[1].plot(AvgSignals(item).times, AvgSignals(item)+Offset[sl], color = monocolors[sl], label = name)
            # Ax[1].text(0.52, (-20+Offset[sl].magnitude), name)
            # Ax[1].text(0.52, 200, "Days")
            StatGold.append(sig)
#
    # Ax[0].axis('off')
    # Ax[1].axis('off')
    # DrawBarScaleDV2(Ax=Ax[1], Location='Top Left', Color = 'black',
    #                      xsize=0.02, ysize=75,LineWidth=2,
    #                      xlabelpad=-0.05, ylabelpad = 0.08,
    #                      xoff = 0.55, yoff = 0.1)

    
#%% Remove light artefact
    
    # plt.rcParams.update({'font.size': 18})   
    # Fig, Ax = plt.subplots(1,2, sharex= True, sharey = True, figsize=(9,4.5))
    # Fig.subplots_adjust(wspace=0.03, left=0.02, right=0.85, bottom = 0.03, top = 0.903)
    # Ax[1].set_xlim((0.008,0.075))
    # Ax[0].set_xlim((0.008,0.075))
    
    ReplaceWindow = (36,67)
        
    ArtefactGrapheneSlot = []
        
    # for ages in SigsGrph[0:3]:
    # dummy = []
    for item in StatGraphene:              
        Cutted = item.GetSignal(None)
        ValueToReplace = np.average(item.GetSignal(None)[ReplaceWindow[0]:ReplaceWindow[1]] )
        Cutted[slice(ReplaceWindow[0], ReplaceWindow[1])] = ValueToReplace
        # dummy.append(Cutted)
        
        ArtefactGrapheneSlot.append(Cutted)
        
    ArtefactGoldSlot = []
            
    # for ages in SigsGold[0:3]:
    #     dummy = []
    for item in StatGold:              
        Cutted = item.GetSignal(None)
        ValueToReplace = np.average(item.GetSignal(None)[ReplaceWindow[0]:ReplaceWindow[1]] )
        Cutted[slice(ReplaceWindow[0], ReplaceWindow[1])] = ValueToReplace
        # dummy.append(Cutted)
        
        ArtefactGoldSlot.append(Cutted)

    
    # for sl, item in enumerate(ArtefactGrapheneSlot):
    #     name=item[0].annotations["Age"] #+ " days"                                       
    #     SlotsGold.append(Rplt.WaveSlot(AvgSignals(item),
    #                                        Position=0,
    #                                        Alpha=1,
    #                                        LineWidth = LineWidth,
    #                                        DispName= name,
    #                                        Ax = 0,
    #                                        Color=colors[sl]
    #                                        ))

    #     Ax[0].fill_between(AvgSignals(item).times, AvgSignals(item).squeeze() + StdSignals(item).squeeze()+Offset[sl],
    #       AvgSignals(item).squeeze() - StdSignals(item).squeeze()+Offset[sl], alpha = 0.37, color = colors[sl])
    #     Ax[0].plot(AvgSignals(item).times, AvgSignals(item)+Offset[sl], color = colors[sl], label = name, linewidth = 2.5)

    # NormalizationFactor = (161+618.0)/(150+597.0)

    # for sl, item in enumerate(ArtefactGoldSlot):
    #     name=item[0].annotations["Stimulus"] #+ " days"                                       
    #     SlotsGold.append(Rplt.WaveSlot(AvgSignals(item),
    #                                        Position=0,
    #                                        Alpha=1,
    #                                        LineWidth = LineWidth,
    #                                        DispName= name,
    #                                        Ax = 0,
    #                                        Color=colors[sl]
    #                                        ))

    #     Ax[1].fill_between(AvgSignals(item).times, NormalizationFactor*AvgSignals(item).squeeze() + NormalizationFactor*StdSignals(item).squeeze()+Offset[sl],
    #       NormalizationFactor*AvgSignals(item).squeeze() - NormalizationFactor*StdSignals(item).squeeze()+Offset[sl], alpha = 0.37, color = colors[sl])
    #     Ax[1].plot(AvgSignals(item).times, NormalizationFactor*AvgSignals(item).magnitude+Offset[sl].magnitude, color = colors[sl], label = name, linewidth = 2.5)
    #     Ax[1].text(0.078, (Offset[sl].magnitude), name)
    #     Ax[1].text(0.078, 30, "cds/m$^{2}$")

    # Ax[0].axis('off')
    # Ax[1].axis('off')
    # DrawBarScaleDV2(Ax=Ax[1], Location='Top Left', Color = 'black',
    #                      xsize=0.02, ysize=30,LineWidth=4,
    #                      xlabelpad=-0.05, ylabelpad = 0.08,
    #                      xoff = 0.55, yoff = 0.03)
    
#    Fig.savefig('Figure_3OPS_Normalized_MRS', dpi=2000, bbox_inches='tight') 
#%% Calculate significance

StatGr = []
StatGold = []

for item in ArtefactGrapheneSlot:
    StatGr.append(max(item))
    
for item in ArtefactGoldSlot:
    StatGold.append(max(item))

SignificanceTest = stats.ttest_ind(ArtefactGrapheneSlot, ArtefactGoldSlot, axis = 1)