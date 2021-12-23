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



def LoadAnalisys(FileName):
    
    Seg = NeoSegment()    
    
    date_format = "%d/%m/%Y" 
    
    SegAnnotations = {}    
    Data = pnd.read_csv(FileName, header=1, sep="\t", keep_default_na=False, engine = "python")
    
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

def AvgSignals(Sigs, name = "Mean"):
    acc = np.array([])
    for sig in Sigs:       
        si = sig.GetSignal(None).magnitude
        acc = np.hstack((acc, si)) if acc.size else si   
  
    return NeoSignal(np.mean(acc, axis=1),
                     name=name,
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
    
def specDV(sig, Ax):
    Fres = 100
    TimeRes = 0.0001
    Fmax = 1e3
    Fmin = Fres
    nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/Fres))+1))
    Ts = sig.sampling_period.magnitude
    noverlap = int((Ts*nFFT - TimeRes)/Ts)
    
    f, t, Sxx = signal.spectrogram(sig, sig.sampling_rate,
                                   window='hanning',
                                   nperseg=nFFT,
                                   noverlap=noverlap,
                                   scaling='spectrum',
                                   axis=0)
    finds = np.where((Fmin < f) & (f < Fmax))[0][1:]
    r, g, c = Sxx.shape
    data = Sxx.reshape((r, c))[finds][:]
    x = t #+ sig.t_start.magnitude
    y = f[finds].magnitude/1000
    MaxPSD = data.max()
    MinPSD = 1e-6
    Norm = mpcolors.LogNorm(MinPSD, MaxPSD)
    img = Ax.imshow(data,
             cmap='viridis',
             norm=Norm,
             interpolation='quadric',
             origin='lower',
             aspect='auto',
             extent=(np.min(x), np.max(x), np.min(y), np.max(y))
             )
#    fig1, ax1 = plt.subplots()
#    cbar = plt.colorbar(img)
#    cbar.ax1.tick_params()
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
            if "flickers" in File:
                continue
            if "OPS" in File:
                continue
            Rat = path.split('\\')[-1]
            FilePath = path + "\\" + File
            print(FilePath + " imported!")
            data = LoadAnalisys(FilePath)
                
            MyData[Rat].append(data)


#%%  
    PeakFilter = [(148, 152), (248, 252), (298, 302), (348, 352)]

    Fband = (0.1,150)
    F1 = [{'function': RPro.Filter, 'args': {'Type':'bandpass',
                                                    'Order':2,
                                                    'Freqs':(Fband[0], Fband[1])}},
#         {'function': RPro.RemoveDC, 'args': {}}
         ]
    
    for item in PeakFilter:
        F1.append({'function': RPro.Filter, 'args': {'Type':'bandstop',
                                                  'Order':2,
                                                  'Freqs':(item[0], item[1])}})
        

    
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

    Grph, Gold = SearchData(MyData, CreateFilter("27", "scoto", "10"))
    SigsGrph.append(Grph), SigsGold.append(Gold)

    Grph, Gold = SearchData(MyData, CreateFilter("61", "scoto", "10"))    
    SigsGrph.append(Grph), SigsGold.append(Gold)
    
    Grph, Gold = SearchData(MyData, CreateFilter("90", "scoto", "10"))    
    SigsGrph.append(Grph), SigsGold.append(Gold)
    
    Grph, Gold = SearchData(MyData, CreateFilter("165", "scoto", "10"))
    SigsGrph.append(Grph), SigsGold.append(Gold)
    
    Grph, Gold = SearchData(MyData, CreateFilter("221", "scoto", "10"))
    SigsGrph.append(Grph), SigsGold.append(Gold)
    
#%%
#    plt.close('all')
    plt.rcParams.update({'font.size': 11})   
    Fig, Ax = plt.subplots(1,2, sharex= True, sharey = True, figsize=(9,4.5))
    Fig.subplots_adjust(wspace=0, left=0.075, right=0.85)
    Ax[1].set_xlim((-0.075,0.525))
    Ax[0].set_xlim((-0.075,0.525))

    
    AvgGold = []
    AvgGraphene = []
    
    SlotsGraphene = []
    SlotsGold = []
    
    Offset = np.array([0,
    -300,
    -450,
    -550,
    -650
    ])*pq.uV
    
    for sl, item in enumerate(SigsGrph):
        # for sig in item:
            # sig.ProcessChain = F1
    #        SlotsGraphene.append(Rplt.WaveSlot(sig,
    #                                           Position=0,
    #                                           Alpha=0.5,
    #                                           ))
        name = item[0].annotations["Age"] + " days"
        SlotsGraphene.append(Rplt.WaveSlot(AvgSignals(item),
                                           Position=0,
                                           Alpha=1,
                                           LineWidth = LineWidth,
                                           DispName= name,
                                           Ax = 0,
                                           Color=colors[sl]
                                           ))
        Ax[0].plot(AvgSignals(item).times, AvgSignals(item)+Offset[sl], color = monocolors[sl], label = None)
        AvgGraphene.append(AvgSignals(item))
    
    for sl, item in enumerate(SigsGold):
        # for sig in item:
            # sig.ProcessChain = F1
    #        SlotsGraphene.append(Rplt.WaveSlot(sig,
    #                                           Position=0,
    #                                           Alpha=0.5,
#                                                ))
        name=item[0].annotations["Age"] #+ " days"                                       
        SlotsGold.append(Rplt.WaveSlot(AvgSignals(item),
                                           Position=0,
                                           Alpha=1,
                                           LineWidth = LineWidth,
                                           DispName= name,
                                           Ax = 0,
                                           Color=colors[sl]
                                           ))
        Ax[1].plot(AvgSignals(item).times, AvgSignals(item)+Offset[sl], color = monocolors[sl], label = name)
        Ax[1].text(0.52, (-20+Offset[sl].magnitude), name)
        Ax[1].text(0.52, 200, "Days")
        AvgGold.append(AvgSignals(item))
#
    Ax[0].axis('off')
    Ax[1].axis('off')
    DrawBarScaleDV2(Ax=Ax[1], Location='Top Left', Color = 'black',
                         xsize=0.1, ysize=250,LineWidth=2,
                         xlabelpad=-0.05, ylabelpad = 0.08,
                         xoff = 0.55, yoff = 0.1)

    
#%% Remove light artefact
    
    ArtefactGrapheneSlot = []
    SlotsGold = []
        
    for ages in SigsGrph:
        dummy = []
        for item in ages:              
            Cutted = item.GetSignal(None)
            ValueToReplace = item.GetSignal(None)[46] 
            Cutted[slice(46, 54)] = ValueToReplace
            dummy.append(Cutted)
            
        ArtefactGrapheneSlot.append(dummy)
        

    plt.rcParams.update({'font.size': 18})   
    Fig, Ax = plt.subplots(1,2, sharex= True, sharey = True, figsize=(9,4.5))
    Fig.subplots_adjust(wspace=0, left=0, right=0.85, bottom = 0.03, top = 0.903)
    Ax[1].set_xlim((-0.075,0.475))
    Ax[0].set_xlim((-0.075,0.475))
    
    ArtefactAverageGrapheneSlot = []
    
    for sl, item in enumerate(ArtefactGrapheneSlot):
        name=item[0].annotations["Age"] #+ " days"                                       
        ArtefactAverageGrapheneSlot.append(Rplt.WaveSlot(AvgSignals(item, name),
                                           Position=0,
                                           Alpha=1,
                                           LineWidth = LineWidth,
                                           name= name,
                                           Ax = 0,
                                           Color=colors[sl]
                                           ))

        Ax[0].fill_between(AvgSignals(item).times, AvgSignals(item).squeeze() + StdSignals(item).squeeze()+Offset[sl],
          AvgSignals(item).squeeze() - StdSignals(item).squeeze()+Offset[sl], alpha = 0.37, color = colors[sl])
        Ax[0].plot(AvgSignals(item).times, AvgSignals(item)+Offset[sl], color = colors[sl], label = name, linewidth = 2)


    for sl, item in enumerate(SigsGold):
        name=item[0].annotations["Age"] #+ " days"                                       
        SlotsGold.append(Rplt.WaveSlot(AvgSignals(item, name),
                                           Position=0,
                                           Alpha=1,
                                           LineWidth = LineWidth,
                                           name = name,
                                           Ax = 0,
                                           Color=colors[sl]
                                           ))

        Ax[1].fill_between(AvgSignals(item).times, AvgSignals(item).squeeze() + StdSignals(item).squeeze()+Offset[sl],
          AvgSignals(item).squeeze() - StdSignals(item).squeeze()+Offset[sl], alpha = 0.37, color = colors[sl])
        Ax[1].plot(AvgSignals(item).times, AvgSignals(item)+Offset[sl], color = colors[sl], label = name, linewidth = 2)
        Ax[1].text(0.52, (-20+Offset[sl].magnitude), name)
        Ax[1].text(0.52, 200, "Days")


    Ax[0].axis('off')
    Ax[1].axis('off')
    DrawBarScaleDV2(Ax=Ax[1], Location='Top Left', Color = 'black',
                         xsize=0.1, ysize=250,LineWidth=4,
                         xlabelpad=-0.05, ylabelpad = 0.08,
                         xoff = 0.55, yoff = 0.1)
    
    # Fig.savefig('Figure_3b', dpi=2500, bbox_inches='tight') 




PSD_and_Noise_Au = Ran.GetNoise(SlotsGold)
PSD_and_Noise_Gr = Ran.GetNoise(ArtefactAverageGrapheneSlot)

FigPSD, AxPSD = plt.subplots()

AxPSD = plt.loglog(PSD_and_Noise_Au["27"]["ff"], PSD_and_Noise_Au["27"]["psd"], label = "Gold")
AxPSD = plt.loglog(PSD_and_Noise_Gr["27"]["ff"], PSD_and_Noise_Gr["27"]["psd"], label = "Graphene")
FigPSD.legend()


#%% Export pickle PSD_and_Noise_Au

import pickle

with open("PSD_Noise_Gr.pickle", "wb") as f:
    pickle.dump(PSD_and_Noise_Gr, f)
    


