# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 12:44:52 2021

@author: jdelacruz
"""


import glob
import OpenEphys
#from openephys.OpenEphys import load as EphyLoad
from PhyREC.NeoInterface import NeoSegment, NeoSignal
import neo
import quantities as pq
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalProcess as RPro
import SignalAnalysis as Ran
from matplotlib import cm
import matplotlib.colors as mpcolors
import matplotlib.pyplot as plt 
import numpy as np
from scipy.interpolate import interp1d
import matplotlib as matplotlib
from scipy import stats
import pickle

Ch = "Así no se hace"

def LoadOpenEphysContinuous(FilePath, IncludeCh=None):
    FileNames = glob.glob(FilePath)
    MySeg = NeoSegment()
    for File in FileNames:
        Data = OpenEphys.load(File)
        ChName=Data['header']['channel'].replace("'",'')
        if ChName.startswith('AUX'):
            continue
        if IncludeCh is not None:
            if ChName not in IncludeCh:
                continue

        Sig = neo.AnalogSignal(Data['data'],
                         name=ChName,
                         sampling_rate=float(Data['header']['sampleRate'])*pq.Hz,
                         units='uV',
                         file_origin=File)
        MySeg.AddSignal(Sig)
        
        
    return MySeg

def LoadOpenEphysEvents(FilePath):
    FileName = FilePath.split("*")[0] + "all_channels.events"
    Trigger = OpenEphys.loadEvents(FileName)
    
    return Trigger

def InterpolatePSD(PSD, Ch, interpmin, interpmax, npoints, ax, color = 'k', linestyle = '-', LineWidth = 2, label = str(Ch)):
        log_ff = np.log10(PSD[Ch]['ff'][1:].squeeze().magnitude)
        log_psd = np.log10(PSD[Ch]['psd'][1:].squeeze().magnitude)
        interpolation = interp1d(log_ff, log_psd, fill_value="extrapolate")
        ff1 = np.linspace(np.log10(interpmin), np.log10(interpmax), num=npoints, endpoint=False)
        psd1 = interpolation(ff1)
        ff2 = 10**ff1
        psd2 = 10**psd1
        
        print(10**ff1)
        ax.loglog(ff2, psd2, color = color, LineWidth = LineWidth, linestyle = linestyle, label = label)
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Power ($\mu$$V^2$/Hz)')
        ax.legend(loc='upper right')
            
        return ff2, psd2
    
#%%  
PeakFilter = [(49.8, 50.2), (59.8, 60.2), (99.8, 100.2), (119.8, 120.2)]

Fband = (5,6000)
F1 = [{'function': RPro.Filter, 'args': {'Type':'bandpass',
                                                'Order':2,
                                                'Freqs':(Fband[0], Fband[1])}},
        {'function': RPro.RemoveDC, 'args': {}}
     ]

for item in PeakFilter:
    F1.append({'function': RPro.Filter, 'args': {'Type':'bandstop',
                                              'Order':2,
                                              'Freqs':(item[0], item[1])}})

#%%  


FlexiblePath = "RigidNoise_Thesis\\*.continuous"
InChans = (
              "CH2",
             'CH3',
             'CH4',
             'CH5',
           #  "CH24",
           )
FlexSeg = LoadOpenEphysContinuous(FlexiblePath)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


FlexSlots = []

for ip, sig in enumerate(FlexSeg.Signals()):
    if sig.name == "CH24":
        sig.ProcessChain = F1
        FlexSlots.append(Rplt.WaveSlot(sig.GetSignal((0,16)*pq.s),
                                   Units='uV',
                                   Color=colors[ip],
                                   Alpha=1,
                                   # label= Labels[sig.name],
                                   ))
    else:
        sig.ProcessChain = F1
        FlexSlots.append(Rplt.WaveSlot(sig,
                                   Units='uV',
                                   Color=colors[ip],
                                   Alpha=1,
                                   # label= Labels[sig.name],
                                   ))

FlexSplt = Rplt.PlotSlots(FlexSlots)
FlexSplt.PlotChannels(None)

PSD = Ran.GetNoise(FlexSlots)

Legend_Dict = {"CH5": "50 $\mu$m",
               'CH3': '25 $\mu$m',
               'CH4': '100 $\mu$m',
               "CH2": "Floor Noise",
               "CH24": "Recording PSD"}


plt.rcParams.update({'font.size': 14})   
fig, ax = plt.subplots(figsize=(7, 5))
ax.tick_params(direction="in", which = 'both')

npoints = 20
InterpRange = (7, 6000)

for sl, Ch in enumerate(PSD.keys()):
    
    InterpolatePSD(PSD, Ch, InterpRange[0], InterpRange[1], npoints, ax, colors[sl], LineWidth=1.5,
                    label = Legend_Dict[Ch]
                   )

ax.legend()