# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 11:36:58 2018

@author: aemdlabs
"""

import glob
import OpenEphys
from PhyREC.NeoInterface import NeoSegment, NeoSignal
import neo
import quantities as pq
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalProcess as RPro
import PhyREC.SignalAnalysis as Ran
from matplotlib import cm
import matplotlib.colors as mpcolors
import matplotlib.pyplot as plt


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


FilePath = 'MEA1\\1to32\\2020-03-06_21-26-14\\*.continuous'

InChans = ("CH1")

#InChans = []
#ExcludeChans = [12, 13]
#
#for i in range(2):
#    if i+1 in ExcludeChans:
#        continue
#    InChans.append("CH{}".format(i+1))
#
Seg = LoadOpenEphysContinuous(FilePath, IncludeCh=InChans)

cmap = cm.ScalarMappable(mpcolors.Normalize(vmin=0, vmax=len(Seg.SigNames.keys())), cm.jet)
colors = [cmap.to_rgba(i) for i in range(len(Seg.SigNames.keys()))]


#TWind = (10*pq.s, 20*pq.s)
Fband = (300, 6000)
F1 = [{'function': RPro.Filter, 'args': {'Type':'bandpass',
                                                'Order':2,
                                                'Freqs':(Fband[0], Fband[1])}},
     {'function': RPro.RemoveDC, 'args': {}}
     ]

#%%

Slots = []
for ip, sig in enumerate(Seg.Signals()):
#    sig.ProcessChain = F1
    Slots.append(Rplt.WaveSlot(sig,
                               Units='uV',
                               Color=colors[ip],
                               Alpha=1,
#                               DispName=sn #+ '-' + Labels[sn],
                               ))

#Splt = Rplt.PlotSlots(Slots,
#                      ShowNameOn='Legend')
#
#Splt.PlotChannels(None)

#%%

#Ran.PlotPSD(Slots,
#            Time=None,
#            FMin=None,
#            scaling='density')

#
SpikeTimes = {}

for signal in Slots:
    SpikeTimes[signal.name] = Ran.threshold_detection(signal.GetSignal(None), -50, sign = "below", RelaxTime = 1*pq.s)
    
#%%
#
Spike = {}

for key, val in SpikeTimes.items():
    
    if len(val)<1:
        continue
    Spike[key]=val
    
#%%

Seg1 = LoadOpenEphysContinuous(FilePath, IncludeCh=Spike.keys())

Slots1 = []
for ip, sig in enumerate(Seg1.Signals()):
#    sig.ProcessChain = F1
    Slots1.append(Rplt.WaveSlot(sig,
                               Units='uV',
                               Color=colors[ip],
                               Alpha=1,
#                               DispName=sn #+ '-' + Labels[sn],
                               ))




Splt2 = Rplt.PlotSlots(Slots1)
#Splt2.PlotEventAvarage((-0.001*pq.s, 0.005*pq.ms),
#                       TimesEvent=SpikeTimes["CH19"],
#                       PlotTrials=True,
#                       TrialsAlpha=0.5)
