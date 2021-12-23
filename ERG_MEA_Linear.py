# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:56:37 2020

@author: jdelacruz
"""

import glob
import OpenEphys
import matplotlib.pyplot as plt
import quantities as pq
from PhyREC.NeoInterface import NeoSegment, NeoSignal
import neo
import PhyREC.SignalProcess as RPro
import numpy as np
import glob
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalAnalysis as Rana
from matplotlib import cm
import matplotlib.colors as mpcolors
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate



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


def get_key(my_dict, val): 
    for key, value in my_dict.items(): 
         if val == value: 
             return key 
  
    return "key doesn't exist"
        

#%%

if __name__ == '__main__':


    # FilePath = '2019-11-28_19-54-40\\*.continuous'
    FilePath = '20200728_open_ephys_erg_l6-5\\R9_2020-07-28_16-53-59\\*.continuous'
    
    # InChans = ("CH15", "CH24", "CH61")
    
    ExcludeChans = []
    
    InChans = []
    
    for i in range(16):
        if i+1 in ExcludeChans:
            continue
        else:
            InChans.append("CH{}".format(i+17))
    
    Seg = LoadOpenEphysContinuous(FilePath, IncludeCh=InChans)
    
    
    
    cmap = cm.ScalarMappable(mpcolors.Normalize(vmin=0, vmax=len(Seg.SigNames.keys())), cm.jet)
    colors = [cmap.to_rgba(i) for i in range(len(Seg.SigNames.keys()))]
    
    
    #TWind = (10*pq.s, 20*pq.s)
    Fband = (.5, 150)
    F1 = [{'function': RPro.Filter, 'args': {'Type':'bandpass',
                                                    'Order':2,
                                                    'Freqs':(Fband[0], Fband[1])}},
          {'function': RPro.Filter, 'args': {'Type':'bandstop',
                                                    'Order':2,
                                                    'Freqs':(48, 52)}},
         {'function': RPro.RemoveDC, 'args': {}}
         ]
    
    
    Slots = []
    
    for ip, sig in enumerate(Seg.Signals()):
        sig.ProcessChain = F1
        Slots.append(Rplt.WaveSlot(sig,
                                   Units='uV',
                                   Color=colors[ip],
                                   Alpha=1,
    #                               DispName=sn #+ '-' + Labels[sn],
                                   ))
            
#%%
    Splt = Rplt.PlotSlots(Slots)
    Splt.PlotChannels(None, Units='uV')
    
    
#%%

AvgSeg = NeoSegment()



ChannToElect = {
                "CH1": "E24", 
                 "CH2": "E25",
                 "CH3": "E26",
                 "CH4": "E27",
                 "CH5": "E28",
                 "CH6": "E29",
                 "CH7": "E30",
                 "CH8": "E31",
                 "CH10": "E17",
                 "CH11": "E18",
                 "CH12": "E19",
                 "CH13": "E20",
                 "CH14": "E21",
                 "CH15": "E22",
                 "CH16": "E23",
                 "CH17": "E9",
                 "CH18": "E10",
                 "CH19": "E11",
                 "CH20": "E12",
                 "CH21": "E13",
                 "CH22": "E14",
                 "CH23": "E15",
                 "CH24": "E16",
                 "CH25": "E1",
                 "CH26": "E2",
                 "CH27": "E3",
                 "CH28": "E4",
                 "CH29": "E5",
                 "CH30": "E6",
                 "CH31": "E7",
                 "CH32": "E8"
                 }

Twind = (-.05, .35)*pq.s
Trig = LoadOpenEphysEvents(FilePath)
TrigTimes = Trig["timestamps"][0::2]/float(Trig["header"]["sampleRate"])*pq.s
CorrectedTrigTimes = TrigTimes - 459.524*pq.s

for sig in Slots:    
    AvgSig = neo.AnalogSignal(sig.CalcAvarage(TimesEvent=CorrectedTrigTimes,
                                              TimeAvg=Twind),
                                t_start=Twind[0],
                                sampling_rate=sig.GetSignal(None).sampling_rate,
                                name=ChannToElect[sig.name]
                               )
    
    AvgSeg.AddSignal(AvgSig)

    

AvgSlots = []   
    
for ip, sig in enumerate(AvgSeg.Signals()):
    AvgSlots.append(Rplt.WaveSlot(sig,
                               Position=ip,
                                Alpha=0.5
                               ))

Splt2 = Rplt.PlotSlots(AvgSlots)
Splt2.PlotChannels(None, Units='uV')




#%%n Calculate bWave // Takes a long processing time

BWaveDict = {}
BWaveMean = [] 

for sig in AvgSeg.signames:
    Max = max(AvgSeg.GetSignal(sig))
    BWaveDict[sig] = Max
    BWaveMean.append(Max) 
    
BWaveMean = np.mean(BWaveMean)


#%%

bWaveList = []
fig = plt.figure()


for key, item in BWaveDict.items():
    plt.scatter(float(key.split("E")[1]), item, color = "b")
    bWaveList.append(item)
    print(float(key.split("E")[1]), item)
    


plt.plot(np.linspace(1, 7), np.linspace(min(bWaveList), min(bWaveList)), color = "r", linewidth = 4)



                 
#%% Plot differences vs average

SigToAvg = []
ExcludeFromAvg = ["E4", "E15"]

for sig in AvgSeg.signames:
    if sig in ExcludeFromAvg:
        continue
    else:
        SigToAvg.append(AvgSeg.GetSignal(sig))
                 
AvgWave = np.mean(SigToAvg, axis = 0)*pq.uV

SubtractSeg = NeoSegment()

for sig in AvgSeg.Signals():
     SubtractSeg.AddSignal(sig.duplicate_with_new_data(signal = sig - AvgWave))



fig, axs = plt.subplots(4, 4, sharex=True, sharey=True)
axs = axs.ravel()

SortedList = [
                "E1",
                "E2",
                "E3",
                "E4",
                "E5",
                "E6",
                "E7",
                "E8",
                "E9",
                "E10", 
                "E11",
                "E12",
                "E13",
                "E14",
                "E15",
                "E16"]



for sl, ch in enumerate(SortedList):
    axs[sl].plot(AvgSeg.GetSignal(ch).times, AvgSeg.GetSignal(ch))
    axs[sl].plot(SubtractSeg.GetSignal(ch).times, SubtractSeg.GetSignal(ch))
    axs[sl].set_title(ch)
    
    
        