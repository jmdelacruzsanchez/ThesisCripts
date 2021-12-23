# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 12:44:52 2021

@author: jdelacruz
"""


import glob
from matplotlib import cm
import matplotlib.colors as mpcolors
import matplotlib.pyplot as plt
import McsPy.McsData as McsData
import quantities as pq
from PhyREC.NeoInterface import NeoSegment, NeoSignal
import neo
import PhyREC.SignalProcess as RPro
import numpy as np
import glob
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalAnalysis as Ran
import pickle
from scipy.interpolate import interp1d
import matplotlib as matplotlib
from scipy import stats


def LoadMCS(FileSearch, InChans = None, Downsampling = 1):
    
    FileNames = glob.glob(FileSearch)
    


    for InFile in FileNames:
        
        print("File " + InFile + " added!")
        
        OutSeg = NeoSegment()
         
        
        Dat = McsData.RawData(InFile)
        Rec = Dat.recordings[0]
        
        NSamps = Rec.duration
        
        
        for AnaStrn, AnaStr in Rec.analog_streams.items():
            if len(AnaStr.channel_infos) == 1:
                continue
        
            for Chn, Chinfo in AnaStr.channel_infos.items():
                
                if InChans == None:
                    
                    print('Analog Stream ', Chinfo.label, Chinfo.sampling_frequency)
                    ChName = str(Chinfo.label)
            
                    Fs = Chinfo.sampling_frequency
                    Var, Unit = AnaStr.get_channel_in_range(Chn, 0, NSamps)
                    sig = neo.AnalogSignal(pq.Quantity(Var, Chinfo.info['Unit']),
                                            t_start=0*pq.s,
                                            sampling_rate=Fs.magnitude*pq.Hz,
                                            name=ChName)
                    DownSamp = RPro.DownSampling(sig, Downsampling)
                    OutSeg.AddSignal(DownSamp)
                    
                
                
                elif Chinfo.label in InChans:

                    
                    print('Analog Stream ', Chinfo.label, Chinfo.sampling_frequency)
                    ChName = str(Chinfo.label)
            
                    Fs = Chinfo.sampling_frequency
                    Var, Unit = AnaStr.get_channel_in_range(Chn, 0, NSamps)
                    sig = neo.AnalogSignal(pq.Quantity(Var, Chinfo.info['Unit']),
                                           t_start=0*pq.s,
                                           sampling_rate=Fs.magnitude*pq.Hz,
                                           name=ChName)
                    DownSamp = RPro.DownSampling(sig, Downsampling)
                    OutSeg.AddSignal(DownSamp)
                    
                else:
                    print("No channels loaded!")
                            
         
        
        
    return OutSeg

Ch = "Así no se hace"

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
        ax.set_ylabel('Power ($\mu V^{2}$/Hz)')
        ax.legend(loc='upper right')
            
        return ff2, psd2
    
#%%  
PeakFilter = [(49.8, 50.2), (59.8, 60.2), (99.8, 100.2), (119.8, 120.2)]

F1 = [{'function': RPro.Filter, 'args': {'Type':'highpass',
                                                'Order':2,
                                                'Freqs':5}},
        {'function': RPro.RemoveDC, 'args': {}}
     ]

# for item in PeakFilter:
#     F1.append({'function': RPro.Filter, 'args': {'Type':'bandstop',
#                                               'Order':2,
#                                               'Freqs':(item[0], item[1])}})

    
    
#%% Load Files

plt.close("all")
Path = 'rGONoise\\*.h5'
FileNames = glob.glob(Path)


BaseChans = [
            "64",
            "52",
            "75"
           ]

RecordChans = ["36"]

BaseSeg = LoadMCS(FileNames[1], InChans = BaseChans)
RecordSeg = LoadMCS(FileNames[0], InChans = RecordChans)

#%%  

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

LegendKwargs = {'fontsize': 'large',
                'ncol': 5,
                'loc': 'upper right',
                'frameon': False}

TWind = (0, 20)*pq.s
BaseSlots = []
RecordSlots = []


for ip, sig in enumerate(BaseSeg.Signals()):
    sig.ProcessChain = F1
    BaseSlots.append(Rplt.WaveSlot(sig,
                               Position=ip,
                                Alpha=0.5
                               ))


BaseSplt = Rplt.PlotSlots(BaseSlots)
BaseSplt.PlotChannels(TWind, Units='uV')
BaseSplt.AddLegend(**LegendKwargs)

for ip, sig in enumerate(RecordSeg.Signals()):
    sig.ProcessChain = F1
    RecordSlots.append(Rplt.WaveSlot(sig,
                               Position=ip,
                                Alpha=0.5
                               ))


RecordSplt = Rplt.PlotSlots(RecordSlots)
RecordSplt.PlotChannels(TWind, Units='uV')
RecordSplt.AddLegend(**LegendKwargs)


#%%
PSDBase = Ran.GetNoise(BaseSlots)
PSDRecord = Ran.GetNoise(RecordSlots)

#%%

plt.close("all")

Legend_Dict = {"75": "50 $\mu$m",
               '64': '25 $\mu$m',
               '52': '100 $\mu$m',
               "36": "Recording PSD"}


plt.rcParams.update({'font.size': 14})   
fig, ax = plt.subplots(figsize=(7, 5))
ax.tick_params(direction="in", which = 'both')

npoints = 20
InterpRange = (7, 4000)

for sl, Ch in enumerate(PSDBase.keys()):
    
    InterpolatePSD(PSDBase, Ch, InterpRange[0], InterpRange[1], npoints, ax, colors[sl], LineWidth=1.5,
                    label = Legend_Dict[Ch]
                    )

InterpolatePSD(PSDRecord, "36", InterpRange[0], InterpRange[1], npoints, ax, colors[3], LineWidth=1.5,
                label = Legend_Dict["36"])

ax.legend()


# for Ch, PSD in PSDBase.items():
#     ax.loglog(PSD["ff"], PSD["psd"], label = Ch)
    
# ax.legend()