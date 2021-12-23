# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:56:37 2020

@author: jdelacruz
"""

import matplotlib.pyplot as plt
import McsPy.McsData as McsData
import quantities as pq
from PhyREC.NeoInterface import NeoSegment
import neo
import PhyREC.SignalProcess as RPro
import numpy as np
import glob
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalAnalysis as Rana
import pickle

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
                            
         
                
# Uncomment to include trigger import
        
        Trigger = Rec.event_streams[0].event_entity[2].get_event_timestamps()[0]*pq.us
        TriggerEvent = neo.Event(times = Trigger, name = "EndPulse", units=Trigger.units)
        OutSeg.Seg.events.append(TriggerEvent)
        OutSeg.UpdateEventDict()
        
        
    return OutSeg
        

def CountValues(my_dict):
    items = 0
    for value in my_dict.values():
        items = items + len(value)
        
    return items
        

#%% Load Files

plt.close("all")


# Path = '8th Culture\\HEX03\DIV12\\'

Path = 'SeanRecords\\Output\\'



File = '2021-11-12T13-53-54Hex01_6_D-00105.h5'

FigSavePath = "Histograms\\Voltage\\"

InChans = [
            "14",
            # "33"
           ]
Slots = []

Seg = LoadMCS(Path+File, InChans = InChans)


F1 = [{'function': RPro.Filter, 'args': {'Type':'highpass',
                                        'Order':2,
                                        'Freqs':200}},
  {'function': RPro.RemoveDC, 'args': {}}
 ]    

LegendKwargs = {'fontsize': 'large',
                'ncol': 5,
                'loc': 'upper right',
                'frameon': False}

#Create Waveslot

for ip, sig in enumerate(Seg.Signals()):
    sig.ProcessChain = F1
    Slots.append(Rplt.WaveSlot(sig,
                               Position=ip,
                                Alpha=0.5
                               ))

Splt = Rplt.PlotSlots(Slots)
Splt.PlotChannels(None, Units='uV')
Splt.AddLegend(**LegendKwargs)

# plt.savefig(Path + FigSavePath + File.split(".")[0] + ".svg")
# plt.savefig(Path + FigSavePath + File.split(".")[0] + ".png", dpi = 300)

#%%

Pulses = Seg.GetEventTimes("EndPulse")

PulseBuffer = 10*pq.ms
WindowLength = 1*pq.s

PulsesWindow = []

for item in Pulses:
    PulsesWindow.append((item + PulseBuffer, item + WindowLength))
    
#slice signal

SliceSeg = NeoSegment()

for sl, window in enumerate(PulsesWindow):
        
        
    SlicedSignal = neo.AnalogSignal(Seg.GetSignal(InChans[0]).GetSignal(window),
                                    t_start = 0*pq.s,
                                    sampling_rate=Seg.GetSignal(InChans[0]).sampling_rate,
                                    name = "Slice " + str(sl)
                                    )
    SliceSeg.AddSignal(SlicedSignal)
                                    
SlicedSlots = []

for ip, sig in enumerate(SliceSeg.Signals()):
    SlicedSlots.append(Rplt.WaveSlot(sig,
                               Position=ip,
                                Alpha=1,
                                label = sig.name
                               ))

SpltSliced = Rplt.PlotSlots(SlicedSlots)
SpltSliced.PlotChannels(None, Units='uV')




#%% Threshold detection

Spikes = {}

Threshold_Slice = "Slice 0"
Threshold_Window = (.1, .175)*pq.s

# Threshold = -1*np.std(SliceSeg.GetSignal(Threshold_Slice).GetSignal(Threshold_Window)).magnitude*5e6*pq.uV

Threshold = -15*pq.uV

RefractoryPeriod = 10*pq.ms


for sig in SliceSeg.Signals():  
    Spikes[sig.name] = Rana.threshold_detection(sig, Threshold, sign = "below", RelaxTime = RefractoryPeriod)
 

Number_spikes = CountValues(Spikes)

print(Number_spikes, len(Spikes.keys()))

for ax in SpltSliced.Axs:
    x_lim = ax.get_xlim()*pq.s
    ax.hlines(Threshold, x_lim[0], x_lim[1], color = "r", Alpha = 0.5)
    
    
for sl, value in enumerate(Spikes.values()):
    ax = SpltSliced.Axs[sl] 
    y_lim = ax.get_ylim()*pq.uV
    for event in value:
        ax.vlines(event, y_lim[0], y_lim[1], color = "r", Alpha = 0.5)

SpltSliced.Axs[0].text(1, 0, str(Number_spikes) + " spikes", fontsize=12)   
SpltSliced.Axs[1].text(1, 0, str(Threshold) + " Threshold", fontsize=12)


# plt.savefig(Path + FigSavePath + File.split(".")[0] + "_Spikes_.svg")
# plt.savefig(Path + FigSavePath + File.split(".")[0] + "_Spikes_.png", dpi = 300)

pickle.dump(Spikes, open(Path + FigSavePath + File.split(".")[0] + "_Ch_" + InChans[0] + ".p", "wb" ))


#%%

ax = Splt.Axs[0]
y_lim = ax.get_ylim()*pq.s

for item in Pulses:
    ax.vlines(item.rescale("s"), y_lim[0], y_lim[1], color = "r", Alpha = 0.5, linestyle = "dashed")
    
#%%  Plot scatter spike timestamps

Fig, AxScat = plt.subplots()


for sl, SpikeTimes in enumerate(Spikes.values()):
    if sl == 5:
        continue
    for spike in SpikeTimes:
        AxScat.scatter(spike.rescale("ms"), sl, color = "black")
        
AxScat.set_xlabel("Time (ms)")
AxScat.set_ylabel("Trials")
AxScat.set_ylim(0, 15)
        

# Fig.savefig(Path + FigSavePath + File.split(".")[0] + "_Ch" + InChans[0] +".svg")
# Fig.savefig(Path + FigSavePath + File.split(".")[0] + "_Ch" + InChans[0] + ".png", dpi = 300)



#%% 


