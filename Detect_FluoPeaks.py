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
import pandas as pnd
from read_roi import read_roi_zip
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
                            
         
                
# Uncomment to include pulses triger
        
        Trigger = Rec.event_streams[0].event_entity[2].get_event_timestamps()[0]*pq.us
        TriggerEvent = neo.Event(times = Trigger, name = "EndPulse", units=Trigger.units)
        OutSeg.Seg.events.append(TriggerEvent)
        OutSeg.UpdateEventDict()
        
# Uncomment to include sync_trigger
        
        TriggerSync = Rec.event_streams[1].event_entity[0].get_event_timestamps()[0]*pq.us
        TriggerSyncEvent = neo.Event(times = TriggerSync, name = "SyncTrigger", units=Trigger.units)
        OutSeg.Seg.events.append(TriggerSyncEvent)
        OutSeg.UpdateEventDict()
        
        
    return OutSeg
        

def CountValues(my_dict):
    items = 0
    for value in my_dict.values():
        items = items + len(value)
        
    return items


def ReadFluo(FileName, Segment, t_start, sampling_rate = 1e3*pq.Hz):

    FluoRaw = pnd.read_csv(FileName, sep=",", keep_default_na=False)
    
    Heads=list(FluoRaw.columns.values)    
    
    print(len(Heads)/4)
    
    for item in Heads:
        
        if "Mean" not in item:
            continue
        else:
                    
            Sig = neo.AnalogSignal(FluoRaw.loc[:, item],
                                   name= item,
                                   sampling_rate = sampling_rate,
                                   t_start=t_start,
                                    units=pq.uV,
                                   annotate=["Area" + item[4:], FluoRaw.loc[0, "Area" + item[4:]]],
                                   file_origin=FileName)
        Segment.AddSignal(Sig)
        
    return Segment
        

#%% Load Files

plt.close()


Path = 'SLCM\\5th Culture\\HEX08\DIV16\\'


File = '2021-03-28T12-13-5310x_VideoRecording_Current_10uA_320us_02.h5'

FluoFile = "VTA\\Current\\10x_VideoRecording_Current10uA_320us_0_1p8579235.csv"

InChans = [
            "16",
            # "17"
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

t_start = Seg.GetEventTimes("SyncTrigger")[0]
Pulses = Seg.GetEventTimes("EndPulse").rescale("s")


# resamp = 1.74282313
# resamp = 1.85752527
# resamp = 1.84819922
# resamp = 1.72
# resamp = 1.89
resamp = 1.8579235

Seg = ReadFluo(Path + FluoFile, Seg, t_start = t_start, sampling_rate= resamp*pq.Hz)


#Create Waveslot

for ip, sig in enumerate(Seg.Signals()):
    if "Mean" not in sig.name:
        sig.ProcessChain = F1
    Slots.append(Rplt.WaveSlot(sig,
                               Position=ip,
                                Alpha=0.5
                               ))

Splt = Rplt.PlotSlots(Slots[0:10])
Splt.PlotChannels(None, Units='uV')
Splt.AddLegend(**LegendKwargs)


#%%

FluoPeaks = Rana.threshold_detection(Seg.GetSignal("Mean1"), 170, sign = "above", RelaxTime = 0)
CurrentPeaks = Rana.threshold_detection(Seg.GetSignal("16"),  -1.9e3*pq.uV, sign = "below", RelaxTime = 1*pq.s)

# CurrentPeaks = [CurrentPeaks[1], CurrentPeaks[4], CurrentPeaks[6],
#                 CurrentPeaks[9], CurrentPeaks[12], CurrentPeaks[-2], CurrentPeaks[-1]]

CurrentPeaks = CurrentPeaks[-7:-1]

#%%

d1 = (Pulses[-1] - Pulses[0])
# d1 = (CurrentPeaks[-2] - CurrentPeaks[0])
d1 = (53 - 7.25)*pq.s
d2 = FluoPeaks[-1] - FluoPeaks[0]
 
resamp = d2/d1
#%%

TimeFluo = np.linspace(Seg.GetSignal("Mean1").times[0], Seg.GetSignal("Mean1").times[-1], len(Seg.GetSignal("Mean1"))).rescale("s")

TimeMEA = np.linspace(Seg.GetSignal(InChans[0]).times[0], Seg.GetSignal(InChans[0]).times[-1], len(Seg.GetSignal(InChans[0]))).rescale("s")

plt.plot(TimeMEA, Seg.GetSignal(InChans[0])*1e5, label = "MEA")
plt.plot(TimeFluo, Seg.GetSignal("Mean1")-165*pq.uV, label = "Fluorescence")
plt.xlabel("Time (s)")
plt.ylabel("Voltage ($\\mu$V)")
plt.legend()


#%% FluoPeak Background substraction and normalization

NormSeg = NeoSegment()
BackgroundWindow = (t_start, 7*pq.s)

for sig in Seg.Signals()[1:-1]:
    # Max = max(sig.GetSignal(None))
    Background = np.mean(sig.GetSignal(BackgroundWindow))
    NormalizedSignal = neo.AnalogSignal((sig.GetSignal(None) - Background),
                            name= sig.name,
                            sampling_rate = sig.sampling_rate,
                            t_start=sig.t_start,
                            units=pq.uV,
                            annotate=sig.annotations["annotate"],
                            file_origin=sig.file_origin)
    NormSeg.AddSignal(NormalizedSignal)
    
    
NormSlots = []
    
for ip, sig in enumerate(NormSeg.Signals()):
    NormSlots.append(Rplt.WaveSlot(sig,
                               Position=ip,
                                Alpha=0.5
                               ))

NormSplt = Rplt.PlotSlots(NormSlots[0:10])
NormSplt.PlotChannels(None, Units='uV')
NormSplt.AddLegend(**LegendKwargs)
    
#%% Detect FluoPeaks

Pulses = Seg.GetEventTimes("EndPulse").rescale("s")
PulsesWindowLength = 1*pq.s
thr_rms = 2.75

FluoPeaks = {}

for sig in NormSeg.Signals()[0:-1]:
    rms = []
    for sl, RMSwindow in enumerate(Pulses):
        if sl == 0:
            continue
        rms.append(np.sqrt(np.mean(sig.GetSignal((Pulses[sl-1]+PulsesWindowLength, RMSwindow-PulsesWindowLength))**2)))
    FluoPeaks[sig.name] = (Rana.threshold_detection(sig, threshold = np.mean(rms)*thr_rms, sign = "above", RelaxTime = (1/resamp)*pq.s), np.mean(rms)*thr_rms)
    

for sl, value in enumerate(list(FluoPeaks.values())[0:10]):
    ax = NormSplt.Axs[sl]
    x_lim = ax.get_xlim()*pq.s
    ax.hlines(value[1], x_lim[0], x_lim[1], color = "r", Alpha = 0.5)
    y_lim = ax.get_ylim()*pq.uV
    for event in value[0]:
        ax.vlines(event, y_lim[0], y_lim[1], color = "r", Alpha = 0.5)
    for item in Pulses:
        ax.vlines(item.rescale("s"), y_lim[0], y_lim[1], color = "black", Alpha = 0.5, linestyle = "dashed")
        
#%% Find FluoPeaks inside stimulation window

PeaksInWindows = {}
PulsesWindow = []

for pulses in Pulses:
   PulsesWindow.append((pulses-PulsesWindowLength, pulses+PulsesWindowLength))

for roi, peaks in FluoPeaks.items():
    dummy = []       
    for peak in peaks[0]:
        for window in PulsesWindow:
            if peak > window[0] and peak < window[1]:
                dummy.append(peak)
    if len(dummy)>0:
        dummy.append(peaks[1])
        PeaksInWindows[roi] = dummy
        
            
#%% Plot FluoPeaks inside window

PeakSlot = []
    
for ip, roi in enumerate(list(PeaksInWindows.keys())[0:10]):
    PeakSlot.append(Rplt.WaveSlot(NormSeg.GetSignal(roi),
                               Position=ip,
                                Alpha=0.5
                               ))

PeakSplt = Rplt.PlotSlots(PeakSlot)
PeakSplt.PlotChannels(None, Units='uV')
PeakSplt.AddLegend(**LegendKwargs)


for sl, value in enumerate(list(PeaksInWindows.values())[0:10]):
    ax = PeakSplt.Axs[sl]
    x_lim = ax.get_xlim()*pq.s
    ax.hlines(value[1], x_lim[0], x_lim[1], color = "r", Alpha = 0.5)
    y_lim = ax.get_ylim()*pq.uV
    for item in Pulses:
        ax.vlines(item.rescale("s"), y_lim[0], y_lim[1], color = "black", Alpha = 0.5, linestyle = "dashed")
    if len(value) == 2:
        ax.vlines(value[0], y_lim[0], y_lim[1], color = "r", Alpha = 0.5)
    else:
        for event in value[0:-1]:
            ax.vlines(event, y_lim[0], y_lim[1], color = "r", Alpha = 0.5)


#%% Save plot and dict

plt.tight_layout()
plt.savefig(Path + FluoFile.split(".")[0] + "_FluoPeaks.svg")
plt.savefig(Path + FluoFile.split(".")[0] + "_FluoPeaks.png", dpi = 300)


PeaksInWindows["Pulses"] = Pulses
pickle.dump(PeaksInWindows, open( Path + FluoFile.split(".")[0] + "_FluoPeaks" +".p", "wb" ))
