# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:56:37 2020

@author: jdelacruz
"""

import matplotlib.pyplot as plt
import quantities as pq
from PhyREC.NeoInterface import NeoSegment
import neo
import PhyREC.SignalProcess as RPro
import numpy as np
import glob
import PhyREC.PlotWaves as Rplt
import PhyREC.SignalAnalysis as Rana
from sklearn.cluster import MiniBatchKMeans, KMeans, SpectralClustering
import pickle
import OpenEphys
from matplotlib import cm
import matplotlib.colors as mpcolors



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
        print(ChName)
        
        
    return MySeg

def LoadOpenEphysEvents(FilePath):
    FileName = FilePath.split("*")[0] + "all_channels.events"
    Trigger = OpenEphys.loadEvents(FileName)
    
    return Trigger

        


def PlotEvents(Signal, Events, TimeAvg, Alpha = .2):
    
    Fig,Ax = plt.subplots()
    
    Ts = Signal.sampling_period
    nSamps = int((TimeAvg[1]-TimeAvg[0])/Ts)
    time = np.arange(nSamps)*Ts + TimeAvg[0]
    Fs = Signal.sampling_rate
    
    OutSeg = NeoSegment()
    
    for event in Events:
       
        OutSig = neo.AnalogSignal(Signal.GetSignal((event + TimeAvg[0], event + TimeAvg[1])),
                                            t_start=TimeAvg[0],
                                            sampling_rate=Fs,
                                            name=event.magnitude)
        OutSeg.AddSignal(OutSig)
       
        Ax.plot(time.rescale("ms"), Signal.GetSignal((event + TimeAvg[0], event + TimeAvg[1])), label = Signal.name, 
                color = "black", alpha = Alpha)
        Ax.set_xlabel("Time (ms)")
        Ax.set_ylabel("Amplitude ($\mu$V)")
        Ax.set_title(Signal.name)
        
    return OutSeg

def PlotAverage(Data, Time, color, alpha = .5):
    
    FigAv, AxAv = plt.subplots()
    
    Mean = np.mean(Data, axis = 0).flatten()
    Std = np.std(Data, axis = 0).flatten()
    
    AxAv.plot(Time, Mean, color = color, lw = 2)
    AxAv.fill_between(Time, Mean + Std, Mean - Std, color = color, alpha = alpha)
    
    AxAv.set_xlabel("Time (ms)")
    AxAv.set_ylabel("Voltage ($\mu$V)")
    AxAv.tick_params(direction="in")
    # AxAv.set_ylim(-55, 45)

    return 


    
    
        
#%%


FilePath = 'MEA1\\1to32\\2020-03-06_21-17-02\\*.continuous'

InChans = (
            # 'CH1',
              # "CH2",
             # 'CH3',
             'CH4',
             # 'CH5',
              'CH6', 
             'CH7',
             # "CH8",
             # 'CH9',
             # 'CH10',
               # 'CH11',
              'CH12', 
               # 'CH13',
              # "CH14",
              # 'CH15',
              # 'CH16',
              # 'CH17',
              "CH18",
              # 'CH19',
              # 'CH20',
              # 'CH21',
              # 'CH22', 
              # 'CH23',
              "CH24",
              'CH25',
              'CH26',
              'CH27',
              # 'CH28', 
              "CH29",
              # 'CH30',
               "CH31",
              # 'CH32'
           )

# InChans = (
#             # 'CH1',
#               # "CH2",
#              # 'CH3',
#              'CH4',
#              # 'CH5',
#               'CH6', 
#              'CH7',
#              # "CH8",
#              # 'CH9',
#              # 'CH10',
#                # 'CH11',
#               'CH12', 
#                # 'CH13',
#               # "CH14",
#               # 'CH15',
#               # 'CH16',
#               # 'CH17',
#               "CH18",
#               # 'CH19',
#               'CH20',
#               # 'CH21',
#               # 'CH22', 
#               'CH23',
#               "CH24",
#               'CH25',
#               'CH26',
#               'CH27',
#               # 'CH28', 
#               "CH29",
#               # 'CH30',
#                "CH31",
#               'CH32'
#            )

# InChans = []

# for i in range(64):
#     InChans.append("CH{}".format(i+1))

Seg = LoadOpenEphysContinuous(FilePath, IncludeCh=InChans)




cmap = cm.ScalarMappable(mpcolors.Normalize(vmin=0, vmax=len(Seg.SigNames.keys())), cm.jet)
colors = [cmap.to_rgba(i) for i in range(len(Seg.SigNames.keys()))]


#TWind = (10*pq.s, 20*pq.s)
Fband = (48, 52)
Fband2 = (98, 102)
Fband3 = (148, 152)
Fband4 = (10, 1000)


F1 = [{'function': RPro.Filter, 'args': {'Type':'bandstop',
                                                'Order':2,
                                                'Freqs':(Fband[0], Fband[1])}},
      {'function': RPro.Filter, 'args': {'Type':'bandstop',
                                                'Order':2,
                                                'Freqs':(Fband2[0], Fband2[1])}},
      {'function': RPro.Filter, 'args': {'Type':'bandstop',
                                                'Order':2,
                                                'Freqs':(Fband3[0], Fband3[1])}},
       {'function': RPro.Filter, 'args': {'Type':'bandpass',
                                                 'Order':2,
                                                 'Freqs':(Fband4[0], Fband4[1])}},
     {'function': RPro.RemoveDC, 'args': {}}
     ]


Slots = []

for ip, sig in enumerate(Seg.Signals()):
    sig.ProcessChain = F1
    Slots.append(Rplt.WaveSlot(sig,
                               Units='uV',
                               Color=colors[ip],
                               Alpha=1,
                               ))
        
LegendKwargs = {'fontsize': "large",
                'ncol': 5,
                'loc': 'upper right',
                'frameon': False}

TWind = (0.2, 15)*pq.s
plt.rcParams.update({'font.size': 18})

Splt = Rplt.PlotSlots(Slots)
Splt.PlotChannels(TWind, Units='uV')
Splt.AddLegend(**LegendKwargs)
plt.tight_layout()
plt.subplots_adjust(hspace=0)    


#%% Threshold detection

# sig = Slots[0].GetSignal((0, 200)*pq.s)

# Threshold = -5*np.std(sig)

ShortSlots = []

for ip, sig in enumerate(Seg.Signals()):
    sig.ProcessChain = F1
    ShortSlots.append(Rplt.WaveSlot(sig.GetSignal(TWind),
                               Units='uV',
                               Color=colors[ip],
                               Alpha=1,
                               ))

Threshold = {
                "CH18": -20*pq.uV,
                "CH24": -15*pq.uV,
                "CH25": -20*pq.uV,
                "CH26": -20*pq.uV,
                "CH27": -20*pq.uV,
                "CH29": -20*pq.uV,
                "CH31": -18*pq.uV,
                "CH4": -25*pq.uV,
                "CH6": -20*pq.uV,
                "CH7": -15*pq.uV
                }


RefractoryPeriod = 2*pq.ms

Spikes = {}

for sl, slot in enumerate(ShortSlots):
    sig = slot.GetSignal(None)
    Spikes[sig.name] = Rana.threshold_detection(sig, Threshold[sig.name], sign = "below", RelaxTime = RefractoryPeriod)
 



    ax = Splt.Axs[sl]
    
    x_lim = ax.get_xlim()*pq.s
    y_lim = ax.get_ylim()*pq.uV
    ax.hlines(Threshold[sig.name], x_lim[0], x_lim[1], color = "r", Alpha = 0.5)
        
        
    for event in Spikes[sig.name]:
        ax.vlines(event, y_lim[0], y_lim[1], color = "r", Alpha = 0.5)




#%% Plot Event Average

TimeAvg = (-1, 2)*pq.ms
TimeAvg = TimeAvg.rescale("s")

SegOut = {}

for sl, slot in enumerate(Slots):
    sig = slot.GetSignal(None)
    SegOut[sig.name] = PlotEvents(sig, Spikes[sig.name], TimeAvg, Alpha = .1)
    
    
FinalSeg = NeoSegment()

for channel in SegOut.values():
    for signal in channel.Signals():
        FinalSeg.AddSignal(signal)
    


#%% Calculate height and width


Height_Width = {}

H_W_list = []

for sig in FinalSeg.signames:
    
    a=FinalSeg.GetSignal(sig)
    
    Max = max(a)
    Min = min(a)
    Height = Max-Min
    
    MaxTime = np.where(a == Max)[0]
    MinTime = np.where(a == Min)[0]
    Width = MinTime-MaxTime

    
    Height_Width[sig] = [Height, Width]
    H_W_list.append((Height, Width[0]))  
    
FigHW, AxHW = plt.subplots()    

for item in H_W_list:
    AxHW.scatter(item[0], item[1])
    

# k_means = KMeans(init='k-means++', n_clusters=4)
k_means = MiniBatchKMeans(init='k-means++', n_clusters=2)
label = k_means.fit_predict(H_W_list)
k_means_cluster_centers = k_means.cluster_centers_

# pickle.dump(label, open( Path + "Label.p", "wb" ))


FigHW_KMeans, AxHW_KM = plt.subplots()
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', "#d62728", "#9467bd", "#8c564b"]
 
for item, sl in zip(H_W_list, label):
    AxHW_KM.scatter(item[0], item[1], color = colors[sl])
    
for sl, item in enumerate(k_means_cluster_centers):
    AxHW_KM.scatter(item[0], item[1], color = colors[sl], marker = '*', edgecolor = "black")
    
AxHW_KM.set_xlabel("Amplitude ($\mu$V)")
AxHW_KM.set_ylabel("Duration (Samples)")
AxHW_KM.tick_params(direction="in")


#Assign cluster to spike

Sorted_HW_spikes = {}

for sl, (spike, params) in enumerate(Height_Width.items()):
    Sorted_HW_spikes[spike] = label[sl]

Time = FinalSeg.GetSignal(spike).times.rescale("ms")

Fig_HW_Colored, Ax_HW_Colored = plt.subplots()

for sig in FinalSeg.signames:
    Ax_HW_Colored.plot(Time, FinalSeg.GetSignal(sig), alpha = .1, color = colors[Sorted_HW_spikes[sig]])
    
Ax_HW_Colored.set_xlabel("Time (ms)")
Ax_HW_Colored.set_ylabel("Voltage ($\mu$V)")
Ax_HW_Colored.tick_params(direction="in")

# FigHW_KMeans.savefig(SavePath + "HW_Cluster.svg")
# FigHW_KMeans.savefig(SavePath + "HW_Cluster.png", dpi = 300)

# Fig_HW_Colored.savefig(SavePath + "HW_Spikes.svg")
# Fig_HW_Colored.savefig(SavePath + "HW_Spikes.png", dpi = 300)

#%% Separate and average spikes

labels = np.unique(label)

Label_0 = []
Label_1 = []
# Label_2 = []
# Label_3 = []
# Label_4 = []
# Label_5 = []

for spike, lb in Sorted_HW_spikes.items():
    if lb == 0:
        Label_0.append(FinalSeg.GetSignal(spike))
    elif lb == 1:
        Label_1.append(FinalSeg.GetSignal(spike))
    # elif lb == 2:
    #     Label_2.append(FinalSeg.GetSignal(spike))
    # elif lb == 3:
    #     Label_3.append(FinalSeg.GetSignal(spike))
    # elif lb == 4: 
    #     Label_4.append(SegOut.GetSignal(spike))
    # elif lb == 5:
    #     Label_5.append(SegOut.GetSignal(spike))
    else:
        print("Watch out, madafacka")
        

# plt.close("all")

PlotAverage(Label_0, Time, color = colors[0])
PlotAverage(Label_1, Time, color = colors[1])
# PlotAverage(Label_2, Time, color = colors[2])
# PlotAverage(Label_3, Time, color = colors[3])
# PlotAverage(Label_4, Time, color = colors[4])
# PlotAverage(Label_5, Time, color = colors[5])


#%%
Dummy = []
for sig in FinalSeg.Signals():
    Dummy.append(sig)
    
PlotAverage(Dummy, Time, "black")
