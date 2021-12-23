# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:56:37 2020

@author: jdelacruz
"""

import glob
import OpenEphys
import matplotlib.pyplot as plt
import quantities as pq
from PhyREC.NeoInterface import NeoSegment
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
import McsPy.McsData as McsData



def LoadMCS(FileSearch):
    
    FileNames = glob.glob(FileSearch)

    for InFile in FileNames:
        
         
        OutSeg = NeoSegment()
        
        Dat = McsData.RawData(InFile)
        Rec = Dat.recordings[0]
        
        NSamps = Rec.duration
        
        for AnaStrn, AnaStr in Rec.analog_streams.items():
            if len(AnaStr.channel_infos) == 1:
                continue
        
            for Chn, Chinfo in AnaStr.channel_infos.items():
                print('Analog Stream ', Chinfo.label, Chinfo.sampling_frequency)
                ChName = str(Chinfo.label)
                print(ChName)
        
                Fs = Chinfo.sampling_frequency
                Var, Unit = AnaStr.get_channel_in_range(Chn, 0, NSamps)
                sig = neo.AnalogSignal(pq.Quantity(Var, Chinfo.info['Unit']),
                                       t_start=0*pq.s,
                                       sampling_rate=Fs.magnitude*pq.Hz,
                                       name=ChName)
                
                # DownSamp = RPro.DownSampling(sig, 25)
                # OutSeg.AddSignal(DownSamp)
                
                OutSeg.AddSignal(sig)

            
         
                
        
        # Trigger = Rec.event_streams[1].event_entity[0].get_event_timestamps()[0][0::2]*pq.us
        # TriggerEvent = neo.Event(times = Trigger, name = "Trigger", units=Trigger.units)
        # OutSeg.Seg.events.append(TriggerEvent)
        # OutSeg.UpdateEventDict()
        
        
    return OutSeg

def get_key(my_dict, val): 
    for key, value in my_dict.items(): 
         if val == value: 
             return key 
  
    return "key doesn't exist"
    
        

#%% Load Files

if __name__ == '__main__':


    Path = '20200630_ERG_LE_V2_C6-4\\R5.h5'

     
    Slots = []
    
    Seg = LoadMCS(Path)
    
    
    F1 = [{'function': RPro.Filter, 'args': {'Type':'bandpass',
                                            'Order':2,
                                            'Freqs':(.5, 150)}},
      {'function': RPro.RemoveDC, 'args': {}}
     ]    
    
#Create horizontal Waveslot

    for ip, sig in enumerate(Seg.Signals()):
        sig.ProcessChain = F1
        Slots.append(Rplt.WaveSlot(sig,
                                   Position=ip,
                                    Alpha=0.5
                                   ))
    
    Splt = Rplt.PlotSlots(Slots)
    Splt.PlotChannels(None, Units='uV')
    
    
#%%

AvgSeg = NeoSegment()

MCStoZiff = {
                "E1": "Ziff16", 
                 "E2": "Ziff1",
                 "E3": "Ziff15",
                 "E4": "Ziff2",
                 "E5": "Ziff14",
                 "E6": "Ziff3",
                 "E7": "Ziff13",
                 "E8": "Ziff4",
                 "E9": "Ziff12",
                 "E10": "Ziff5",
                 "E11": "Ziff11",
                 "E12": "Ziff6",
                 "E13": "Ziff10",
                 "E14": "Ziff7",
                 "E15": "Ziff9",
                 "E16": "Ziff8",
                 "E17": "Ziff24",
                 "E18": "Ziff25",
                 "E19": "Ziff23",
                 "E20": "Ziff26",
                 "E21": "Ziff22",
                 "E22": "Ziff27",
                 "E23": "Ziff21",
                 "E24": "Ziff28",
                 "E25": "Ziff20",
                 "E26": "Ziff29",
                 "E27": "Ziff19",
                 "E28": "Ziff30",
                 "E29": "Ziff18",
                 "E30": "Ziff31",
                 "E31": "Ziff17",
                 "E32": "Ziff32"
                 }

ZiffToElect = {
                "Ziff1": "E24", 
                 "Ziff2": "E25",
                 "Ziff3": "E26",
                 "Ziff4": "E27",
                 "Ziff5": "E28",
                 "Ziff6": "E29",
                 "Ziff7": "E30",
                 "Ziff8": "E31",
                 "Ziff10": "E17",
                 "Ziff11": "E18",
                 "Ziff12": "E19",
                 "Ziff13": "E20",
                 "Ziff14": "E21",
                 "Ziff15": "E22",
                 "Ziff16": "E23",
                 "Ziff17": "E9",
                 "Ziff18": "E10",
                 "Ziff19": "E11",
                 "Ziff20": "E12",
                 "Ziff21": "E13",
                 "Ziff22": "E14",
                 "Ziff23": "E15",
                 "Ziff24": "E16",
                 "Ziff25": "E1",
                 "Ziff26": "E2",
                 "Ziff27": "E3",
                 "Ziff28": "E4",
                 "Ziff29": "E5",
                 "Ziff30": "E6",
                 "Ziff31": "E7",
                 "Ziff32": "E8"
                 }

Twind = (-.05, .35)*pq.s
# Trig = LoadOpenEphysEvents(FilePath)
# TrigTimes = Trig["timestamps"][0::2]/float(Trig["header"]["sampleRate"])*pq.s
# CorrectedTrigTimes = TrigTimes - 158.1064*pq.s

Trigger = []
k = 0

for i in range(20):
    Trigger.append(k)
    k = k+2

TrigTimes = Trigger*pq.s + 13.415*pq.s


for sig in Slots:
    if sig.name == "E15": #Exclude not connected ziff pad
        continue
    else:        
        AvgSig = neo.AnalogSignal(sig.CalcAvarage(TimesEvent=TrigTimes,
                                                  TimeAvg=Twind),
                                    t_start=Twind[0],
                                    sampling_rate=sig.GetSignal(None).sampling_rate,
                                    # name=sig.name
                                    name=ZiffToElect[MCStoZiff[sig.name]]
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



#%%

BWaveDict = {}
BWaveMean = [] 


for sig in AvgSeg.signames:
    Max = max(AvgSeg.GetSignal(sig))
    BWaveDict[sig] = Max
    BWaveMean.append(Max) 

    
BWaveMean = np.mean(BWaveMean)


ElectrodeGrid = {"E1": (0, 2600), #The order of the Eannels is flipped because the Ziff adapter was solderer upside down
                 "E2": (-1500, 2600),
                 "E3": (-1000, 1700),
                 "E4": (-2300, 1300),
                 "E5": (-400, 800),
                 "E6": (-1500, 800),
                 "E7": (-3000, 0),
                 "E8": (-2000, 0),
                  "E9": (-800, 0),
                 "E10": (-2300, -1300),
                 "E11": (-1500, -800),
                 "E12": (-400, -800),
                 "E13": (-1500, -2600),
                 "E14": (-1000, -1700),
                 "E15": (0, -1700),
                 "E16": (0, -2600),
                 "E17": (0, 0),
                 "E18": (1000, -1700),
                 "E19": (1500, -2600),
                 "E20": (400, -800),
                 "E21": (1500, -800),
                 "E22": (2300, -1300),
                 "E23": (800, 0),
                 "E24": (2000, 0),
                 "E25": (3000, 0),
                 "E26": (1500, 800),
                 "E27": (400, 800),
                 "E28": (2300, 1300),
                 "E29": (1000, 1700),
                 "E30": (1500, 2600),
                 "E31": (0, 1700)
                 }


X, Y, bWave3D = [], [], []

for Channel, Coordinate in ElectrodeGrid.items():
    X.append(Coordinate[0])
    Y.append(Coordinate[1])
    bWave3D.append(BWaveDict[Channel])
                


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, bWave3D, zdir='z', s=20, c=None, depthshade=True)

ax.set_xlabel("X Axis (um)")
ax.set_ylabel("Y Axis (um)")
ax.set_zlabel("bWave (uV)")

# Plot injury quadrant

xLine = np.linspace(0, -3000, 50)
yLine = np.linspace(0, 2600, 50)
lol = np.linspace(0, 0, 50)

ax.plot(xLine, lol, min(bWave3D), color = "r")
ax.plot(lol, yLine, min(bWave3D), color = "r")
ax.plot(xLine, 2600+ lol, min(bWave3D), color = "r")
ax.plot(-3000+lol, yLine, min(bWave3D), color = "r")
             
                 
#%% Create grid for 3D plot
                 
GridNodes = [-3000, -2600, -2300, -2000, -1700, -1500, -1300, -1000, -800, -400, 0, 400, 800, 1000, 1300, 1500, 1700, 2000, 2300, 2600, 3000]


XGrid = np.zeros((len(GridNodes), len(GridNodes)))
YGrid = np.zeros((len(GridNodes), len(GridNodes)))

for sl, i in enumerate(GridNodes):
    XGrid[:, sl] = i
    YGrid[sl, :] = i



#%% Y interpolation


yForYInter = [
              [3000, 0, -3000],
              [3000, 1300, -1300, -3000],
              [3000, 0, -3000],
              [3000, 2600, 800, -800, -2600, 3000],
              [3000, 1700, -1700, -3000],
              [3000, 0, -3000],
              [3000, 800, -800, -3000],
              [3000, 2600, 1700, 0, -1700, -2600, -3000],
              [3000, 800, -800, -3000],
              [3000, 0, -3000],
              [3000, 1700, -1700, -3000],
              [3000, 2600, 800, -800, -2600, 3000],
              [3000, 0, -3000],
              [3000, 1300, -1300, -3000],
              [3000, 0, -3000]      
              ]


XBlueprint = [-3000, -2300, -2000, -1500, -1000, -800, -400, 0, 400, 800, 1000, 1500, 2000, 2300, 3000]


bWaveForYInter = []

for sl, x in enumerate(XBlueprint):
    dummy = [BWaveMean]
    
    for y in yForYInter[sl]:
        if y == -3000 or y == 3000:
            continue
        key = get_key(ElectrodeGrid,(x,y))
        if key == "key doesn't exist":
            continue
        else:
            dummy.append(BWaveDict[key])
            # print(key , BWaveDict[key])
    
    dummy.append(BWaveMean)
    bWaveForYInter.append(dummy)
    
InterpolationObjects = []


for sl, i in enumerate(yForYInter):
    InterpolationObjects.append(interpolate.interp1d(i, bWaveForYInter[sl], fill_value = "extrapolate"))
    
Z = np.zeros((21,21))

i = 0
for sl, column in enumerate(YGrid.T):
    if sl in (1, 4, 6, 14, 16, 19):
        Z[:, sl] = BWaveMean
    else:
        Z[:, sl] = InterpolationObjects[i](column)
        i = i+1
        
#%% 1D Inter plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


# Plot a basic wireframe.
ax.plot_wireframe(XGrid, YGrid, Z, rstride=1, cstride=1)
ax.scatter(X, Y, bWave3D, zdir='z', s=20, color = "r", depthshade=True)


    

ax.set_xlabel("X Axis (um)")
ax.set_ylabel("Y Axis (um)")
ax.set_zlabel("bWave (uV)")


# Plot injury quadrant

xLine = np.linspace(0, -3000, 50)
yLine = np.linspace(0, 2600, 50)
lol = np.linspace(0, 0, 50)

ax.plot(xLine, lol, min(bWave3D), color = "r")
ax.plot(lol, yLine, min(bWave3D), color = "r")
ax.plot(xLine, 2600+ lol, min(bWave3D), color = "r")
ax.plot(-3000+lol, yLine, min(bWave3D), color = "r")

#%% Export data

file = open(Path.split(".")[0] + "_Z", "w")
np.savetxt(file, Z.T)
file.close()

file = open(Path.split(".")[0] + "_bWave3D", "w")
np.savetxt(file, bWave3D)
file.close()

# #%% 2D Interpolation

# X2D, Y2D, Z2D = [], [], []

# for k, v in ElectrodeGrid.items():
#     if k == "CH9":
#         X2D.append(v[0])
#         Y2D.append(v[1])
#         Z2D.append(BWaveMean)
        
#     else:        
#         X2D.append(v[0])
#         Y2D.append(v[1])
#         Z2D.append(BWaveDict[k].magnitude.tolist()[0])
#         # print(k)
    
    
# # for item in GridNodes:
# #     if item == 0:
# #         continue
# #     X2D.append(-3000)
# #     Y2D.append(item)
# #     Z2D.append(BWaveMean)
    
# #     X2D.append(3000)
# #     Y2D.append(item)
# #     Z2D.append(BWaveMean)
    
# #     X2D.append(item)
# #     Y2D.append(-3000)
# #     Z2D.append(BWaveMean)
    
# #     X2D.append(item)
# #     Y2D.append(3000)
# #     Z2D.append(BWaveMean)

# # #Manual points to improve stability 

# # X2D.append(2300), Y2D.append(400), Z2D.append(BWaveMean)
# # X2D.append(1500), Y2D.append(-400), Z2D.append(BWaveMean)



# f = interpolate.interp2d(X2D, Y2D, Z2D) #, kind='quintic')


# Z2Dnew = f(GridNodes, GridNodes)
    
    
# #%% 2D Inter plot

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')


# # Plot a basic wireframe.
# ax.plot_wireframe(XGrid, YGrid, Z2Dnew, rstride=1, cstride=1)
# ax.scatter(X, Y, bWave3D, zdir='z', s=20, color = "r", depthshade=True)

# # for k, v in ElectrodeGrid.items():
# #     ax.scatter(v[0], v[1], BWaveMean, color = "r")    
    

# ax.set_xlabel("X Axis (um)")
# ax.set_ylabel("Y Axis (um)")
# ax.set_zlabel("bWave (uV)")
    
    
#%% Calculate differences vs average

SigToAvg = []
ExcludeFromAvg = [
                    "E1"
                    ]


for sig in AvgSeg.signames:
    if sig in ExcludeFromAvg:
        continue
    else:
        SigToAvg.append(AvgSeg.GetSignal(sig))
                 
AvgWave = np.mean(SigToAvg, axis = 0)*pq.uV

SubtractSeg = NeoSegment()

for sig in AvgSeg.Signals():
     SubtractSeg.AddSignal(sig.duplicate_with_new_data(signal = sig - AvgWave))

GridDict = {
                "E1": 4, 
                 "E2": 2,
                 "E3": 11,
                 "E4": 10,
                 "E5": 21,
                 "E6": 20,
                 "E7": 27,
                 "E8": 28,
                  "E9": 30,
                 "E10": 37,
                 "E11": 38,
                 "E12": 39,
                 "E13": 56,
                 "E14": 47,
                 "E15": 49,
                 "E16": 58,
                 "E17": 31,
                 "E18": 51,
                 "E19": 60,
                 "E20": 41,
                 "E21": 42,
                 "E22": 43,
                 "E23": 32,
                 "E24": 34,
                 "E25": 35,
                 "E26": 24,
                 "E27": 23,
                 "E28": 16,
                 "E29": 15,
                 "E30": 6,
                 "E31": 13
                 }

#%% Plot differences vs average

# fig, axs = plt.subplots(4, 8, sharex=True, sharey=True)
fig, axs = plt.subplots(7, 9, sharex=True, sharey=True)
axs = axs.ravel()


# for sl, ch in enumerate(AvgSeg.signames):
#     axs[sl].plot(AvgSeg.GetSignal(ch).times, AvgSeg.GetSignal(ch))
#     axs[sl].plot(SubtractSeg.GetSignal(ch).times, SubtractSeg.GetSignal(ch))
#     axs[sl].set_title(sl)


for ch in AvgSeg.signames:
    axs[GridDict[ch]].plot(AvgSeg.GetSignal(ch).times, AvgSeg.GetSignal(ch), linewidth=.6)
    # axs[GridDict[ch]].plot(SubtractSeg.GetSignal(ch).times, SubtractSeg.GetSignal(ch))
    axs[GridDict[ch]].set_title(ch)
    
    
    
for k in range(63):
    if k not in GridDict.values():
        axs[k].axis("off")
        
        
#%% Single trial average

plt.plot(Slots[1].GetSignal(TrigTimes[0]+Twind))

SingleTrialAverage = []

for signal in Slots:
    if signal.name == "E18":
        continue
    else:        
        SingleTrialAverage.append(signal.GetSignal(TrigTimes[0]+Twind))
        

plt.plot(np.array(SingleTrialAverage).mean(0))
        
    