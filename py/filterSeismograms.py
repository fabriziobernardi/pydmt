#########################################################################
#                                                                       #
# Set of python functions to remove instrument and filtering            #
#                                                                       #
#########################################################################

import os,sys
import obspy.signal
from obspy.sac import attach_paz
from obspy.signal import seisSim,cosTaper
from scipy.signal import detrend, resample
from processData import purgeListStation

def removeInstrument(st,args):

    if(args.sim == 'PZs'):

       # prefilters
       f = args.flim.split()
       f0 = eval(f[0])
       f1 = eval(f[1])
       f2 = eval(f[2])
       f3 = eval(f[3])
       toPurge=  []   # station to purge if no Paz found

       for i in range(len(st)):
   
           # attach poles and zeros instrument
           if(args.dva=='1'):
              try:
                attach_paz(st[i], st[i].stats.PZs_file,todisp=False)
              except:
                print "No appropriate PZs file found for station " + st[i].stats.station,st[i].stats.channel,st[i].stats.network
                toPurge.append(st[i].stats.station)
           else:
              try:
                attach_paz(st[i], st[i].stats.PZs_file,tovel=True)
              except:
                print "No appropriate PZs file found for station " + st[i].stats.station,st[i].stats.channel,st[i].stats.network
                toPurge.append(st[i].stats.station)

                
       # remove stations if len(toPurge>0)
       if len(toPurge) > 0:
           st = purgeListStation(st,toPurge,'r')
           print "Check if station/channel/network/location of the PZs files and the same string within loaded binary files "
           print "do correspond. It may occour for instance that the headers strings of the waveform files (e.g. sac, fseed) "
           print "do not agrees with the same strings of the PZs name files. For instance the name of the network. "
           print "If these strings do not correspond, modify the name of the PZs files or the header values of the waveforms"
           print "You may also choose to remove this station using the option --purge (see help for details)"

       # now do remove
       for i in range(len(st)):

           # remove instrument to displacement
#          st[i].data=detrend(st[i].data)
           st[i].data = seisSim(st[i].data,st[i].stats.sampling_rate,paz_remove=st[i].stats.paz, \
                        taper=True, taper_fraction=0.050, pre_filt=(f0,f1,f2,f3)) #,water_level=60.0) 

           # from meters to centimeters
           st[i].data = st[i].data * 100
       
    
    return st


def filtering(st,ty,args):

    hip = args.highpass 
    lop = args.lowpass 
    bdp = args.bandpass 

    if hip != "0":
       elements = hip.split()
       cors = int(elements[0])
       freq = eval(elements[1])
       for i in range(len(st)):
           st[i].data =  obspy.signal.highpass(st[i].data, freq, \
           df=st[i].stats.sampling_rate, corners=cors, zerophase=args.zeroph)

    if lop != "0":
       elements = lop.split()
       cors = int(elements[0])
       freq = eval(elements[1])
       for i in range(len(st)):
           st[i].data =  obspy.signal.lowpass(st[i].data, freq, \
           df=st[i].stats.sampling_rate, corners=cors, zerophase=args.zeroph)

    if bdp != "0":
       elements = bdp.split()
       cors = int(elements[0])
       freq_min = eval(elements[1])
       freq_max = eval(elements[2])
       for i in range(len(st)):
           st[i].data =  obspy.signal.bandpass(st[i].data, freq_min, freq_max, \
           df=st[i].stats.sampling_rate, corners=cors, zerophase=args.zeroph)

    return st



def FilterData(tr,bdp,hip,lop):

    if hip != "0":
       elements = hip.split(' ')
       cors = int(elements[0])
       freq = eval(elements[1])
       for i in range(len(tr)):
          tr[i].filter("highpass", freq=freq, corners=cors, zerophase="False")
    elif lop != "0":
       elements = lop.split(' ')
       cors = int(elements[0])
       freq = eval(elements[1])
       #print lop,elements,cors,freq
       for i in range(len(tr)):
          tr[i].filter("lowpass", freq=freq, corners=cors, zerophase="False")
    elif bdp != "0":
       elements = bdp.split(' ')
       cors = int(elements[0])
       frem = eval(elements[1])
       freM = eval(elements[2])
       for i in range(len(tr)):
          tr[i].filter("bandpass", freqmin=frem, freqmax=freM, corners=cors, zerophase="False")
    else:
       pass

    return tr

def decimate(st,spin,args):

    Df = eval(args.DeltaInv)
    # if data is green (g), then find here decimation factor
    if (spin == 'g'):
       frs = 1/st[0].stats.delta
       frd = 1/eval(args.delta)
       fri = 1/eval(args.DeltaInv)
       c  = int(1/eval(args.delta))/(1/eval(args.DeltaInv))
       ny = 1.0/(2*eval(args.delta))
       c = int(c)
       cn = (frs/fri)
       for i in range(len(st)):
           st[i].filter("lowpass", freq=ny*0.5, corners=1, zerophase="False")
           st[i].decimate(c,strict_length=False, no_filter=True)

    elif (spin== 'd'):
       for i in range(len(st)):
           c  = int(1/st[i].stats.delta)/(1/eval(args.DeltaInv))
           ny = 1.0/(2*st[i].stats.delta)
           c = int(c)
           st[i].filter("lowpass", freq=ny*0.5, corners=1, zerophase="False")
           st[i].decimate(c,strict_length=False, no_filter=True)

    else:
       print "Decimation option not recognized. EXIT!"
       sys.exit()

     # check if correct sampling rate, else interpolate
    if(args.inter == 'Y' and spin== 'd'):
      for i in range(len(st)):
          if (st[i].stats.npts != int(args.len)):
             if args.war == 'Y':
                print "Station " + st[i].stats.station + "." + st[i].stats.channel + " resampled from " + st[i].stats.delta + " to" + args.DeltaInv
             st[i].data = resample(st[i].data,int(args.len))
             st[i].stats.delta = float(args.DeltaInv) 
           

    return st


def decimateStream(self, factor):

    c=int(factor)
    for i in range(len(self)):
       self[i].downsample(c,strict_length=False, no_filter=True)

    return self

