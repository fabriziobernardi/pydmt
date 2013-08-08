#########################################################################
#                                                                       #
# Set of python function to load observed data into a stream            #
#                                                                       #
#########################################################################

from dlaz import dlaz
from obspy.core import Stream
from obspy.core import UTCDateTime
from scipy.signal import detrend
from collections import Counter
from math import fabs
from obspy.signal import rotate_NE_RT
import sys
import numpy as np

def getSTationslist(self):

    list=[]
    n=0

    for i in range(len(self)):
       list.append(self[i].stats['station'])

    # to be shure are really ordered
    list.sort()

    # list of distances
    out  = []
    for i in range(len(self)-1):
       if(self[i].stats['station'] != self[i+1].stats['station']):
          out.append(self[i].stats['station'])
    out.append(self[-1].stats['station'])

    return (out)

def selectUniqueTraces(tr,args):

    # Test on orizontal component, since if only vertical component exists, 
    # no xcorr on horiz allowed --> crash
  
    st = Stream()
    List = []
    ST = []
    CleanList = []

    for i in range(len(tr)):
        if(tr[i].stats.channel[2:3] == "N"):
           List.append(tr[i].stats.station)

    for i in range(len(tr)):
       a = List.count(tr[i].stats.station)
       if(a > 1):
           ST.append(tr[i].stats.station)

    d = Counter(ST)
    for key in d:
        CleanList.append(key) 
    
    for i in range(len(tr)):
        if CleanList.count(tr[i].stats.station) == 0:
           st.append(tr[i]) 

    return st 

def rotateToGCP(tr):


    #begin loop over data stream
    for i in range(len(tr)-1):
      # split channel
      li0 = list(tr[i+0].stats['channel'])
      li1 = list(tr[i+1].stats['channel'])

      # chech if station and part 1 of channel is identical and location
      if li0[0] == li1[0] and li0[1] == li1[1] \
         and tr[i+0].stats['station']  == tr[i+1].stats['station']\
         and tr[i+0].stats['location'] == tr[i+1].stats['location']:

         rch = li0[0] + li0[1] + 'R'
         tch = li0[0] + li0[1] + 'T'

         # if yes 3 possibility: EN, NE , pass
         if li0[2]=="E" and li1[2]=="N":
            #baz
            baz = tr[i].stats['baz']
            if tr[i+0].stats['npts'] == tr[i+1].stats['npts']:
               # rotate 0-1
               (tr[i+1].data,tr[i+0].data) = rotate_NE_RT(tr[i+1].data,tr[i+0].data,baz)
               tr[i+0].stats['channel']=tch
               tr[i+1].stats['channel']=rch
               i=i+1
            else:
               print "Can't rotate ",tr[i+0].stats['station'],tr[i+0].stats['channel'], " and ", \
                      tr[i+1].stats['station'],tr[i+1].stats['channel']

         elif li0[2]=="N" and li1[2]=="E":
            #baz
            baz = tr[i].stats['baz']
            if tr[i+0].stats['npts'] == tr[i+1].stats['npts']:
#              # rotate 1-0
               (tr[i+0].data,tr[i+1].data) = rotate_NE_RT(tr[i+0].data,tr[i+1].baz)
               tr[i+1].stats['channel']=tch
               tr[i+0].stats['channel']=rch
               i=i+1
            else:
               print "Can't rotate ",tr[i+0].stats['station'],tr[i+0].stats['channel'], " and ", \
                      tr[i+1].stats['station'],tr[i+1].stats['channel']

         else:
            pass

    return tr


def removeMeanTrend(tr):

    for i in range(len(tr)):
       #remove mean
       tr[i].data=tr[i].data - tr[i].data.mean()
       #remove trend
       tr[i].data=detrend(tr[i].data)

    return tr


def cutWindow(st,ty,args):

    if(ty == 'd'):

       #set begin and end time for data extraction
       # define t1 and t2
       t1 = UTCDateTime(args.ori)
       t1 = t1-int(args.pre)
       t2 = t1+int(args.len)
       t1 = str(t1)[:19]
       t2 = str(t2)[:19] #this must be a string for
       Tb = UTCDateTime(t1)
       Te = UTCDateTime(t2)
       #trim data
       st.trim(starttime=Tb, endtime=Te)

    if(ty == 's'):

      #here add Zcor
      for i in range(len(st)):
         #set begin and end time for data extraction
         # define t1 and t2
         #time2shift=st[i].stats.Zcor*st[i].stats.delta
         time2shift=0.

         
         t1 = UTCDateTime(args.ori)
         t1 = t1-int(args.pre)
         BB = t1+int(time2shift)
         bb = UTCDateTime(BB)
         st[i].stats.starttime=bb
         t2 = t1+int(args.len)
         t1 = str(t1)[:19]
         t2 = str(t2)[:19] #this must be a string for
         Tb = UTCDateTime(t1)
         Te = UTCDateTime(t2)

         st[i].trim(starttime=Tb, endtime=Te)

    return (st,Tb,Te)

def cleanStream(st,Tb,Te,args):

    # remove traces with gaps
    if (args.verbose==1):
       st = removeGaps(st, 0, 0, verbose="true")
    else:
       st = removeGaps(st, 0, 0, verbose="false")

    # remove shorter traces
    if len(st) != 0:
       st = removeShortTraces(st, 100,Tb,Te)
    else:
       print "No data left in the pool while cleaning the data-set"
       sys.exit()

    # remove stations with only 2 or less components
    if len(st) != 0:
       st = remove21Comp(st)
    else:
       print "No data left in the pool while cleaning the data-set"
       sys.exit()

    return st


def dlazStream(st):


    for i in range(len(st)):
 
        # value not set as initialized
        if(st[i].stats.evla != -1000.0):

           (distD,Az,Baz,distkm)=dlaz(st[i].stats.evla,st[i].stats.evlo,st[i].stats.stla,st[i].stats.stlo)
           st[i].stats.dist=distkm
           st[i].stats.gcarc=distD
           st[i].stats.baz=Baz
           st[i].stats.az=Az
       
        else:
           print "No station coordinates found for station ",st[i].stats.station,"! Check station file. Exit!!"

    return st
   

def selectData(st,args):


    # Only if a station set of data is required.
    # e.g.: --set "AQU MUGIO VLC TRI"
    if(args.set != 'All'):
      st = purgeListStation(st,args,'s')


    # purge stations (list)
    # remove specific stations from data set
    # e.g.: --purge "AGU TRI"
    if(args.purge != 'None'):
       st = purgeListStation(st,args,'p')

    # purge stations (distance)
    # accept only stations lying into a pecific range of distance in km
    # e.g.: --dist "300 400"
    if(args.range != 'None'):
      st = purgeListStation(st,args,'d')

    # purge stations (azimuth)
    # accept only stations lying into a pecific range of distance in km
    # e.g.: --dist "300 400"
    if(args.azi != 'None'):
      st = purgeListStation(st,args,'a')

    return st

def cleanZeroMaxTraces(st, args):

    #       1. One or more with max(data)= 0 
    toPurge = []
    for j in range(len(st)):
        if(max(st[j].data) == 0.0):

            if(len(toPurge) == 0):
               toPurge.append(st[j].stats.station)
            else:
               if (st[j].stats.station != toPurge[-1]):
                   toPurge.append(st[j].stats.station)

    st = purgeListStation(st,toPurge,'r')
    
    return st

def cleanMultipleFragmentTraces(st, args):

    #remove 1. Multiple fragment traces
    #       2. less than 3 components
    #       3. components different than Z,E,N

    #Instert into out stations to pugre
    toPurge = []

    # get station list
    a = getStationslist(st)

    # count nr components and select with comp >= 4
    for i in range(len(a)):

        count = 0
        for j in range(len(st)):

            if (st[j].stats.station == a[i]):
               count += 1
            
        if (count != 3):
           toPurge.append(a[i])

    #remove stations from st
    st = purgeListStation(st,toPurge,'r')

    # now remove stations with components differen from Z,N,E
    toPurge = []
    for j in range(len(st)):
        if (st[j].stats.channel[2:3] != 'Z' and \
            st[j].stats.channel[2:3] != 'E' and \
            st[j].stats.channel[2:3] != 'N'): 
          
            if(len(toPurge) == 0):
               toPurge.append(st[j].stats.station)
            else:
               if (st[j].stats.station != toPurge[-1]):
                   toPurge.append(st[j].stats.station)

    st = purgeListStation(st,toPurge,'r')

    return st



def purgeListStation(st,args,ty):
     
    new=Stream()
     
    if(ty != 'r'):
      ra=args.range.split()
      li=args.purge.split()
      se=args.set.split()
      de=args.azi.split()
      lii = []
    else:
#     ra=args.split()
#     li=args.split()
#     se=args.split()
#     de=args.split()
      li = args
      lii = []

    # select for distances
    if(ty=='d'):
      for i in range(len(st)):
          if(st[i].stats.dist >= float(ra[0]) and st[i].stats.dist <= float(ra[1])):
            new.append(st[i])

    # select station from azimuth
    if(ty=='a'):
       if(float(de[1]) >= float(de[0])):
         for i in range(len(st)):
             if(st[i].stats.az >= float(de[0]) and st[i].stats.az <= float(de[1])):
               new.append(st[i])
       else:
         for i in range(len(st)):
             if(st[i].stats.az >= float(de[0]) and st[i].stats.az <= 360):
               new.append(st[i])  
         for i in range(len(st)):
             if(st[i].stats.az >= 0 and st[i].stats.az <= float(de[1])):
               new.append(st[i])  

    # purge station from list 
    if(ty=='p' or ty=='r'):
      for i in range(len(st)):
          nn=0
          for j in range(len(li)):
              if(st[i].stats.station == li[j]):
                nn=1
          if(nn==0):
              new.append(st[i])
    
    # select station from list
    if(ty=='s'):
      for i in range(len(se)):
          c=0
          for j in range(len(st)):
              if (st[j].stats.station == se[i]):
                 c=c+1
                 new.append(st[j])
                 if(c==3):
                    break
          if(c==0):
            lii.append(se[i])
      print "    Stations",lii,"not in stream"
 


    return new


def sortStream(st,args):

    new=Stream()
    keys = args.sort.split()
    for i in range(len(keys)):
        if(keys[i] == 'A'):
           keys[i] = 'az'
        if(keys[i] == 'D'):
           keys[i] = 'dist' 
        if(keys[i] == 'C'):
           keys[i] = 'station' 

    for _i in keys[::-1]:
       st.traces.sort(key=lambda x: x.stats[_i], reverse=False) 

    return st
    

def removeGaps(self, min_gap, max_gap, verbose="TRUE"):

    """
    Returns the Stream object without trace gaps/overlaps.
    :param min_gap: All gaps smaller than this value will be omitted. The
          value is assumed to be in seconds. Defaults to None.
    :param max_gap: All gaps larger than this value will be omitted. The
          value is assumed to be in seconds. Defaults to None.
    :param verbose: stdout traces removed. Default verbose=False
    """

    new=Stream()
    self.sort()
    gap_list = []

    # happend first element at the end of stream to consider also last element of 
    # self stream. This is neaded beacuase loop of len(self.traces) - 1, which 
    # ignore the last element of stream
    self.append(self[0])

    for _i in xrange(len(self.traces) - 1):

       # skip traces with different network, station, location or channel
       if self.traces[_i].id != self.traces[_i + 1].id:
          new.append(self.traces[_i])
          continue
       # different sampling rates should always result in a gap or overlap
       if self.traces[_i].stats.delta == self.traces[_i + 1].stats.delta:
          flag = True
       else:
          flag = False
       stats = self.traces[_i].stats
       stime = stats['endtime']
       etime = self.traces[_i + 1].stats['starttime']
       delta = etime.timestamp - stime.timestamp

       # Check that any overlap is not larger than the trace coverage
       if delta < 0:
             temp = self.traces[_i + 1].stats['endtime'].timestamp - \
                    etime.timestamp
             if (delta * -1) > temp:
                 delta = -1 * temp
       # Check gap/overlap criteria
       if min_gap and delta < min_gap:
             new.append(self.traces[_i])
             continue
       if max_gap and delta > max_gap:
             new.append(self.traces[_i])
             continue
       # Number of missing samples
       nsamples = int(round(fabs(delta) * stats['sampling_rate']))
       # skip if is equal to delta (1 / sampling rate)
       if flag and nsamples == 1:
             new.append(self.traces[_i])
             continue
       elif delta > 0:
             nsamples -= 1
       else:
             nsamples += 1

       gap_list.append([_i,stats['network'], stats['station'],
                             stats['location'], stats['channel'],
                             stime, etime, delta, nsamples])
       if verbose == "True" or verbose == "TRUE" or verbose == "true":
          print  "Removed because of gap: ",stats['network'],stats['station'],stats['channel'],stime,etime,delta, nsamples

    return new

def removeShortTraces(st,tolerance,Tb,Te):

    expected = Te - Tb

    noGap  = []
    yesGap = []
    nrGap  = []

    for i in range(len(st)):

        npts  = st[i].stats['npts']
        delta = st[i].stats['delta']
        length = npts * delta
        obtained   = length * 100 / expected

        if obtained <= tolerance:
           yesGap.append(st[i])
           nrGap.append(i)
        else:
           noGap.append(st[i])

    return noGap

def remove21Comp(st):


    nn = Stream()
    StazOb = getStationslist(st)
    lis    = []

    for i in range(len(StazOb)):
        c=0
        for j in range(len(st)):
            if(st[j].stats.station==StazOb[i]):
              c=c+1
            if(c==3):
              lis.append(st[j].stats.station)
              break

    for i in range(len(lis)):
        for j in range(len(st)):
           if(st[j].stats.station==lis[i]):
             nn.append(st[j])

    return nn

def getStationslist(self):

    list=[]
    n=0

    for i in range(len(self)):
       list.append(self[i].stats['station'])

    # to be shure are really ordered
    list.sort()

    # list of distances
    out  = []
    for i in range(len(self)-1):
       if(self[i].stats['station'] != self[i+1].stats['station']):
          out.append(self[i].stats['station'])
    out.append(self[-1].stats['station'])

    return (out)

