####################################################################
#
# Script tdmt: Apply Time Domain MT convolution/inversion
#
####################################################################

import sys,os,copy
from loadData import LoadObserved
from processData import cleanStream,removeMeanTrend,cutWindow,rotateToGCP
from processData import selectData,dlazStream,getStationslist,sortStream
from processData import purgeListStation
from processData import cleanMultipleFragmentTraces,cleanZeroMaxTraces
from filterSeismograms import removeInstrument,filtering,decimate
from makeGreens import generateGreens
from readGreens import aquireGreens,upGrStats,reorderGreen
from makeSynt import makeSynt
from plotMTfit import plotMTfit
from invTDMT import invTDMT
from cleanSetMT import cleanSpike,cleanNoise,purgeStream
from myParser import applyAutoSettings
from processData import getSTationslist

def tdmt(args):

   # Begin data aquisition
   # ---------------------
   #
   # Load observed data: fseed, mseed, RESP, PAZ, Station_file
   print "... load data"
   dataStream = LoadObserved(args)
   #
   # compute distances and angles between event and station
   dataStream = dlazStream(dataStream)
   #
   # select data stations
   dataStream = selectData(dataStream,args) 
   if(len(dataStream)==0):
     print "no data found "
     sys.exit()
   #
   # Remove stations with multiple segment windows or less than 2 components
   dataStream = cleanMultipleFragmentTraces(dataStream,args) 
   #
   #here check if empty 0 line:
   dataStream = cleanZeroMaxTraces(dataStream,args)
   # cut windwow
   (dataStream,Tb,Te) = cutWindow(dataStream,'d',args)
   if(len(dataStream)==0):
     print "no data in stream"
     sys.exit()
   #
   # clean stream from traces with gaps, incomplete trace and
   # stations with 1 or 2 components only, or trace of zeros
   dataStream = cleanStream(dataStream,Tb,Te,args)
   if(len(dataStream)==0):
     print "Data pool empty. Exit."
     sys.exit()
   #
   # sort stream for station/distance/azimuth
   dataStream = sortStream(dataStream,args)
   #
   # remove trend and mean
   dataStream = removeMeanTrend(dataStream)
   #
   # End data aquisition

   # Begin removing instrument and filtering
   # ---------------------------------------
   #
   if(args.deco == "Y"):
     print "... removing instrument" 
     dataStream = removeInstrument(dataStream,args)
   #
   # here check if empty 0 line: this maybe because of zero line and/or wron PZ file
   dataStream = cleanZeroMaxTraces(dataStream,args)
   #
   print "... rotate to GCP, filtering and decimation"
   datastream = rotateToGCP(dataStream)
   #
   dataOrigin = copy.deepcopy(dataStream)
   dataStream = filtering(dataStream,'d',args)
   #
   dataStream = decimate(dataStream,'d',args)
   #
   (dataStream, noiseList) = cleanNoise(dataStream, args)
   #
   # check if dataSteram empty before to continue
   if(len(dataStream)==0):
     print "\n !! Empty dataStrem after cleaning !!"
     sys.exit()

   # Check if too many stations
   StList = getStationslist(dataStream)
   stationToPurge = []
   listToPurge    = []
   if(len(StList)>100):
       for l in range(100,len(StList)):
         listToPurge.append(l)
         stationToPurge.append(StList[l])
       dataStream  = purgeStream(dataStream,listToPurge)
       print "Only 100 stations allowed (FKRPROG!!). Take only first 100 station of the sorted stream.!"
       print "  Stations removed: ", stationToPurge,"\n"
     

   # Begin Greensfunction
   # --------------------
   # generate earth model
   #
   print "... generate greens"
   greenSpecFile = generateGreens(dataStream,args)
   #
   # greens (spectral domain) -> greenStream (time domain)
   print "... spectral to time domain"
   greenStream   = aquireGreens(greenSpecFile,args) 
    
   # update stats.station of greenStream
   print "... update stream header and sync with data"
   greenStream   = upGrStats(greenStream,dataStream)
   #
   # make rigth order of greens elements
   StazList = getStationslist(greenStream)
   greenStream = reorderGreen(greenStream,StazList)

   #
   # filtering
   print "... greens filtering"
   GreenOrigin = copy.deepcopy(greenStream)
   greenStream = filtering(greenStream,'g',args)
   greenStream = decimate(greenStream,'g',args)
   
   ################################################################
   # ---- MAKE TDMT
   # Run first MT RUN 
   print "\n... RUN TDMT\n"
   (greenStream, dataStream, synStream, MTx, metaMT, spin, rmStaz) = invTDMT(greenStream,dataStream,args)
   #
   ################################################################

   ################################################################
   # BEGIN CLEAN SPIKES
   # Clean dataset from spikes
   if (args.clean == "1" and len(dataStream) <21):
      print "Spike detection not possible, stations not enougth"

   elif(args.clean == "1" and len(dataStream) >=21):
      print "\n\nRun spike detection and clean\n"
   else:
      pass

   SpikeList = []
   while (args.clean == "1" and len(dataStream) >=21):

      Spikes = cleanSpike(dataStream, synStream, args)
      SpikeList.extend(Spikes)
      
      if (len(Spikes) > 0):
         
         listToRemove = ' '.join(['%d:%s' % (Spikes[n],dataStream[Spikes[n]*3].stats.station) for n in xrange(len(Spikes))])
         print "Station to clean because of spikes: ",listToRemove

         # Purge station from observed
         dataStream  = purgeStream(dataStream,Spikes)
         if(len(dataStream)==0):
            print args.title," No data left. No solution found. Exit"
            sys.exit()
         # Purge station from greens
         greenStream = purgeStream(greenStream,Spikes)

      else:
         args.clean = "0"

   # If spikes found, reapply mt inversion for correct MT estimation
   if (len(SpikeList)>0):
      (greenStream, dataStream, synStream, MTx, metaMT, spin, rmStaz) = invTDMT(greenStream,dataStream,args)
   #
   # END CLEAN SPIKES
   ################################################################


   ################################################################
   # BEGIN ITERATIONS
   if(args.nrIter != "0"):
     i=1
     while(spin==1):


        # ---------------------------------#
        # BEGIN SYNC FOR AUTOMATIC SETTINGS
        # Apply automatic settings. resync greens and data sets as after first iteration and Spike detection if applied
        # Apply only at first iteration !!!!!
        if(args.auto == "Y" and i ==1):

           print "\n\nApply Automatic settings"
           # APPLY NEW AUTOMATIC MW-DEPENDENT SETTINGS
           args = applyAutoSettings(args, metaMT)
       
           # BEGIN SYNC data and synt: sync the data and the green as original data set
           # DATA
           dataStream = copy.deepcopy(dataOrigin)          
           # remove noisy list removed up
           if(len(noiseList)>0):
              dataStream  = purgeStream(dataStream,noiseList)
              if (len(listToPurge) >0 ):
                 dataStream  = purgeStream(dataStream,listToPurge)
           # SYN
           greenStream = copy.deepcopy(GreenOrigin)
           #
           # Remove Spikes if the case
           if(len(SpikeList)>0):
             dataStream  = purgeStream(dataStream,SpikeList)
             greenStream = purgeStream(greenStream,SpikeList)
           #
           #
           # Now select station range
           dataStream  = purgeListStation(dataStream,args,'d')
           greenStream = purgeListStation(greenStream,args,'d')
           # Maybe with new distance range, no more data
           if(len(dataStream)==0):
             print args.title," No data left with new distance range. Exit"
             sys.exit()
           #
           #
           # Now filtering and decimation
           # DATA
           dataStream = filtering(dataStream,'d',args)
           dataStream = decimate(dataStream,'d',args)
           # SYN
           greenStream = filtering(greenStream,'g',args)
           greenStream = decimate(greenStream,'g',args)
           #
           # Apply first iteration MT
           # Run MT inv
           print "\n\nIteration Nr: 0"
           (greenStream, dataStream, synStream, MTx, metaMT, spin, rmStaz) = invTDMT(greenStream,dataStream,args)
        # END SYNC FOR AUTOMATIC SETTINGS           
        # ---------------------------------#
        
        print "\n\nIteration Nr: ",i
        # --------------------------
        # Output stations to remove from dataset
        listToRemove = ' '.join(['%d:%s' % (rmStaz[n],dataStream[rmStaz[n]*3].stats.station) for n in xrange(len(rmStaz))])
        print "Remove Station:",listToRemove

        # Purge station from observed
        dataStream  = purgeStream(dataStream,rmStaz)
        if(len(dataStream)==0):
          print args.title,"No data left. No solution found. Exit"
          sys.exit()
        # Purge station from greens
        greenStream = purgeStream(greenStream,rmStaz)

        #
        # Run MT inv
        (greenStream, dataStream, synStream, MTx, metaMT, spin, rmStaz) = invTDMT(greenStream,dataStream,args)

        # check if stop
        if(i>=int(args.nrIter) and args.nrIter != '-1'):
           spin = 0
        i=i+1

   ################################################################
   # ---- Plot
   print "\n\n",args.title," Solution OK"
   print "\n... PLOT MT fit"
   print "... cut window and realign for Zcor"
   synStream  = upGrStats(synStream,dataStream)
   #
   (synStream,Tb,Te) = cutWindow(synStream,'s',args)


   # Make plot
   FinalLog = plotMTfit(dataStream,synStream,metaMT,Tb,Te,args)
   #


   # Plot observed and synthetics in sac binary format
   if(args.wsac == "Y"):
     for i in range(len(synthStream)):
        name2 = args.outdir + os.sep + synStream[i].stats.station + "." + synStream[i].stats.channel + ".sac"
        synStream[i].write(name2, format="SAC")
     for i in range(len(dataStream)):
        name2 = args.outdir + os.sep + dataStream[i].stats.station + "." + dataStream[i].stats.channel + ".sac"
        dataStream[i].write(name2, format="SAC")


