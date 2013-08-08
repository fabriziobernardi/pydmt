#!/usr/bin/env python
# encoding: utf-8

import argparse,sys,os.path

def parseMyLine():

  parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='TDMT utility\n------------\n')

  parser.add_argument('--fseed', default='None',help='fseed file name inclusive of path if path different than \. Default None')
  parser.add_argument('--mseed', default='None',help='mseed file name inclusive of path if path different than \. Default None')
  parser.add_argument('--sac', default='None',help='Path for sac binary files. Default None')
  parser.add_argument('--asc', default='None',help='Path for sac alpha files. Default None')
  parser.add_argument('--Paz', default='None',help='Path for \'Paz\' poles and zeros (sac file format) files. Default None')
  parser.add_argument('--Resp', default='None',help='Path for \'RESP\' posponse files. Default None')
  parser.add_argument('--stsC', default='None',help='File for station file coordinates. Default None.')
  parser.add_argument('--deco', default='Y',help='Take Paz and resp for deconvolution of observed.')
  ############################
  # ---- optionals
  #
  # data type, begin and length
  parser.add_argument('--ori', default='None', help='Event Origin Time (e.g.: 2011-07-25T12:30:00).')
  parser.add_argument('--len', default='None', help='Length of signal in seconds After event origin Time.')
  parser.add_argument('--pre', default='0', help='Length of signal in seconds Before event origin Time.')
  parser.add_argument('--dva',default='1', help='Data type: 1(displacement); 2(velocity); Default = 1')
  parser.add_argument('--area',default='Earthquake', help='Earthquake epicenter area string. Default string ="Earthquake"')

  # event information
  parser.add_argument('--epi', default='None', help='Epicenter coordinate and depth: Lat Lon Depth.')
  parser.add_argument('--title', default='Event', help='Title of the event to diplay in the plot. Default=Event.')
  

  # directories and data
  parser.add_argument('--outdir',default='data', help='directory for data extraction. Default=data')
  parser.add_argument('--wsac',default='N', help='Allowed N/Y. Write stations and component raw sac files before remove instrument into --outdir directory. Default=N')
 
  # earth model and greens parameters
  parser.add_argument('--model',default='None',help='Earth model for greenfunctions. Default=None. See README_MODEL.txt for details')
  parser.add_argument('--npts',default='1024',help='Number of points for greens. Power of 2. Default=1024')
  parser.add_argument('--delta',default='0.5',help='Sampling interval in seconds for greens. Default=0.5')
  parser.add_argument('--cpus',default='1',help='Number of CPU available for greens computation. Min=1, Max=4. Default=1')
  parser.add_argument('--rvel',default='8',help='Reduction velocity. Recommendet value = 8km/s. Default=8')

  # Filters and similar
  parser.add_argument('--bandpass',default='0', help='Bandpass filter "corners fimn fmax". No Defaults. E.g.: "2 0.01 0.1"')
  parser.add_argument('--highpass',default='0', help='Highpass filter "corners freq". No Defaults. E.g.: "2 0.01"')
  parser.add_argument('--lowpass',default='0', help='Lowpass filter "corners freq". No Defaults. E.g.: "2 0.1"')
  parser.add_argument('--zeroph',default='False', help='Zerophase for high, low and bandpass. True/False. Defaul:False')
  parser.add_argument('--taper',default='0.1', help='cos Taper. If taper=-1 no taper is applied. Defaults=0.1')
  parser.add_argument('--sim',default='PZs', help='Remove instrument method: PZs for poles and zeros, RESP, for RESP_ files. Default=PZs')
  parser.add_argument('--flim',default='0.002 0.005 0.5 1', help='Corner frequency for deconvolution filtering. Defaults 0.002 0.005 0.5 1')
  parser.add_argument('--deci',default='None', help='Decimation factor for sampling rate. Only integer decimation factor allowed. Default=None')
  parser.add_argument('--inter', default='Y', help='Interpolate data to correct samplingrate if sampling not correct after decimation [Y]/N. Artifact may occour. Warning about decimation use --war Y. Default=Y')
  parser.add_argument('--war', default='N', help='Warnings Y/[N]. Default=N')
  

  # Analysis 
  parser.add_argument('--DeltaInv',default='1.0',help='Delta for data and greens for MT inversion. Default=1.0 Hz')
  parser.add_argument('--mti',default='1 1 1 0 0 0 1e20',help='Mxx Myy Mzz Mxy Mxz Myz Mo. Default = 1 1 1 0 0 0 1e20') 
  parser.add_argument('--maxw',default='2.5',help='Maximum variance pro station allowed. maxw=1 means perfect fit. Default=2.5')
  parser.add_argument('--range',default='None',help='Min and Max distance range for stations to use in km. Default=None')
  parser.add_argument('--purge',default='None',help='Station list to purege. Default=None')
  parser.add_argument('--set',default='All',help='Station list to use. Default=All')
  parser.add_argument('--azi',default='0 360',help='Azimuth range to select station. Default=0 360')
  parser.add_argument('--nrIter',default='0',help='Max number of iteration allowed. nrIter=0->No iteration, -1=auto,just 1 shot. Default=0. See tutorial for details')
  parser.add_argument('--iso',default='0',help='iso=0 [isotropic component set to 0, else iso=1. Default=0')
  parser.add_argument('--vr',default='50',help='min VR (variance reduction) pro station allowed.Default=50')
  parser.add_argument('--zcor',default='uniso',help='Fix Zcor for all 3 components (iso) or for each component (uniso), Default=uniso')
  parser.add_argument('--clean',default='0', help='Quality data check and dataset clean. 0: none; 1: Remove stations with spikes. See webpage for method details. Default=0')
  parser.add_argument('--noise',default='-1', help='Remove noisy traces [0-1]. -1: off; 0: pure noise; 1: pure signal (autocross corrr). General 0.4 is a good value fom small and moderate earthquake. See webpage for method details. Default=-1')

  # Auto calibration
  # when --auto Y enabled, distance and filter enabled with respect to the magnitude.
  # this start after removed spikes
  parser.add_argument('--auto', default='N', help='Auto calibration of the main parameters: station distance and filtering depending on first magnitude iteration. Automatically set nrIter = -1. Run after removing spike. See options auto_Dmin, auto_Dmax, auto_lowC and auto_highC. Default [N]/Y')
  parser.add_argument('--auto_mw', default='3.6 4.0 4.8 5.5', help='Magnitude Mw limits for automatic inversion. Two Mw -> 3 automatic settings. Defaults =\"3.6 4.0 4.8 5.5\"')
  parser.add_argument('--auto_Dmin', default='40 50 80 100 200', help='Min distances for Mw[-inf 3.6[ = 50, Mw[3.6 4.0[ = 75, Mw[4.0 4.8[ = 100 , Mw[5.5 inf] = 250 when --auto enabled. Defaults =\"40 50 80 100 200\"')
  parser.add_argument('--auto_Dmax', default='150 180 250 400 800', help='Max distances for Mw[-inf 3.6[ = 180, Mw[3.6 4.0[ = 250, Mw[4.0 4.8[ = 400 , Mw[5.5 inf] = 800 when --auto enabled. Defaults =\"150 180 250 400 800\"')
  parser.add_argument('--auto_lowC', default='0.033 0.033 0.030 0.017 0.010', help='Low corners filter for Mw[-inf 3.6[ = 0.030, Mw[3.6 4.0[ = 0.033, Mw[4.0 4.8[ = 0.020 , Mw[5.5 inf] = 0.010 when --auto enabled. Defaults =\".033 0.033 0.030 0.017 0.010\"')
  parser.add_argument('--auto_highC', default='0.066 0.066 0.050 0.025 0.020', help='High corners filter for Mw[-inf 3.6[ = 0.100, Mw[3.6 4.0[ = 0.066, Mw[4.0 4.8[ = 0.050 , Mw[5.5 inf] = 0.020 when --auto enabled. Defaults =\"0.066 0.066 0.050 0.025 0.020\"')
  parser.add_argument('--auto_vr', default='15 30 45 55 65', help='VR for Mw[-inf 3.6[ = 20, Mw[3.6 4.0[ = 30, Mw[4.0 4.8[ = 40 , Mw[5.5 inf] = 60 when --auto enabled. Defaults =\"15 30 45 55 65\"')
  

  # plot
  parser.add_argument('--plt', default='Y',help='Plot waveformfit [Y]/N.')
  parser.add_argument('--pltname', default='plot',help='Plot name file. Default = plot.pdf')
  parser.add_argument('--pltndc', default='Y',help='Plot non douple couple[Y] / douple couple [N] focal mechanism. Default = Y')
  parser.add_argument('--sort', default='A D C',help='sort stations with respect to C:code D:distance A:azimuth. Default = A D C')


# # Verbose
  parser.add_argument('--verbose',default='0',help='Verbose option. [0]/1')



  ###################################
  # ---- parse
  args=parser.parse_args()


  # deci must be integer
  if  args.deci != "None":
      args.deci = int(args.deci)

  # apply automati setting for interations
  if(args.auto == "Y"):
     args.nrIter = '-1'
     args.clean  = '1'

  return args


def checkConsistency(self):

  if(self.zcor != "iso" and self.zcor != "uniso"):
    print "Set correct value for --zcor option"
    sys.exit()

  # legth synt:
  tot = eval(self.delta) * eval(self.npts)
  mel = eval(self.len) + eval(self.pre)
  if(tot < mel):
    print "Greens too short for your data windows. Increase npts(2048/4096/..) / decrease delta / cut window"
    sys.exit()

  if(self.bandpass != '0'):
     li = self.bandpass.split()
     if(li[2] <= li[1]):
        print "bandpass filter: Use nrCorners min_freq max_freq"
        sys.exit()

  if(self.fseed == "None" and self.stsC == "None"):
     print "\nYou need to provide a file with name and station coordinates"
     print "Use --stsC NameFile"
     print "File format: Code Net lat lon elevation"
     print "        e.g: DIX CH 46.0801 7.4082 0"
     sys.exit()

  if(self.epi == "None"):
     print "\nYou need to provide the epicenter location of the earthquake"
     print "use --epi \"lat lon depth\"   e.g:  --epi \"36.50 12.56 10\""
     sys.exit()
     
  if(self.ori == "None"):
     print "\nYou need to provide the earthquake epicenter time"
     print "use --ori YYY-MM-DDThh:mm:ss   e.g: --ori 2011-07-25T12:30:00"
     sys.exit()

  if(self.ori == "None"):
     print "\nYou need to provide the length of your window in seconds after origin time"
     print "use --len nrSeconds .g:  e.g.: --len 200"
     sys.exit()

  if self.model != "None":
     aa=os.path.exists(self.model)
     if aa == False:
        print "Model --model", self.model, "does not exists"
        sys.exit()
  else:
     print "Model file unknown. Use --model option"
     sys.exit()

  # test if rdseed exists
  if os.system("which rdseed") == 256:
     print "rdseed not found"
     sys.exit()
       
  # check if outdir exists
  if self.outdir != ".":
     if os.path.exists(self.outdir) == False:
        os.mkdir(self.outdir)

  # Write args_log
  args_list=str(self).replace('(',',').replace(')',',').split(',')
  name = 'args_options.log'
  fou = open(self.outdir + os.sep + name,'w+')
  for i in range(1,len(args_list)-1):
      fou.write(args_list[i] +"\n") 
  fou.close()

def applyAutoSettings(self, metaMT):

    self.nrIter = '-1'
    self.bandpass = '0'

    mw = ('%.2f' % (metaMT[1]))
   
    # split 
    vrs = self.auto_vr.split()
    mws = self.auto_mw.split()
    hic = self.auto_highC.split()
    loc = self.auto_lowC.split()
    mnD = self.auto_Dmin.split()
    mxD = self.auto_Dmax.split()

    if(len(mws) < 2 or len(mws) > 5):
       print "Only 2 (3 option segments) or 3 magnitude (4 option segments) or 4 magnitudes (5 option segments) limits allowed with this version. Exit"
       sys.exit()

    # 3 different settings
    if(len(mws)==2):

       if(mw < mws[0]):
         self.vr       = str(vrs[0])
         self.lowpass  = str('2 ' + hic[0])
         self.highpass = str('2 ' + loc[0])
         self.range    = str(mnD[0] + ' ' + mxD[0])

       
       if(mw > mws[0] and mw < mws[1]):
         self.vr       = str(vrs[1])
         self.lowpass  = str('2 ' + hic[1])
         self.highpass = str('2 ' + loc[1])
         self.range    = str(mnD[1] + ' ' + mxD[1])

       
       if(mw > mws[1]):  
         self.vr       = str(vrs[2])
         self.lowpass  = str('2 ' + hic[2])
         self.highpass = str('2 ' + loc[2])
         self.range    = str(mnD[2] + ' ' + mxD[2])

    # 4 different settings
    if(len(mws)==3):

       if(mw < mws[0]):
         self.vr       = str(vrs[0])
         self.lowpass  = str('2 ' + hic[0])
         self.highpass = str('2 ' + loc[0])
         self.range    = str(mnD[0] + ' ' + mxD[0])

       
       if(mw > mws[0] and mw < mws[1]):
         self.vr       = str(vrs[1])
         self.lowpass  = str('2 ' + hic[1])
         self.highpass = str('2 ' + loc[1])
         self.range    = str(mnD[1] + ' ' + mxD[1])

       if(mw > mws[1] and mw < mws[2]):
         self.vr       = str(vrs[2])
         self.lowpass  = str('2 ' + hic[2])
         self.highpass = str('2 ' + loc[2])
         self.range    = str(mnD[2] + ' ' + mxD[2])
       
       if(mw > mws[2]):  
         self.vr       = str(vrs[3])
         self.lowpass  = str('2 ' + hic[3])
         self.highpass = str('2 ' + loc[3])
         self.range    = str(mnD[3] + ' ' + mxD[3])

    # 4 different settings
    if(len(mws)==4):

       if(mw < mws[0]):
         self.vr       = str(vrs[0])
         self.lowpass  = str('2 ' + hic[0])
         self.highpass = str('2 ' + loc[0])
         self.range    = str(mnD[0] + ' ' + mxD[0])

       
       if(mw > mws[0] and mw < mws[1]):
         self.vr       = str(vrs[1])
         self.lowpass  = str('2 ' + hic[1])
         self.highpass = str('2 ' + loc[1])
         self.range    = str(mnD[1] + ' ' + mxD[1])

       if(mw > mws[1] and mw < mws[2]):
         self.vr       = str(vrs[2])
         self.lowpass  = str('2 ' + hic[2])
         self.highpass = str('2 ' + loc[2])
         self.range    = str(mnD[2] + ' ' + mxD[2])
       
       if(mw > mws[2] and mw < mws[3]):
         self.vr       = str(vrs[3])
         self.lowpass  = str('2 ' + hic[3])
         self.highpass = str('2 ' + loc[3])
         self.range    = str(mnD[3] + ' ' + mxD[3])
       
       if(mw > mws[3]):  
         self.vr       = str(vrs[4])
         self.lowpass  = str('2 ' + hic[4])
         self.highpass = str('2 ' + loc[4])
         self.range    = str(mnD[4] + ' ' + mxD[4])


    print "---------------------------------------------"
    print "Automatic settings for magnitude Mw: ",mw
    print "Filter LowPass:   ",self.lowpass
    print "Filter HighPass:  ",self.highpass
    print "VR min:           ",self.vr
    print "Distance Range:   ",self.range
    print "---------------------------------------------"
    print "\n"

   
    return self
