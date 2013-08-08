#########################################################################
#                                                                       #
# Set of python function to load observed data into a stream            #
#                                                                       #
#########################################################################

import sys,os,glob
from obspy.core import Stream,read


def LoadObserved(args):

  ###################################################################
  # Simple load of data from different sources into different streams
  # fseed
  if(args.fseed != "None"):
    try:
       fse = read(args.fseed)
    except:
       print args.fseed + " Not found"
       sys.exit()

  # mseed
  if(args.mseed != "None"):
    try:
       mse = read(args.mseed)
    except:
       print args.mseed + " Not found"
       sys.exit()

  # sac binary
  if(args.sac != "None"):
    try:
       sac = read(os.path.join(args.sac + os.sep , '*.sac')) 
    except:
       try:
          sac = read(os.path.join(args.sac + os.sep , '*.SAC')) 
       except:
          print "sac/SAC files not found"
          sys.exit()
   
  # sac alpha
  if(args.asc != "None"):
    try:
       asc = read(os.path.join(args.asc + os.sep , '*.asc')) 
    except:
       try:
          asc = read(os.path.join(args.asc + os.sep , '*.ASC')) 
       except:
          print "asc/ASC files not found"
          sys.exit()

  # Extract resp and paz and associate names to stream 

  # initialize resfiles, pazfiles, ans statsions
  resfiles = [-1]
  pazfiles = [-1]
  stations = [-1]

  
  # FSEED
  if(args.fseed != "None"):
    # explode fseed and extract Paz and Resp
    (resfseed,pazfseed,stations) = extractResponse(args.fseed, 3,args)
    # associate names to each trace of stream
    # insert stats names PAZ_file, RESP_file
    # this is then neaded for later deconvolution
    # maybe no PZs or RESP file. So initialize stats.value with None
    # Epicenter coordinates are also here initialized into stats
    fse = initStats(fse,args)
    fse = associateNamesFilesToTraces(fse,resfseed,pazfseed,stations)

      
  # MSEED
  if(args.mseed != "None"):
    # aquire resp, pzs and stations files
    if(args.Resp != "None"):
      resfiles = findPazFiles(args.Resp,1)
    if(args.Paz != "None"):
      pazfiles = findPazFiles(args.Paz,2)
    if(args.stsC != "None"):
      stations = args.stsC 
    # associate
    mse = initStats(mse,args)
    mse = associateNamesFilesToTraces(mse,resfiles,pazfiles,stations)

  # SAC bin
  if(args.sac != "None"):
    # aquire resp, pzs and stations files
    if(args.Resp != "None"):
      resfiles = findPazFiles(args.Resp,1)
    if(args.Paz != "None"):
      pazfiles = findPazFiles(args.Paz,2)
    if(args.stsC != "None"):
      stations = args.stsC 
    # associate
    sac = initStats(sac,args)
    sac = associateNamesFilesToTraces(sac,resfiles,pazfiles,stations)

  # SAC alpha
  if(args.asc != "None"):
    # aquire resp, pzs and stations files
    if(args.Resp != "None"):
      resfiles = findPazFiles(args.Resp,1)
    if(args.Paz != "None"):
      pazfiles = findPazFiles(args.Paz,2)
    if(args.stsC != "None"):
      stations = args.stsC 
    # associate
    asc = initStats(asc,args)
    asc = associateNamesFilesToTraces(asc,resfiles,pazfiles,stations)

   
  # Combine all streams
  all = Stream()
  if(args.fseed != "None"):
    all = all + fse
  if(args.mseed != "None"):
    all = all + mse
  if(args.sac != "None"):
    all = all + sac
  if(args.asc != "None"):
    all = all + asc

  
  return all


def initStats(st,args):

    # here already initialize epicentral information
    a = args.epi.split(' ')

    for i in range(len(st)):
        st[i].stats.RESP_file = "None"
        st[i].stats.PZs_file  = "None"
        st[i].stats.stla      = -1000.0
        st[i].stats.stlo      = -1000.0
        st[i].stats.stev      = -1000.0
        st[i].stats.Zcor      = -1234.5
        st[i].stats.Rcor      = -1234.5
        st[i].stats.Tcor      = -1234.5
        st[i].stats.Vcor      = -1234.5
        st[i].stats.VR        = -1234.5
        st[i].stats.evla      = eval(a[0]) 
        st[i].stats.evlo      = eval(a[1]) 
        st[i].stats.depth     = eval(a[2])

    return st
        

def associateNamesFilesToTraces(st,res,paz,sts):

 
    for i in range(len(st)):
        sta = st[i].stats.station
        loc = st[i].stats.location
        com = st[i].stats.channel
        net = st[i].stats.network

        if(paz[0] != -1):
           paz_file = findFile('PAZ',sta,net,loc,com,paz)
           st[i].stats.PZs_file  = paz_file
        if(res[0] != -1):
           res_file = findFile('RES',sta,net,loc,com,res)
           st[i].stats.RESP_file = res_file
        if(sts[0] != -1):
           coord    = findcoor(sta,net,loc,com,sts)
           st[i].stats.stla      = coord[0]
           st[i].stats.stlo      = coord[1]
           st[i].stats.stev      = coord[2]
    
    return st

def findFile(mode,sta,net,loc,com,res):


    if(mode=="PAZ"):
      for i in range(len(res)):
        b=res[i].split(os.sep)
        a=b[-1].split('_')
        if(len(a)>1):
          if(a[0] == "SAC" and a[1] == "PZs" and a[2] == net and \
           a[3] == sta and a[4] == com):
          
#         sac files do not have location header value
#         if(a[0] == "SAC" and a[1] == "PZs" and a[2] == net and \
#          a[3] == sta and a[4] == com and a[5] == loc):
           return res[i] 
        else:
          a=b[-1].split('.')
          if(a[0] == "SAC" and a[1] == "PZs" and a[2] == net and \
           a[3] == sta and a[4] == com):
           return res[i]

    if(mode=="RES"):
      for i in range(len(res)):
        b=res[i].split(os.sep)
        a=b[-1].split('.')
        if(a[0] == "RESP" and a[1] == net and a[2] == sta and \
           a[3] == loc and a[4] == com):
           return res[i]
         
    return 0


def findcoor(sta,net,loc,com,sfi):

    sts = open(sfi)
    for line in sts:
      a = line.split( ) 
      if(len(a)>=1):
        if(a[0] == sta and a[1] == net):
           la = eval(a[2])
           lo = eval(a[3])
           el = eval(a[4])
           return(la,lo,el)
 
    return (0,0,0)
      
    
def extractResponse(seedFile,extractMode,args):

    # Call rdseed and extract station file, Respo and PAZ file
    # - seedFile must be a full seed
    # - extractMode = 0 --> nothoing to do (this func not called) 
    #                 1 --> RESP    
    #                 2 --> PAZ 
    #                 3 --> RESP and PAZ 
    #   string = rdseed  -f %s -S -R -p" %(data_dir, os.sep,data_id,data_id)

    resfiles = []
    pazfiles = []
    log = args.outdir + os.sep + 'rdseed.log'
    RES = args.outdir + os.sep + 'RESP' 
    PAZ = args.outdir + os.sep + 'PAZ'
    LOC = args.outdir + os.sep + 'LOC'

    # RESP FILES
    if extractMode == 1 or extractMode == 3:
       if os.path.exists(RES) == False:
          os.mkdir(RES)
       if os.path.exists(LOC) == False:
          os.mkdir(LOC)
       cmd_rdseed = "rdseed  -f %s -R  -q %s > %s" %(seedFile,RES,log)
       os.system(cmd_rdseed)

    # PAZs FILES
    if extractMode == 2 or extractMode == 3:
       if os.path.exists(LOC) == False:
          os.mkdir(LOC)
       if os.path.exists(PAZ) == False:
          os.mkdir(PAZ)
       cmd_rdseed = "rdseed  -f %s -p -q %s > %s" %(seedFile,PAZ,log)
       os.system(cmd_rdseed)

    # STATION COORDINATES
    cmd_rdseed = "rdseed  -f %s -S -q %s > %s" %(seedFile,LOC,log)
    stations = LOC + os.sep + 'rdseed.stations'
    os.system(cmd_rdseed)

    # load resp and paz file names
    resfiles = findPazFiles(RES,1)
    pazfiles = findPazFiles(PAZ,2)
    
  
    return(resfiles,pazfiles,stations)

def findPazFiles(pdir,mode):

    paz = []
    if(mode==1):
      for infile in glob.glob( os.path.join(pdir + os.sep ,'RESP.*') ):
        paz.append(infile)
    else:
      for infile in glob.glob( os.path.join(pdir + os.sep ,'SAC_PZs_*') ):
        paz.append(infile)
      if(len(paz) == 0):
        for infile in glob.glob( os.path.join(pdir + os.sep ,'SAC.PZs.*') ):
          paz.append(infile)


    return paz


def listStream(st):

    for i in range(len(st)):
        print st[i].stats.station, st[i].stats.channel, st[i].stats.location, \
        "   ",st[i].stats.starttime, st[i].stats.endtime
