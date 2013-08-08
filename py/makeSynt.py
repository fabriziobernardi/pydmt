def makeSynt(st,mti,mo,args):

    from obspy.core import Stream,Trace
    from math import cos,sin,radians
    import numpy as np
    import sys

    mxx = mti[0] 
    myy = mti[1] 
    mxy = mti[2]
    mxz = mti[3]
    myz = mti[4] 
    mzz = mti[5]
    mo  = 1.0*mo 
    MM  = mo/(1e20)
    MM  = 1.0

   
    # initialize
    synt = Stream()
    staz_list = []

    # npts
    npts  = st[0].stats.npts
    delta = st[0].stats.delta


    # aquire station list
    for i in range(len(st)):
        if (st[i].stats.channel == 'tss'): 
           staz_list.append(st[i].stats.station)
    
    # for each station of list: create 3 traces: ver,rad,tan
    # compute time elements for each one and update stats


    # make synt. Repeated 3 time the component loops over npts on syn.data 
    # because of some pointer confusion (mine)
    k=1
    for l in range(len(staz_list)):


        dati  = np.arange(npts)*0.0  #!!! initialize with floats
        az = radians(st[l*10+k].stats.az)

        ###############
        # TAN Component
        syn   = Trace(dati)
        for i in range(npts):

            syn.data[i] =    mxx*0.5*st[l*10+k*0].data[i]*sin(2*az) \
                           - myy*0.5*st[l*10+k*0].data[i]*sin(2*az) \
                           - mxy*1.0*st[l*10+k*0].data[i]*cos(2*az) \
                           - mxz*1.0*st[l*10+k*1].data[i]*sin(1*az) \
                           + myz*1.0*st[l*10+k*1].data[i]*cos(1*az) 

        # apply Mo
        syn.data=syn.data*MM*(-1)

        # update stats
        syn.stats.station = st[l*10].stats.station    
        syn.stats.channel = 'TAN'   
        syn.stats.az      = st[l*10].stats.az    
        syn.stats.baz     = st[l*10].stats.baz    
        syn.stats.dist    = st[l*10].stats.dist    
        syn.stats.gcarc   = st[l*10].stats.gcarc    
        syn.stats.evla    = st[l*10].stats.evla    
        syn.stats.evlo    = st[l*10].stats.evlo    
        syn.stats.stlo    = st[l*10].stats.stlo    
        syn.stats.stla    = st[l*10].stats.stla    
        syn.stats.delta   = delta

        # add to synt stream
        synt.append(syn) 


        ###############
        # RAD Component
        syn   = Trace(dati)
        for i in range(npts):

            syn.data[i] =    mxx*1/6*st[l*10+k*4].data[i]*(+1) \
                           - mxx*0.5*st[l*10+k*2].data[i]*cos(2*az) \
                           + mxx*1/3*st[l*10+k*8].data[i] \
                           + myy*1/6*st[l*10+k*4].data[i]*(+1) \
                           + myy*0.5*st[l*10+k*2].data[i]*cos(2*az) \
                           + myy*1/3*st[l*10+k*8].data[i] \
                           + mzz*1/3*st[l*10+k*8].data[i] \
                           - mzz*1/3*st[l*10+k*4].data[i]*(+1) \
                           - mxy*1.0*st[l*10+k*2].data[i]*sin(2*az) \
                           + mxz*1.0*st[l*10+k*3].data[i]*cos(1*az) \
                           + myz*1.0*st[l*10+k*3].data[i]*sin(1*az)

        # apply Mo
        syn.data=syn.data*MM*(-1)

        # update stats
        syn.stats.station = st[l*10].stats.station    
        syn.stats.channel = 'RAD'   
        syn.stats.az      = st[l*10].stats.az    
        syn.stats.baz     = st[l*10].stats.baz    
        syn.stats.dist    = st[l*10].stats.dist    
        syn.stats.gcarc   = st[l*10].stats.gcarc    
        syn.stats.evla    = st[l*10].stats.evla    
        syn.stats.evlo    = st[l*10].stats.evlo    
        syn.stats.stlo    = st[l*10].stats.stlo    
        syn.stats.stla    = st[l*10].stats.stla    
        syn.stats.delta   = delta

        # add to synt stream
        synt.append(syn) 


        ###############
        # VER Component
        syn   = Trace(dati)
        for i in range(npts):

            syn.data[i] =    mxx*1/6*st[l*10+k*7].data[i] \
                           - mxx*0.5*st[l*10+k*5].data[i]*(+1)*cos(2*az) \
                           + mxx*1/3*st[l*10+k*9].data[i] \
                           + myy*1/6*st[l*10+k*7].data[i] \
                           + myy*0.5*st[l*10+k*5].data[i]*(+1)*cos(2*az) \
                           + myy*1/3*st[l*10+k*9].data[i] \
                           + mzz*1/3*st[l*10+k*9].data[i] \
                           - mzz*1/3*st[l*10+k*7].data[i] \
                           - mxy*1.0*st[l*10+k*5].data[i]*(+1)*sin(2*az) \
                           + mxz*1.0*st[l*10+k*6].data[i]*(+1)*cos(1*az) \
                           + myz*1.0*st[l*10+k*6].data[i]*(+1)*sin(1*az)

        # apply Mo
        syn.data=syn.data*MM*(+1)

        # update stats
        syn.stats.station = st[l*10].stats.station   
        syn.stats.channel = 'VER'   
        syn.stats.az      = st[l*10].stats.az    
        syn.stats.baz     = st[l*10].stats.baz    
        syn.stats.dist    = st[l*10].stats.dist    
        syn.stats.gcarc   = st[l*10].stats.gcarc    
        syn.stats.evla    = st[l*10].stats.evla
        syn.stats.evlo    = st[l*10].stats.evlo
        syn.stats.stlo    = st[l*10].stats.stlo
        syn.stats.stla    = st[l*10].stats.stla
        syn.stats.delta   = delta

        # add to synt stream
        synt.append(syn) 
        
    return synt
    

def generateGreens(st,args):

    import os
    import sys
    import shutil
    from subprocess import call

    delta = eval(args.delta)
    nrcpu = eval(args.cpus)
    npts  = eval(args.npts)
    rvel  = eval(args.rvel)

    # Load earth model
    # infile
    fin = open(args.model,'r')
    # outfilename
    outModel = args.outdir + os.sep + 'MODEL1'

    # read source earth model
    a=fin.readlines()
    z=[]
    for line in a:
        g=line.rsplit(' ')
        z.append(g[0])
    fin.close()

    # get depth value from args
    epi = args.epi.rsplit(' ')
    depth = eval(epi[2])

    # get lyers from model
    (nr_lyers)   = findLyersInModel(z,depth)

    # rewrite lyers model with new depth
    (nr_lyers,g) = rewriteLyersModel(z,a,depth,nr_lyers)

    # get number of stations to generate Greens
    (nrStaz,Dkms) = getNrSTationsIntoStream(st)

    # Write first 13 lines of fkprog lInput file
    newFileModel = writeNewfileModel(nr_lyers,g,depth,nrcpu,npts,delta,nrStaz)

    # add station lines for each station to compute the Greens 
    for i in range(nrStaz):
         dKm = st[i].stats['dist']
         newFileModel.append("%8.2f%9.1f%10.1f" % (Dkms[i],0,rvel))

    # Write new model file for fkprog
    fou = open(outModel,'w+')
    for i in range(len(newFileModel)):
        fou.write(newFileModel[i]+"\n")
    fou.close()

    # remove old Greens file if exists and call fkprog to generate new greens
    # remove if exists
    if os.path.exists(args.outdir + os.sep + 'GREEN.1') == True:
       os.remove(args.outdir + os.sep + 'GREEN.1')
    if os.path.exists('GREEN.1') == True:
       os.remove('GREEN.1')

    # commands to esexute 1:cd outdir 2: fkprog
    os.system('FKRPROG < ' + outModel + ' > green.log')
    
    # move greens and log file
    if os.path.exists('green.log') == True:
       shutil.move('green.log',args.outdir + os.sep + 'green.log')
    if os.path.exists('GREEN.1') == True:
       shutil.move('GREEN.1',args.outdir + os.sep + 'GREEN.1')

    outgreen = args.outdir + os.sep + 'GREEN.1'

    return outgreen 


def getNrSTationsIntoStream(self):

    list=[]
    n=0

    for i in range(len(self)):
       list.append(self[i].stats['station'])

    # to be shure are really ordered
    list.sort()

    # list of distances
    Dkms = []
    for i in range(len(self)-1):
       if(self[i].stats['station'] != self[i+1].stats['station']):
          Dkms.append(self[i].stats['dist'])
          n=n+1

    return (n,Dkms)

def findLyersInModel(z,depth):

    sum=0
    nrs=0
    nrb=0
    for i in range(len(z)):
        if sum <= depth:
           sum = sum + eval(z[i])      
           nrs=nrs+1
        else:
           nrb=nrb+1
       
    return [nrs,nrb]

def rewriteLyersModel(z,l,depth,nrs):

    # z: array of depth from earth model (first colun)
    # l: lines of models
    # depth: depth of source
    # nrs: array of nr_lyers [0]: above source [1]: below of source

    for i in range(len(l)):
        l[i]=l[i].rstrip()

    out_l = []
  
    sum = 0
    for i in range(nrs[0]):
         sum = sum + eval(z[i])

    a = sum - depth 
    b = eval(z[i]) - a 
    if(b==0):
       out_l = l
       nrs[1]=nrs[1]+0
  
    else:  
       # reconstruct model
       # if b!=0
       # line 0 - nrs[0] ok
       # line nrs[0] split: top below with respect to the source depth 
       # add nrs[0]-nrs[1] lyers
       # update nrs[1] 

       new_l = []
       sum = 0
       for i in range(nrs[0]-0):
           sum = sum + eval(z[i])

       for i in range(nrs[0]-1):
           new_l.append(l[i])
    
       # extract depth from lyer nrs[0]-1
       # and splito into two lyers according to depth source
       foe = l[nrs[0]-1].rsplit(' ')
      
       # new thinkes top lyer
       lyer_bel_z = sum - depth
       lyer_top_z = eval(foe[0]) - lyer_bel_z

       lyer_top_tmp = l[nrs[0]-1]
       lyer_bel_tmp = l[nrs[0]-0]

       foe_top = lyer_top_tmp.rsplit(' ')
       foe_bel = lyer_top_tmp.rsplit(' ')

       foe_top[0] = "%.4E" % (lyer_top_z)
       foe_bel[0] = "%.4E" % (lyer_bel_z)

       foe_top = ' '.join(foe_top) 
       foe_bel = ' '.join(foe_bel) 

       new_l.append(foe_top)
       new_l.append(foe_bel)
  
       
       for i in range(nrs[1]):  
           new_l.append(l[nrs[0]+i])
           
       out_l=new_l
       nrs[1]=nrs[1]+1

    return nrs,out_l

def writeNewfileModel(nrs,l,depth,nrcpu,npts,delta,nr_staz):

    # Station parameters
    ph_vel=[10000,30]
    ra_vel=[2.9,2.5]    # less than ralay veleocity of the model
    x=6
    beg_freq=1
    nr_l = nrs[0]+nrs[1]
    nr_fr=int(npts/2)
    out = []

    # start writing model list
    out.append('.F.')                            # line 1
    out.append('    0   64')                     # line 2
    out.append("GREEN.%1d" % nrcpu)              # line 3
    # line 4
    out.append("%7.1f%10.2f%8d%5d%5d%9.3f%6d%5d" % (x,depth,beg_freq,nr_fr,npts,delta,nr_l,nrcpu))
    # line 5
    out.append("%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d" % (1,1,1,1,1,1,1,1,1,1,0))
    # lyers
    for i in range(len(l)):
        out.append(" %s" % (l[i]))
    out.append("%5d" % (nrs[1]))                 # line M+1
    out.append("%15.7E%14.6E%10d" % (400,1.5,0)) # line M+2
    # line M+3
    out.append("%5d%9.1f%9.1f%9.1f%10.1f" % (nr_staz,ph_vel[0],ph_vel[1],ra_vel[0],ra_vel[1]))

    return out
