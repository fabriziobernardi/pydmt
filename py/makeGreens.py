def generateGreens(st,args):

    import os
    import sys
    import shutil
    from subprocess import call

    delta = float(args.delta)
    nrcpu = float(args.cpus)
    npts  = float(args.npts)
    rvel  = float(args.rvel)

    # Load earth model
    # infile
    fin = open(args.model,'r')
    # outfilename
    outModel = args.outdir + os.sep + 'MODEL1'

    # read source earth model
    a=fin.readlines()
    z=[]
    max_dpth_model = 0
    for line in a:
        g=line.rsplit(' ')
        max_dpth_model += float(g[0])
        z.append(g[0])
    fin.close()
    max_dpth_model -= float(g[0])

    # get depth value from args
    epi = args.epi.rsplit(' ')
    depth = float(epi[2])

    # if requiede depth too deep for model
    if(depth > max_dpth_model):
      print "Source too deep for your earth model. Exit!"
      sys.exit()

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
         newFileModel.append("%8.4f%9.1f%10.1f" % (Dkms[i],0,rvel))

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
    diss=[]
    n=0

    for i in range(len(self)):
       list.append(self[i].stats['station'])
       diss.append(self[i].stats['dist'])

    list.append(list[0])
    diss.append(diss[0])

    # list of distances
    Dkms = []
    for i in range(len(list)-1):
        if(list[i] != list[i+1]):
           Dkms.append(diss[i])
           n = n+1
        
    if(n==0):
       n=1
       Dkms.append(diss[0])
    return (n,Dkms)

def findLyersInModel(z,depth):

    sum=0
    nrs=0
    nrb=0
    for i in range(len(z)):
        if sum <= depth:
           sum = sum + float(z[i])      
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
         sum = sum + float(z[i])

    a = sum - depth 
    b = float(z[i]) - a 
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
           sum = sum + float(z[i])

       for i in range(nrs[0]-1):
           new_l.append(l[i])
    
       # extract depth from lyer nrs[0]-1
       # and splito into two lyers according to depth source
       foe = l[nrs[0]-1].rsplit(' ')
      
       # new thinkes top lyer
       lyer_bel_z = sum - depth
       lyer_top_z = float(foe[0]) - lyer_bel_z

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
#   ra_vel=[0.6,0.5]    # less than ralay veleocity of the model
    ra_vel=[2.6,2.5]    # less than ralay veleocity of the model
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
    out.append("%5d" % (nrs[0]+1))               # line M+1
                                                 # Here to correct previous nrs[1] which was made for
                                                 # "Line 11 gives the layer number below the source"
                                                 # in the manual, but ...
                                                 # "Line 11 gives the layer number above the source + 1 layer"
    out.append("%15.7E%14.6E%10d" % (400,1.5,0)) # line M+2
    # line M+3
    out.append("%5d%9.1f%9.1f%9.1f%10.1f" % (nr_staz,ph_vel[0],ph_vel[1],ra_vel[0],ra_vel[1]))

    return out
