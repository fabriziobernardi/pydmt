#!/usr/bin/env python
# encoding: utf-8

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from obspy.imaging.beachball import Beachball, Beach
from obspy.imaging.beachball import PrincipalAxis
from matplotlib.backends.backend_pdf import PdfPages
from processData import getSTationslist
from stressRegime import stressRegime
from obspy.core.event import Event
from obspy.core.event import Catalog
from obspy.core.event import ResourceIdentifier
from obspy.core.event import CreationInfo
from obspy.core.event import EventDescription
from obspy.core.event import FocalMechanism
from obspy.core.event import NodalPlane, NodalPlanes
from obspy.core.event import PrincipalAxes, Axis
from obspy.core.event import MomentTensor, Tensor
from obspy.core.event import DataUsed
from obspy.core.event import Magnitude
from obspy.core.event import Origin
from obspy.core.event import CreationInfo
import datetime



def plotMTfit(st,gr,mt,Tb,Te,args):

    # length 
    tim=int(float(args.len)/float(args.DeltaInv))


    # set axis length
    tBegin = 0
    tEnd   = st[0].stats.endtime.timestamp \
                - st[0].stats.starttime.timestamp
    t      = np.arange(Tb,Te,st[0].stats.delta)
    nt     = len(t)
    ntr=len(st[0].data)

    # define number of stations:
    StazOb = getSTationslist(st)

    # subplot positions
    Xsize = 8.27  # default A4 size in inc.
    Ysize = 11.69

    # set np1 and title and catalog and get Stress regime
    title,np1,mkl = setNP1Title(StazOb, mt, '1e20', args)
    mtCatalog = makeCatalog(StazOb, mt, '1e20', args)

    # number of figures
    F = int(len(StazOb)/10) + 1
    i=0

    # loop over figures F
    for f in range(F):

        p = 3
        u = 12
        nrPage = 'Page ' + str(f+1) + '/' + str(F)
        page   = str(f+1) + '-' + str(F)
        name = args.pltname + page
 
        # number of stations for each plot
        if((f+1) < F):
           nr = 10
        else:
           nr = len(StazOb) - (F-1)*10
        nr *= 3


#       #  Initialize figure
        fig=plt.figure(f,figsize=(Xsize,Ysize),facecolor='w')

#       # write info
        plt.figtext(0.45, 0.98, nrPage, fontsize=8)
        plt.figtext(0.20, 0.96, 'Tan', fontsize=12)
        plt.figtext(0.50, 0.96, 'Rad', fontsize=12)
        plt.figtext(0.80, 0.96, 'Ver', fontsize=12)
        plt.figtext(0.05, 0.06, title, fontsize=10)

        p=3
               
        for k in range(0,nr):
             p=p+1
             d = "%s,%s,%s" % (u,3,p)
             Mmo = max(abs(st[i].data))
             Mmg = max(abs(gr[i].data))

             tt=gr[i].data[0:tim]
             ss=st[i].data[0:tim]
             ax1=fig.add_subplot(u,3,p)
             ax1=plt.plot(t,tt,color='r',linestyle='--')
             ax1=plt.plot(t,ss,color='k')
             if(args.zcor == 'uniso'):
              info = gr[i].stats.station + '  ' + \
               '$\Delta$ = ' + str(int(gr[i].stats.dist)) + ' [km]   ' + \
               '$\phi$ = '+ str(int(gr[i].stats.az)) + '$^\circ$   ' + \
               'Tcor = ' + str(int(st[i].stats.Tcor)) + '   ' + \
               'Rcor = ' + str(int(st[i].stats.Rcor)) + '   ' + \
               'Vcor = ' + str(int(st[i].stats.Vcor)) + '   ' + \
               'VR = ' + str(int(st[i].stats.VR))
             else:
              info = gr[i].stats.station + '  ' + \
               '$\Delta$ = ' + str(int(gr[i].stats.dist)) + ' [km]   ' + \
               '$\phi$ = '+ str(int(gr[i].stats.az)) + '$^\circ$   ' + \
               'Zcor = ' + str(int(st[i].stats.Zcor)) + '   ' + \
               'VR = ' + str(int(st[i].stats.VR))
             if(gr[i].stats.channel=='TAN'):
                ax1=plt.title(info,fontsize=8,position=(0.82,0.95)) 
             ax1=plt.axis('off')
             i=i+1
       
           # plot focal mechanism
        u = 10 
        np11 = [float(np1[0]),float(np1[1]),float(np1[2])+360.]
        ax0=plt.subplots_adjust(top=1.0)
        ax0=fig.add_subplot(u,3,0)
        ax0 = pylab.axis('off')
        if(args.pltndc == 'Y'):
          try:
            beach=Beach(mkl,  xy=(0.5,0.5), width=0.97)  
          except:
            beach=Beach(np11,  xy=(0.5,0.5), width=0.97)  
        else:
          beach=Beach(np11,  xy=(0.5,0.5), width=0.97)  
        ax0 = plt.gca()
        ax0.add_collection(beach,autolim=True)
        ax0.set_aspect("equal")

        # save fig
        fig.savefig(args.outdir + os.sep + name + '.pdf')


    # output result
    print "\n\nMT SOLUTION:\n"
    print title,"\n"

    # output mt_inv.out, mt_inv.log and mt_inv.xml
    namelast = 'mt_inv.out'
    namehist = 'mt_inv.log'
    nameqxml = 'mt_inv.xml'

    # out 
    outMTlast = open(args.outdir + os.sep + namelast,'w')
    outMTlast.write(title + "\n")
    outMTlast.close()
    # log
    outMThist = open(args.outdir + os.sep + namehist,'a+')
    outMThist.write(title + "\n\n ------ \n")
    outMThist.close()
    # xml
    mtCatalog.write(args.outdir + os.sep + "mt_inv.xml", format="QUAKEML")
 
    

    if(args.plt == 'Y'):
      plt.show()

    return title

def makeCatalog(StazList, mt, scale, args):

    epi   = args.epi.rsplit()
    model = args.model.split(os.sep)
    NrSt  = len(StazList)
    NrCo  = NrSt*3
    (Fmin,Fmax) = getFreq(args)
    Tmin  = ('%.0f' % (1/Fmax))
    Tmax  = ('%.0f' % (1/Fmin))
    mo    = ('%.3e' % (mt[0]))
    mw    = ('%.2f' % (mt[1]))
    Pdc   = ('%.2f' % (float(mt[2])/100))
    Pclvd = ('%.2f' % (float(mt[3])/100))
   
    Tval  = ('%10.3e' % (mt[22]))
    Tplg  = ('%4.1f' % (mt[23]))
    Tazi  = ('%5.1f' % (mt[24]))
    Nval  = ('%10.3e' % (mt[25]))
    Nplg  = ('%4.1f' % (mt[26]))
    Nazi  = ('%5.1f' % (mt[27]))
    Pval  = ('%10.3e' % (mt[28]))
    Pplg  = ('%4.1f' % (mt[29]))
    Pazi  = ('%5.1f' % (mt[30]))

    STp1  = ('%5.1f' % (mt[31]))
    DPp1  = ('%4.1f' % (mt[32]))
    RAp1  = ('%6.1f' % (mt[33]))
    STp2  = ('%5.1f' % (mt[34]))
    DPp2  = ('%4.1f' % (mt[35]))
    RAp2  = ('%6.1f' % (mt[36]))
    var   = ('%.2f' % (mt[37]))
    qua   = ('%d'   % (mt[38]))
    mij   = [mt[4],mt[5],mt[6],mt[7],mt[8],mt[9]]

    mm0   = str('%10.3e' % (mij[0]))
    mm1   = str('%10.3e' % (mij[1]))
    mm2   = str('%10.3e' % (mij[2]))
    mm3   = str('%10.3e' % (mij[3]))
    mm4   = str('%10.3e' % (mij[4]))
    mm5   = str('%10.3e' % (mij[5]))
    # Aki konvention
    Mrr   = mm5
    Mtt   = mm0 
    Mff   = mm1
    Mrt   = mm3
    Mrf   = mm4
    Mtf   = mm2

    # stress regime
    A1 = PrincipalAxis(val=mt[22], dip=mt[23], strike=mt[24])
    A2 = PrincipalAxis(val=mt[25], dip=mt[26], strike=mt[27])
    A3 = PrincipalAxis(val=mt[28], dip=mt[29], strike=mt[30])

    (regime, sh) = stressRegime(A1, A2, A3)
    sh = ('%5.1f' % (sh))

    #### Build classes #################################
    #
    #Resource Id is the event origin time for definition
    
    res_id      = ResourceIdentifier(args.ori)
    nowUTC      = datetime.datetime.utcnow()
    info        = CreationInfo(author="pytdmt", version="2.4", creation_time=nowUTC)
    evOrigin    = Origin(resource_id    = res_id,
                         time           = args.ori,
                         latitude       = epi[0],
                         longitude      = epi[1],
                         depth          = epi[2],
                         earth_model_id = model[-1],
                         creation_info  = info)
    # Magnitudes
    magnitude = Magnitude(mag=mw, magnitude_type = "Mw")
    # Nodal Planes
    np1      = NodalPlane(strike=STp1, dip=DPp1, rake=RAp1)
    np2      = NodalPlane(strike=STp2, dip=DPp2, rake=RAp2)
    planes   = NodalPlanes(nodal_plane_1=np1, nodal_plane_2=np2)
    # Principal axes
    Taxe     = Axis(azimuth=Tazi, plunge=Tplg, length=Tval)
    Naxe     = Axis(azimuth=Nazi, plunge=Nplg, length=Nval)
    Paxe     = Axis(azimuth=Pazi, plunge=Pplg, length=Pval)
    axes     = PrincipalAxes(t_axis=Taxe, p_axis=Paxe, n_axis=Naxe)
    # MT elements
    MT       = Tensor(m_rr=Mrr, m_tt=Mtt, m_pp=Mff, 
                      m_rt=Mrt, m_rp=Mrf, m_tp=Mtf)
    # Stress regime
    regStr   = 'Stress regime: ' + regime + ' -  SH = ' + sh
    strDes   = EventDescription(regStr)
    # MT dataset
    dataInfo = DataUsed(wave_type          = "combined",
                      station_count        = NrSt,
                      component_count      = NrCo,
                      shortest_period      = Tmin,
                      longest_period       = Tmax)
    source = MomentTensor(data_used        = dataInfo,
                      scalar_moment        = mo,
                      tensor               = MT,
                      variance_reduction   = var,
                      double_couple        = Pdc,
                      clvd                 = Pclvd,
                      iso                  = 0)
    focMec      = FocalMechanism(moment_tensor        = source,
                      nodal_planes         = planes,
                      principal_axes       = axes,
                      azimuthal_gap        = -1)

    #Initialize Event Catalog
    mtSolution = Event(creation_info=info)
    mtSolution.origins.append(evOrigin)
    mtSolution.magnitudes.append(magnitude)
    mtSolution.focal_mechanisms.append(focMec)
    mtSolution.event_descriptions.append(strDes)

    cat = Catalog()
    cat.append(mtSolution)

    return cat

def setNP1Title(NrSt, mt, scale, args):

    #cale = '1e+20'
    epi = args.epi.rsplit()
    mo    = ('%.3e' % (mt[0]))
    mw    = ('%.2f' % (mt[1]))
    Pdc   = ('%.2f' % (mt[2]))
    Pclvd = ('%.2f' % (mt[3]))
   
    Tval  = ('%10.3e' % (mt[22]))
    Tplg  = ('%4.1f' % (mt[23]))
    Tazi  = ('%5.1f' % (mt[24]))
    Nval  = ('%10.3e' % (mt[25]))
    Nplg  = ('%4.1f' % (mt[26]))
    Nazi  = ('%5.1f' % (mt[27]))
    Pval  = ('%10.3e' % (mt[28]))
    Pplg  = ('%4.1f' % (mt[29]))
    Pazi  = ('%5.1f' % (mt[30]))

    STp1  = ('%5.1f' % (mt[31]))
    DPp1  = ('%4.1f' % (mt[32]))
    RAp1  = ('%6.1f' % (mt[33]))
    STp2  = ('%5.1f' % (mt[34]))
    DPp2  = ('%4.1f' % (mt[35]))
    RAp2  = ('%6.1f' % (mt[36]))
    var   = ('%.2f' % (mt[37]))
    qua   = ('%d'   % (mt[38]))
    mij   = [mt[4],mt[5],mt[6],mt[7],mt[8],mt[9]]

    mm0   = str('%10.3e' % (mij[0]))
    mm1   = str('%10.3e' % (mij[1]))
    mm2   = str('%10.3e' % (mij[2]))
    mm3   = str('%10.3e' % (mij[3]))
    mm4   = str('%10.3e' % (mij[4]))
    mm5   = str('%10.3e' % (mij[5]))

    # Define axis for stress regime
    A1 = PrincipalAxis(val=mt[22], dip=mt[23], strike=mt[24])
    A2 = PrincipalAxis(val=mt[25], dip=mt[26], strike=mt[27])
    A3 = PrincipalAxis(val=mt[28], dip=mt[29], strike=mt[30])
 

    (regime, sh) = stressRegime(A1, A2, A3)
    sh = ('%5.1f' % (sh))

    (Fmin,Fmax) = getFreq(args)

    tim=int(float(args.len)/float(args.DeltaInv))

    s1  = args.title + '\n'
    s2  = args.ori + '  ' + 'lat: ' + str(epi[0]) + '  lon: ' + str(epi[1]) + '  depth: ' + str(epi[2]) + '\n' 
    s3  = 'Mxx = ' + mm0 + '   Myy = ' + mm1 + '   Mzz = ' + mm5 + '\n'
    s4  = 'Mxy = ' + mm2 + '   Mxz = ' + mm3 + '   Myz = ' + mm4 + '\n'
    s5  = 'Planes (Str/Dip/Slp): ' + str(STp1) + ' / ' + str(DPp1) + ' / ' + str(RAp1) + \
          ' --- ' + str(STp2) + ' / ' + str(DPp2) + ' / ' + str(RAp2) + '\n' 
#   s6  = 'Plane2              : ' + str(STp2) + ' / ' + str(DPp2) + ' / ' + str(RAp2) + '\n' 
    s7  = 'T-Axe  (Val/Plg/Azi): ' + str(Tval) + ' / ' + str(Tplg) + ' / ' + str(Tazi) + '\n'
    s8  = 'N-Axe               : ' + str(Nval) + ' / ' + str(Nplg) + ' / ' + str(Nazi) + '\n'
    s9  = 'P-Axe               : ' + str(Pval) + ' / ' + str(Pplg) + ' / ' + str(Pazi) + '\n'
    s10 = 'Mo =' + str(mo) + '[dyn*cm]' +  '  Mw = ' + str(mw) + '  Pdc: ' + Pdc + '  Pclvd: ' + Pclvd + '\n'
    s11 = 'Weighted reduced Variance[0-100] = ' + str(var) + '    Quality[0-4]: ' + str(qua) + '\n' 
    s12 = 'Nr stations: ' + str(len(NrSt)) + '   Frequency width: ' + str(Fmin) + '-' + str(Fmax) + ' [Hz]' + \
         '   Signal length: ' + str(tim) + ' [s]' + '\n'
    s13 = 'Stress regime: ' + regime + '      Max. horizontal stress orientation (SH): ' + sh

    inf = s1 + s2 + s3 + s4 + s5 + s7 + s8 + s9 + s10 + s11 + s12 + s13

    np1 = [float(STp1), float(DPp1), float(RAp1)]
    # xyz -> rtf convetion
    Mxyz = [float(mm5), float(mm0), float(mm1), float(mm3), -1.0*float(mm4), -1.0*float(mm2)]

    return (inf,np1,Mxyz)

def getFreq(args):

     hip = args.highpass
     lop = args.lowpass
     bdp = args.bandpass


     if hip != "0":
       elementh = hip.split()
       fmin = float(elementh[1])

     if lop != "0":
       elementl = lop.split()
       fmax = float(elementl[1])

     if bdp != "0":
       elements = bdp.split()
       fmin = float(elements[1])
       fmax = float(elements[2])
 
     return (fmax,fmin)

def realignZcor(gr,st,tim,args):


    for i in range(len(gr)):

        if  (gr[i].ststs.channel == "RAD"):
            zc=int(st[i].stats['Rcor'])
        elif(gr[i].ststs.channel == "TAN"):
            zc=int(st[i].stats['Tcor'])
        else:
            zc=int(st[i].stats['Zcor'])

        ze=np.zeros(zc)
       
        # if zc>=1 add before, remove end
        # if zc<=1 add end, renmove befor
        if(zc>=1):
          gr[i].data=np.concatenate((ze,gr[i].data), axis=0)
          gr[i].data=gr[i].data[0:tim]
      
    return gr
        
 


def plotWaves(tr,mode,kmrad,area,chan,nert,aziplot,cft,slta,outdir):

  if slta!="None":
    tas = slta.split(' ')
    sta = eval(tas[0])
    lta = eval(tas[1])
    thrOn=eval(tas[2])
    thrOff=eval(tas[3])
  else:
    thrOn=0
    thrOff=0


  if len(tr) == 0:
     print "No traces to plot"
     sys.exit() 
    
  # normalization to (kmrad[1]-kmrad[0])/len(tr)   
  normlization=(kmrad[1]-kmrad[0]) / (len(tr)/2)
  if mode ==2:
     for i in range(len(tr)):
         tr[i].data=(tr[i].data/max(abs(tr[i].data)))*normlization
  else:
     for i in range(len(tr)):
         tr[i].data=(tr[i].data/max(abs(tr[i].data)))

  #number of subplots
  nr_subpl = len(tr)

  #starting level for y_axe
  axe_y_lev=0
  xpoM1 = 0
  xpoM2 = 0
  xpoM3 = 0

  # create empty list
  mykmList=[None]*len(tr)

  (Amin,Amax)=aziplot.split(' ')
  Amin=eval(Amin)
  Amax=eval(Amax)


  #loop over nr_subpl
  # define figures to save
  fig1=plt.figure(1)
  if nert == "RT" or nert == "NE":
     fig2=plt.figure(2)
     fig3=plt.figure(3)
  for i in range(nr_subpl):

    # when --redo==Y: header valued of dist and az do not exist
    try:
       refaz = tr[i].stats['az']
    except:
       refaz = 0
       tr[i].stats['az']=refaz

    # define controller for plotting azimuth. This allows to plot
    # for example 20-50 or 330-20 cases
    pltAziController=0
    if mode >= 1:
       if Amin<=Amax:
          if refaz >= Amin and refaz <= Amax:
              pltAziController=1 
       else:
          if refaz <= Amin and refaz >= Amax:
             pltAziController=0
          else:
             pltAziController=1

    if pltAziController==1:

         # y-axe
         # when --redo==Y: header valued of dist and az do not exist
         try:
           dist=tr[i].stats['dist']
         except:
           dist= 10
           tr[i].stats['dist']=dist
           
         # x- axe
         tBegin = 0
         tEnd   = tr[i].stats.endtime.timestamp \
                - tr[i].stats.starttime.timestamp
         t      = np.arange(tBegin,tEnd,tr[i].stats.delta)
         nt     = len(t)
         ntr=len(tr[i].data)
         if nt == ntr:
            pass
         elif nt < ntr:
            tr[i].data=tr[i].data[:nt]
         elif ntr < nt:
            t=t[:ntr]

         # Labels and info
         xPosText = 0
         kPosText = t[len(t)-1]+0 
         distanceToPlot = "%8.2f" % (dist)
         nametoPlot     = "%-9s"   % (tr[i].stats['station'])
         li = list(tr[i].stats['channel'])
         if mode == 2:
            axe_y_lev = dist
            dista = int(tr[i].stats['dist'])
            azimu = int(tr[i].stats['az'])
            textdist = r"$\Phi$ = " + `azimu` + u"\u00b0"
            xlab = "Time [s]"
            ylab = "Distance [km]"
         else:
            textdist = ""
            xlab = "Time [s]"
            ylab = "Nr Traces"

#        try:
         if slta!="None":
            onOff = np.array(triggerOnset(cft[i].data, thrOn, thrOff))
#        except IndexError:
#               pass

         df = tr[i].stats.sampling_rate

         if mode == 1:
            pik=0.5
         else:
            pik=5.0

         if li[2] == "Z" and li[0] == chan:
              if mode == 1:
                 axe_y_lev = xpoM1
              plt.figure(1)
              plt.plot(t,tr[i].data+axe_y_lev,color='k')
              plt.text(xPosText, axe_y_lev, tr[i].stats['station'],fontsize=10,backgroundcolor='#99FFFF') #99FFFF
              plt.text(kPosText, axe_y_lev,textdist ,fontsize=8, backgroundcolor='#F0FFFF')
              plt.xlabel(xlab)
              plt.ylabel(ylab)
              plt.title(tr[i].stats['channel'])
              i=axe_y_lev-pik
              j=axe_y_lev+pik
              if slta!="None":
                 try:
                    plt.vlines(onOff[:, 0] / df, i, j, color='r', lw=2, label="Trigger On")
                 except IndexError:
                    pass
              xpoM1 += 1
         if nert == "RT":
           if li[2] == "T" and li[0] == chan:
              if mode == 1:
                 axe_y_lev = xpoM2
              plt.figure(2)
              plt.plot(t,tr[i].data+axe_y_lev,color='k')
              plt.text(xPosText, axe_y_lev, tr[i].stats['station'],fontsize=10,backgroundcolor='#99FFFF')
              plt.text(kPosText, axe_y_lev,textdist ,fontsize=8,backgroundcolor='#F0FFFF')
              plt.xlabel(xlab)
              plt.ylabel(ylab)
              plt.title(tr[i].stats['channel'])
              i=axe_y_lev-pik
              j=axe_y_lev+pik
              if slta!="None":
                 try:
                    plt.vlines(onOff[:, 0] / df, i, j, color='r', lw=2, label="Trigger On")
                 except IndexError:
                    pass
              xpoM2 += 1
           elif li[2] == "R" and li[0] == chan:
              if mode == 1:
                 axe_y_lev = xpoM3
              plt.figure(3)
              plt.plot(t,tr[i].data+axe_y_lev,color='k')
              plt.text(xPosText, axe_y_lev, tr[i].stats['station'],fontsize=10,backgroundcolor='#99FFFF')
              plt.text(kPosText, axe_y_lev,textdist ,fontsize=8,backgroundcolor='#F0FFFF')
              plt.xlabel(xlab)
              plt.ylabel(ylab)
              plt.title(tr[i].stats['channel'])
              i=axe_y_lev-pik
              j=axe_y_lev+pik
              if slta!="None":
                 try:
                    plt.vlines(onOff[:, 0] / df, i, j, color='r', lw=2, label="Trigger On")
                 except IndexError:
                    pass
              xpoM3 += 1
           else:
              pass
         elif nert == "NE":
           if li[2] == "N" and li[0] == chan:
              if mode == 1:
                 axe_y_lev = xpoM2
              plt.figure(2)
              plt.plot(t,tr[i].data+axe_y_lev,color='k')
              plt.text(xPosText, axe_y_lev, tr[i].stats['station'],fontsize=8,backgroundcolor='#99FFFF')
              plt.text(kPosText, axe_y_lev,textdist ,fontsize=8,backgroundcolor='#F0FFFF')
              plt.xlabel(xlab)
              plt.ylabel(ylab)
              plt.title(tr[i].stats['channel'])
              i=axe_y_lev-pik
              j=axe_y_lev+pik
              if slta!="None":
                 try:
                    plt.vlines(onOff[:, 0] / df, i, j, color='r', lw=2, label="Trigger On")
                 except IndexError:
                    pass
              xpoM2 += 1
           elif li[2] == "E" and li[0] == chan:
              if mode == 1:
                 axe_y_lev = xpoM3
              plt.figure(3)
              plt.plot(t,tr[i].data+axe_y_lev,color='k')
              plt.text(xPosText, axe_y_lev, tr[i].stats['station'],fontsize=8,backgroundcolor='#99FFFF')
              plt.text(kPosText, axe_y_lev,textdist ,fontsize=8,backgroundcolor='#F0FFFF')
              plt.xlabel(xlab)
              plt.ylabel(ylab)
              plt.title(tr[i].stats['channel'])
              i=axe_y_lev-pik
              j=axe_y_lev+pik
              if slta!="None":
                 try:
                    plt.vlines(onOff[:, 0] / df, i, j, color='r', lw=2, label="Trigger On")
                 except IndexError:
                    pass
              xpoM3 += 1
           else:
              pass

  fig1.savefig(outdir + os.sep + 'plotWavesZ.pdf')
  if nert == "NE":
     fig2.savefig(outdir + os.sep + 'plotWavesNS.pdf')
     fig3.savefig(outdir + os.sep + 'plotWavesEW.pdf')
  if nert == "RT":
     fig2.savefig(outdir + os.sep + 'plotWavesT.pdf')
     fig3.savefig(outdir + os.sep + 'plotWavesR.pdf')
  plt.show()
