from obspy.signal.cross_correlation import xcorr
from obspy.core import Stream
from obspy.imaging.beachball import MT2Axes,MomentTensor,MT2Plane,AuxPlane
from processData import getSTationslist
from makeSynt import makeSynt
from gaussj import gaussj
import math
import numpy as np
import sys,copy

def invTDMT(gr,st,args):

    # gf scaling factor
    gfscale = 1.0e+20; # Dyn * cm

    #trim data and synt for required legth
    (gr,st) = trimStreams(gr,st,args)

    # here start with cross correlation. 
    # Zcore returned into stats
    (gr,st) = correlatePhaseShift(gr,st,args)

    # Each trace (gr or st) includes Zcore values for relignment
    # using args.zcor = iso   -> realign with Zcor
    #       args.zcor = uniso -> realign with Tcor, Rcor, Vcor
    (gr,st) = realignPhaseSt(gr,st,args)

    # Initialize station weigths, and update
    # with respect station distance. The closer has 1
    # W=dist[km]/mindist[km]
    (gr,st) = setStationWeigths(gr,st)

    # maybe multiply gr6,7,8 | gr8,9 for iso=6 => -1

    # Set and normalize AtA and AIV matrix
    (AIV,B,AJ,W) = setATA(gr,st,args)
    # OK

    # Calculate Righthand Side
    B = righthandSide(st,B,AJ,W,args)
    # OK as B value

    # Solve for MT
    M = gaussj(AIV,B)

    # if isotropic constrain, set Mzz
    if(args.iso == "0"):
       zz = -1*(M[0]+M[1])
       M[5]=zz

    M *= -1.

    # *Convert deviatoric moment tensor to AKI convention*
    MTx = toXyz(gfscale,M)
    MTr = toRtf(gfscale,M)

    # Here compute Planes, Axes, Mo, Mw
    (Mo, Mw, Pdc, Pclvd, EigVal, EigVec, T, N, P, np1, np2) = mt2parameters(MTx,M,gfscale)

    # make synt
    sy = makeSynt(gr,M,Mo,args)
   
    # compute Variance
    (gr,st,sy,var,qual) = checkFit(gr,st,sy,args) 
 
    # check stations to purge 
    (nrToPurge, spin)=checkAndpurgeVR(st,gr,args)

    # PackMetadata
    metaMT = packMetaMT(M, Mo, Mw, Pdc, Pclvd, EigVal, EigVec, T, N, P, np1, np2, var, qual, gfscale)
    
    return (gr, st, sy, MTx, metaMT, spin, nrToPurge)
 

def packMetaMT(M, Mo, Mw, Pdc, Pclvd, EigVal, EigVec, T, N, P, np1, np2, var, qual, k):

    meta = np.zeros(39.)
    meta[0] = Mo 
    meta[1] = Mw 
    meta[2] = Pdc 
    meta[3] = Pclvd

    for i in range(0,6):
        meta[i+4] = M[i]*k 

    for i in range(0,3):
        meta[i+10] = EigVal[i] 

    flat = [x for sublist in EigVec for x in sublist]
    for i in range(0,9):
        meta[i+13] = flat[i] 

    meta[22] = T.val
    meta[23] = T.dip
    meta[24] = T.strike
    meta[25] = N.val
    meta[26] = N.dip
    meta[27] = N.strike
    meta[28] = P.val
    meta[29] = P.dip
    meta[30] = P.strike

    meta[31] = np1[0]
    meta[32] = np1[1]
    meta[33] = np1[2]
    meta[34] = np2[0]
    meta[35] = np2[1]
    meta[36] = np2[2]

    meta[37] = var
    meta[38] = qual

    return meta

def checkAndpurgeVR(st,gr,args):

    # get station list
    ls_Sta = getSTationslist(st)

    # Find max VAR
    max_W = float(args.maxw)
    min_V = float(args.vr)
    spin  = 0
    NOKc   = []    #station not ok with corr
    NOKe   = []    #station not ok with E

    (W,V) = getDistrErr(st,max_W,min_V,args)

    if(W>=max_W):
      max_W=W
    if(V<=min_V):
      min_V=V

    for i in range(len(st)/3):
        if (float(st[i*3+0].stats.VR) <= min_V):
            NOKe.append(i) 

    if len(NOKe)>=1:
       spin = 1
    else:
        for i in range(len(st)/3):
            if (float(st[i*3+0].stats.VR) <= float(args.vr)):
               NOKe.append(i)
    
    if len(NOKe)>=1:
       spin = 1
    else:
       spin = 0

    return NOKe,spin

def getDistrErr(st,maxW,minV,args):

    List_fW = []
    List_fZ = []
    List_fV = []
     
    for i in range(len(st)/3):
        List_fW.append(st[i*3+0].stats.W)
        List_fZ.append(st[i*3+0].stats.Zcor)
        List_fV.append(st[i*3+0].stats.VR)

    np.array(List_fW)
    np.array(List_fZ)
    np.array(List_fV)
 
    #mean and std
    meanW=np.mean(List_fW)
    meanZ=np.mean(List_fZ)
    meanV=np.mean(List_fV)

    stdW =np.std(List_fW)
    stdZ =np.std(List_fZ)
    stdV =np.std(List_fV)

    # if more than 7 station used make statistics, else use fix max values of arges
    if(len(List_fW) >= 7):
      # 68% larger
      maxW = meanW+stdW
      minW = meanW-stdW
      maxZ = int(meanZ+stdZ)
      minV = meanV-stdV

    vr_list = ' '.join(['%.2f' % (List_fV[n]) for n in xrange(len(List_fV))])
    print 'VRs: ',vr_list

    table = {'Mean':meanV,'Std':stdV,'Threshold':minV}
    print ('Mean: {0[Mean]:.2f}; Std: {0[Std]:.2f}; Threshold: {0[Threshold]:.2f};'.format(table))

    return (maxW,minV)

def checkFit(gr,st,sy,args):

    npts   = len(st[0].data)
    degree = float(args.iso)
    WSUM   = 0.0
    Etot   = 0.0
    VAR    = 0.0
    DVAR   = 0.0
    Dtot   = 0.0

    for i in range(len(st)/3):

        Dpower = 0.0;
        Etmp   = 0.0;
        E      = 0.0;

        # For 1 station
        for j in range(len(st[i*3+0].data)):
           Etmp = st[i*3+0].data[j] - sy[i*3+0].data[j]
           E += Etmp*Etmp;
           Etmp = st[i*3+1].data[j] - sy[i*3+1].data[j]
           E += Etmp*Etmp;
           Etmp = st[i*3+2].data[j] - sy[i*3+2].data[j]
           E += Etmp*Etmp;

           Dpower += st[i*3+0].data[j]**2
           Dpower += st[i*3+1].data[j]**2
           Dpower += st[i*3+2].data[j]**2

        # resume 1 station info
        WSUM += float(st[i*3+0].stats.W)
        Etot += E
        VAR  += float(st[i*3+0].stats.W) * E
        Dtot += Dpower
        DVAR += float(st[i*3+0].stats.W)*Dpower;
        E    /= Dpower;
        VR    = (1.0 - E)*100.0
        st[i*3+0].stats.VR = VR
        st[i*3+1].stats.VR = VR
        st[i*3+2].stats.VR = VR
         
    # end single stations checkfit loop
    # begin general solution
    vred  = (1.0-Etot)*100.0
    VAR  /= WSUM
    DVAR /= WSUM;
    VAR  /= DVAR;

    # THIS is he VR quantity of the orgininal code.
    # Keep this value for compatibility
    # The same if for station VR in the loop above
    VAR   = (1.0-VAR)*100.0

    # Set UQality
    if (VAR < 20.0):
       qual = 0
    elif (VAR >= 20.0 and VAR < 40.0):
       qual = 1
    elif (VAR >= 40.0 and VAR < 60.0):
       qual = 2
    elif (VAR >= 60.0 and VAR < 80.0):
       qual = 3
    else:
       qual = 4
   
    return (gr,st,sy,VAR,qual)
    

def mt2parameters(MTx,M,gfscale):

    #Iso Mo
    MoIso = (MTx[0][0] + MTx[1][1] + MTx[2][2])/3

    c = MomentTensor(MTx[2][2],MTx[0][0],MTx[1][1],MTx[0][2],-MTx[1][2],-MTx[0][1],gfscale) 

    # Principal axes
    (T, N, P) = MT2Axes(c)

    # Nodal planes
    np0 = MT2Plane(c)
    np2 = AuxPlane(np0.strike,np0.dip,np0.rake)
    # Convention rake: up-down
    if (np0.rake>180):
       np0.rake= np0.rake-360
    if (np0.rake<-180):
       np0.rake= np0.rake+360
    np1 = np.zeros(3.)
    np1 = [np0.strike,np0.dip,np0.rake]
    
    # Compute Eigenvectors and Eigenvalues
    # Seismic Moment and Moment Magnitude
    (EigVal, EigVec) = np.linalg.eig(MTx)
    b = copy.deepcopy(EigVal)
    b.sort()
    Mo = (abs(b[0]) + abs(b[2]))/2.
    Mw = math.log10(Mo)/1.5-10.73

    # Compute double-couple, CLVD & iso
    d = copy.deepcopy(EigVal)
    d[0]=abs(d[0]) 
    d[1]=abs(d[1]) 
    d[2]=abs(d[2]) 
    d.sort()
    eps=abs(d[0])/abs(d[2])
    pcdc=100.0*(1.0-2.0*eps)
    pcclvd=200.0*eps
    pcdc=pcdc/100.0
    pcclvd=pcclvd/100.0
    pciso=abs(MoIso)/Mo
    pcsum=pcdc+pcclvd+pciso
    pcdc=100.0*pcdc/pcsum
    pcclvd=100.0*pcclvd/pcsum
    pciso=100.0*pciso/pcsum
   
    Pdc   = pcdc
    Pclvd = pcclvd

    return (Mo, Mw, Pdc, Pclvd, EigVal, EigVec, T, N, P, np1, np2)
    
    
def getAxis(Val,Vec):

    radian=57.29577951
    hsq2=.70710678

    pl = np.zeros(3.)
    az = np.zeros(3.)

    for i in range(0,3):
       pl[i] = abs(math.degrees(math.atan2 (-1.*Vec[0][i], math.sqrt(Vec[1][i]**2 + Vec[2][i]**2)) ))
       az[i] = math.degrees(math.atan2 (Vec[2][i], Vec[1][i]))+180.

    return 1


def toRtf(k,M):
    
    Aki = np.zeros(shape=(3,3))
    Aki[0][0] = ( 1.0* k * M[0]) 
    Aki[0][1] = (-1.0* k * M[2]) 
    Aki[0][2] = ( 1.0* k * M[3]) 
    Aki[1][0] = (-1.0* k * M[2]) 
    Aki[1][1] = ( 1.0* k * M[1]) 
    Aki[1][2] = (-1.0* k * M[4]) 
    Aki[2][0] = ( 1.0* k * M[3]) 
    Aki[2][1] = (-1.0* k * M[4]) 
    Aki[2][2] = ( 1.0* k * M[5]) 

    return Aki
    
def toXyz(k,M):
    
    Aki = np.zeros(shape=(3,3))
    Aki[0][0] = ( 1.0* k * M[0]) 
    Aki[0][1] = ( 1.0* k * M[2]) 
    Aki[0][2] = ( 1.0* k * M[3]) 
    Aki[1][0] = ( 1.0* k * M[2]) 
    Aki[1][1] = ( 1.0* k * M[1]) 
    Aki[1][2] = ( 1.0* k * M[4]) 
    Aki[2][0] = ( 1.0* k * M[3]) 
    Aki[2][1] = ( 1.0* k * M[4]) 
    Aki[2][2] = ( 1.0* k * M[5]) 

    return Aki
    
def righthandSide(st,B,AJ,W,args):

    if (args.iso == "0"):
       isoflag = 5
    else:
       isoflag = 6

    cnt1=cnt2=cnt3=0
    tmp = np.zeros(10*st[0].stats.npts*st[0].stats.npts)
    for i in range(len(st)/3):
        Np=st[i].stats.npts
        Z =int(st[i].stats.Zcor)
        Z = 0
        cnt1 = cnt2 = cnt3
#       cnt2 = cnt3
        cnt2 += Np
        cnt3 += 2*Np
#       for j in range(0,Np):
        for j in range(Z,Np+Z):
            if j < len(st[i*3+0].data):
              tmp[cnt1] = st[i*3+0].data[j]
              tmp[cnt2] = st[i*3+1].data[j]
              tmp[cnt3] = st[i*3+2].data[j]
#             print("%d %d %g %g %g" % (i,j,st[i*3+0].data[j], st[i*3+1].data[j], st[i*3+2].data[j]))
            else:
              tmp[cnt1] = 0.0
              tmp[cnt2] = 0.0
              tmp[cnt3] = 0.0

            cnt1 +=1
            cnt2 +=1
            cnt3 +=1

    # To build B must realign j to Z
    nn=0
    for i in range(0,isoflag):
        for j in range(0+0,cnt3+0):
            B[i][0] += AJ[i][j-0] * tmp[j] * W[j]
            nn+=1

    return B
   

def setATA(gr,st,args):

    if (args.iso == "0"):
       isoflag = 5
    else:
       isoflag = 6

    #set W matrix
    W = np.zeros(len(st)*st[0].stats.npts)
    n = 0
    for i in range(len(st)):
        for j in range(st[i].stats.npts):
            W[n]=st[i].stats.W
            n+=1

    #set ATA matrix
    AIV = np.zeros(shape=(isoflag,isoflag))
    
    #set B matrix
    B = np.zeros(shape=(isoflag,2))

    #set AJ matrix
    # -- Initialize
    # number of elements for each row: 3* Nrstaz * npts
    Np = st[0].stats.npts 
    nsta = int(len(st)/3)
    nrElementsInRow = float(3 * Np * nsta)
    AJ = np.zeros(shape=(isoflag,nrElementsInRow))
    #
    # -- allocate gr and st
    cnt1=cnt2=cnt3=0
    for i in range(0,nsta):
        Np=st[i].stats.npts 
        Z = 0
        cnt1 = cnt2 = cnt3
#       cnt2 = cnt3 
        cnt2 += Np
        cnt3 += 2*Np
#       print i,Np,Z,cnt1,cnt2,cnt3
        for j in range(0,Np):
 
            alpha = st[i*3].stats.az * math.pi / 180.

            #Mxx term
            AJ[0][cnt1] = 0.5 * math.sin(2 * alpha) * gr[i*10+0].data[j]
            if(isoflag ==6):
              AJ[0][cnt2] = (1./6.) * (-1.0) * gr[i*10+4].data[j] - 0.5*math.cos(2. * alpha) * gr[i*10+2].data[j] + (1./3.) * (-1.0) * gr[i*10+8].data[j]
              AJ[0][cnt3] = (1./6.) * (-1.0) * gr[i*10+7].data[j] - 0.5*math.cos(2. * alpha) * (-1.0) * gr[i*10+5].data[j] + (1./3.) * (-1.0) * gr[i*10+9].data[j]
            else:
              AJ[0][cnt2] = (1./2.) *          gr[i*10+4].data[j] - 0.5*math.cos(2. * alpha) * gr[i*10+2].data[j]
              AJ[0][cnt3] = (1./2.) * (-1.0) * gr[i*10+7].data[j] - 0.5*math.cos(2. * alpha) * (-1.0) * gr[i*10+5].data[j]
               
            #Myy term
            AJ[1][cnt1] = -1.0 * 0.5 * math.sin(2. * alpha) * gr[i*10+0].data[j]
            if(isoflag ==6):
              AJ[1][cnt2] = (1./6.) *          gr[i*10+4].data[j] + 0.5*math.cos(2. * alpha) * gr[i*10+2].data[j] + (1./3.) * (-1.0) * gr[i*10+8].data[j]
              AJ[1][cnt3] = (1./6.) * (-1.0) * gr[i*10+7].data[j] + 0.5*math.cos(2. * alpha) * (-1.0) * gr[i*10+5].data[j] + (1./3.) * (-1.0) * gr[i*10+9].data[j]
            else:
              AJ[1][cnt2] = (1./2.)          * gr[i*10+4].data[j] + 0.5*math.cos(2. * alpha) * gr[i*10+2].data[j]
              AJ[1][cnt3] = (1./2.) * (-1.0) * gr[i*10+7].data[j] + 0.5*math.cos(2. * alpha) * (-1.0) * gr[i*10+5].data[j]
            
            #Mxy term
            AJ[2][cnt1] = (-1.0) * math.cos(2. * alpha) * gr[i*10+0].data[j]
            AJ[2][cnt2] = (-1.0) * math.sin(2. * alpha) * gr[i*10+2].data[j]
            AJ[2][cnt3] = (-1.0) * math.sin(2. * alpha) * (-1.0) * gr[i*10+5].data[j]
#           print("----   %d %d %g" % (i,j,gr[i*10+5].data[j]))

            #Mxz term
            AJ[3][cnt1] = (-1.0) * math.sin(alpha) * gr[i*10+1].data[j]
            AJ[3][cnt2] =          math.cos(alpha) * gr[i*10+3].data[j]
            AJ[3][cnt3] =          math.cos(alpha) * (-1.0) * gr[i*10+6].data[j]

            #Myz term*/
            AJ[4][cnt1] = math.cos(alpha) * gr[i*10+1].data[j]
            AJ[4][cnt2] = math.sin(alpha) * gr[i*10+3].data[j]
            AJ[4][cnt3] = math.sin(alpha) * (-1.0) * gr[i*10+6].data[j]

            #Mzz term*/
            if(isoflag==6):
              AJ[5][cnt1] = 0.0
              AJ[5][cnt2] = 1/3 * (-1.0) * gr[i*10+8].data[j] - 1/3 * (1.0) * gr[i*10+4].data[j]
              AJ[5][cnt3] = 1/3 * (-1.0) * gr[i*10+9].data[j] - 1/3 * (-1.0) * gr[i*10+7].data[j]

            # increment counters
            cnt1 += 1
            cnt2 += 1
            cnt3 += 1

    # Compute ATA
    for i in range(0,isoflag):
        for j in range(0,isoflag):
            for k in range(0,cnt3):
             AIV[i][j] += AJ[i][k]* AJ[j][k] * W[k];

    return AIV,B,AJ,W


def setStationWeigths(gr,st):

    # find closer distance
    mindist=1000000
    for i in range(len(st)):
        if(st[i].stats.dist <= mindist):
           mindist=st[i].stats.dist
   
    # set station weights
    for i in range(len(st)):
        st[i].stats.W = st[i].stats.dist/mindist
    for i in range(len(gr)):
        gr[i].stats.W = gr[i].stats.dist/mindist

    return (gr,st)


def realignPhaseSt(gr,st,args):

    #tim=int(int(args.len)/eval(args.DeltaInv))
    tim=int((float(args.len)+float(args.pre))/float(args.DeltaInv))
    if(args.zcor == "iso"):

      for i in range(len(st)):

          ph = int(st[i].stats.Zcor)

          if(ph >=1):
             st[i].data = np.roll(st[i].data, -ph)
             for l in range(1,ph):
                 st[i].data[-l] = 0

          elif(ph <=1):
             st[i].data = np.roll(st[i].data, ph)
             for l in range(0,ph):
                 st[i].data[l] = 0

          else:
             pass

    else:
      for i in range(len(st)):

          if  (gr[i].stats.channel == "VER") :
                ph = int(st[i].stats.Vcor)
          elif(gr[i].stats.channel == "RAD"):
                ph = int(st[i].stats.Rcor)
          else:
                ph = int(st[i].stats.Tcor)

          if(ph >=1):
             st[i].data = np.roll(st[i].data, -ph)
             for l in range(1,ph):
                 st[i].data[-l] = 0

          elif(ph <=1):
             st[i].data = np.roll(st[i].data, ph)
             for l in range(0,ph):
                 st[i].data[l] = 0

          else:
             pass
#     pass 

    return (gr,st)


def correlatePhaseShift(gr,st,args):

    # Correlation for Zcor for alla and for each component
    # T --> TSS,TDS         | S1 - u1,u2
    # R --> RSS,RDS,RDD     | S2 - u3,u4,u5
    # V --> ZSS,ZDS,ZDD     | S3 - u6,u7,u8
   
    # width for time cross correlation
    wid = int(float(st[0].stats.npts)/4)


    for i in range(len(st)/3):

          x = np.arange(16.).reshape(8,2)
          t = np.arange(4.).reshape(2,2)
          r = np.arange(6.).reshape(3,2)
          v = np.arange(6.).reshape(3,2)

          # Tangential
          a,b = xcorr(st[3*i+0], gr[10*i+0], wid)
          x[0][0]=abs(b)
          x[0][1]=a
          t[0][0]=abs(b)
          t[0][1]=a
          c,d = xcorr(st[3*i+0], gr[10*i+1], wid)
          x[1][0]=abs(b)
          x[1][1]=a
          t[1][0]=abs(b)
          t[1][1]=a
          

          # Radial
          a,b = xcorr(st[3*i+1], gr[10*i+2], wid)
          x[2][0]=abs(b)
          x[2][1]=a
          r[0][0]=abs(b)
          r[0][1]=a
          a,b = xcorr(st[3*i+1], gr[10*i+3], wid)
          x[3][0]=abs(b)
          x[3][1]=a
          r[1][0]=abs(b)
          r[1][1]=a
          a,b = xcorr(st[3*i+1], gr[10*i+4], wid)
          x[4][0]=abs(b)
          x[4][1]=a
          r[2][0]=abs(b)
          r[2][1]=a

          # Vertical
          a,b = xcorr(st[3*i+2], gr[10*i+5], wid)
          x[5][0]=abs(b)
          x[5][1]=a
          v[0][0]=abs(b)
          v[0][1]=a
          a,b = xcorr(st[3*i+2], gr[10*i+6], wid)
          x[6][0]=abs(b)
          x[6][1]=a
          v[1][0]=abs(b)
          v[1][1]=a
          a,b = xcorr(st[3*i+2], gr[10*i+7], wid)
          x[7][0]=abs(b)
          x[7][1]=a
          v[2][0]=abs(b)
          v[2][1]=a

          # sort for zcor
          X=np.array(sorted(sorted(x,key=lambda e:e[1]),key=lambda e:e[0]))
          T=np.array(sorted(sorted(t,key=lambda e:e[1]),key=lambda e:e[0]))
          R=np.array(sorted(sorted(r,key=lambda e:e[1]),key=lambda e:e[0]))
          V=np.array(sorted(sorted(v,key=lambda e:e[1]),key=lambda e:e[0]))
          Zco = X[-1][1]
          Tco = T[-1][1]
          Rco = R[-1][1]
          Vco = V[-1][1]

          # Update stats
          for l in range(0,3):
             st[3*i+l].stats.Zcor = Zco
             st[3*i+l].stats.Tcor = Tco
             st[3*i+l].stats.Rcor = Rco
             st[3*i+l].stats.Vcor = Vco
          for l in range(0,10):
             gr[10*i+l].stats.Zcor = Zco
             gr[10*i+l].stats.Tcor = Tco
             gr[10*i+l].stats.Rcor = Rco
             gr[10*i+l].stats.Vcor = Vco


    return (gr,st)

def trimStreams(gr,st,args):

    
    #tim=int(int(args.len)/eval(args.DeltaInv))
    tim=int((float(args.len)+float(args.pre))/float(args.DeltaInv))
    for i in range(len(gr)):
        gr[i].data=gr[i].data[0:tim]
    for i in range(len(st)):
        st[i].data=st[i].data[0:tim]

    return (gr,st)
 
