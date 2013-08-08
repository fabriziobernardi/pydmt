#!/usr/bin/env python -W ignore::DeprecationWarning
# encoding: utf-8

from obspy.core import Stream,Trace
from FortranFile import FortranFile
from math import acos,sin,exp
from processData import getStationslist
import numpy,sys



def upGrStats(gr,st):
    
    for j in range(len(gr)):
        Dgr=float(gr[j].stats.dist)

        for i in range(len(st)):
            Dga=float(st[i].stats.dist)
            difDist=abs(Dgr-Dga)

     
            if (difDist<=0.0100):
                # station metadata
                gr[j].stats.station   = st[i].stats.station
                gr[j].stats.dist      = st[i].stats.dist
                gr[j].stats.az        = st[i].stats.az
                gr[j].stats.baz       = st[i].stats.baz
                gr[j].stats.gcarc     = st[i].stats.gcarc
                gr[j].stats.stlo      = st[i].stats.stlo
                gr[j].stats.stla      = st[i].stats.stla
                gr[j].stats.evlo      = st[i].stats.evlo
                gr[j].stats.evla      = st[i].stats.evla
                # begin time
                gr[j].stats.starttime = st[i].stats.starttime
                # fit
                gr[j].stats.Zcor      = st[i].stats.Zcor
                gr[j].stats.Tcor      = st[i].stats.Tcor
                gr[j].stats.Rcor      = st[i].stats.Rcor
                gr[j].stats.Vcor      = st[i].stats.Vcor
                break

    return gr

def reorderGreen(gr,Sta):

    temp = Stream()
    out  = Stream()

    for i in range(len(Sta)):
        for n in range(len(gr)):
            if(gr[n].stats.station == Sta[i]):
               temp.append(gr[n])

        out.append(temp[7])
        out.append(temp[4])
        out.append(temp[6])
        out.append(temp[3])
        out.append(temp[1])
        out.append(temp[5])
        out.append(temp[2])
        out.append(temp[0])
        out.append(temp[8])
        out.append(temp[9])
        temp = Stream()

    return out


def aquireGreens(greenFile,args):


   # define strem to load the greens in time
   green=Stream()

   # displacement o velocity
   ivel = args.dva

   # max number of stations (i.e.: distance) allowed. 
   # when modify this value, modify also the same value
   # range, vred and t0 into FKRPROG and recomplie
   Max_nr_dists = 100
   # dim(n2) (4097)

   # read twice the Green.1 file. One to access integers and one for reals 
   F = FortranFile(greenFile)
   I = FortranFile(greenFile)

   # First Line of header. This line is an mix float-int array 
   h_1f   = F.readReals()
   h_1i   = I.readInts()
   alpha  = h_1f[0]
   depth  = h_1f[1]
   fl     = h_1f[2]
   fu     = h_1f[3]
   dt     = h_1f[4]
   n1     = h_1i[5]
   n2     = h_1i[6]
   df     = h_1f[7]
   nyq    = h_1i[8]
   nrange = h_1i[9]
   nskip  = h_1i[10]

   # Second Line of header. This line is a 10 integer vector (isrc in code) 
   # If h_2i[j] == 1 then compute green for this source element
   # else do not.
   h_2f   = F.readReals()  # Usealess for this line but neaded 
                          # to place pointers at the end of this array
                          # for next header line
   h_2i   = I.readInts()
   isrc   = h_2i

   # Third Line of header. Mix float-int array. All float and one int
   # this heder inludes 6 float arry with 70 elemnts and one integer(nmax)
   # after the first 4 vectors.
   # Name of vectors and constant:
   # d(70),a(70),b(70),rho(70),nmax(1),qa(70),qb(70)
   h_3f   = F.readReals()
   h_3i   = I.readInts()
   beg    = 0
   inc    = 70
   d      = h_3f[beg:beg+inc]   
   beg    = beg+inc
   a      = h_3f[beg:beg+inc]
   beg    = beg+inc
   b      = h_3f[beg:beg+inc]
   beg    = beg+inc
   rho    = h_3f[beg:beg+inc]
   beg    = beg+inc
   mmax   = h_3i[beg:beg+1] 
   beg    = beg+1
   qa     = h_3f[beg:beg+inc]
   beg    = beg+inc
   qb     = h_3f[beg:beg+inc]
   
   # 4th line: range of distances, vred and t0
   # this refers to stations lines into earth model
   # maximum station allowed in 100 and is hardcoded into
   # FKRPROG.f To extend the number of stationallowed, modify the code
   # recomplie and modify parameters at the beginn of this procedure
   # Vectors of float
   h_4f   = F.readReals()
   h_4i   = I.readInts()
   beg    = 0
   inc    = Max_nr_dists
   Range  = h_4f[beg:beg+inc]
   beg    = beg+inc
   vred   = h_4f[beg:beg+inc]
   beg    = beg+inc
   t0     = h_4f[beg:beg+inc]

   # 5th line --> loop over distances to access cmplx spectral values
   # 3 loops:
   # 1 i: 1->n2 (n2=npts/2) -  omega(float),nkk(int)
   #   2  j: 1->nrange
   #      3  k: 1->10
   #           read(aa),(bb) --> gg(j,i,k)=cmplx(sngl(aa),sngl(bb))
   omega = [0.0 for x in range(0,n2+0)]
   freq  = [0.0 for x in range(0,n2+0)]
   # (complex spectral matrix) 
   # z --- k: foundamental mts
   # y --- i: (npts/2). nyq = npts/2 + 1
   # x --- j: distsances
   gg    = [ [ [ 0.0 for z in range(0,10+0)] \
                   for y in range(0,n2+0)]     \
                   for x in range(0,nrange+0)] 

   for i in range(0,n2+0):
       h_5f   = F.readReals()
       omega[i] = h_5f[0]
       
       for j in range(0,nrange+0):
           
           for k in range(0,10+0):
               
               FF = F.readDouble()
               aa = complex(FF[0],FF[1])
               gg[j][i][k] = aa
 
       
   # here we have the resulting matrix gg with complex spectra values
   # for each distance, each mtfound, take vector of complex (generate data(j), 
   # make conjg, add complex for nyq, spectra->time, damping, 
   # integration if required(disp)
   pi   = acos(-1.0)
   twpi = 2.*pi
   n    = 2 * (nyq-1)
   nm     = n2
   npoint = n
   rep    = 'n'
   tau = dt
   fmax = fu
   inst = 0

   for j in range(0,nrange+0):
       t0x = (Range[j])/(vred[j])
       yr  = 0.0

       for k in range(0,10+0):

           if (isrc[k] ==1):

              # inizialize data
              data  = numpy.array(numpy.zeros(n2+0+n2),dtype=complex)
#             data  = [0.0 for x in range(0,n2+0+n2)]
              for i in range(0,n2+0):
                  # arrange data
                  data[i] = gg[j][i][k]
                  # arrang frequency
                  freq = i*df
                  if(freq < df):
                     freq = 0.01*df
                  om = twpi * freq
              for i in range(n2+0,nyq):
                  data[i] = complex(0.0,0.0)
 
              # conjug
              for i in range(1,n2+0):
                  data[n+0-i] = data[i].conjugate()

              data[0]     = complex(0.0,0.0)
              data[nyq-1] = complex(0.0,0.0)


              # From spectraToTime
              data = four1(data,n,1,dt,df)

              # Apply damping factor
              fac = exp(alpha*t0x)
              dfac = exp(alpha*dt)
              for i in range(len(data)):
                data[i]=(data[i])*fac
                fac = fac * dfac
              
              # velocity to displavement if required
              if(ivel=='1'):
                 data=velTodisData(data,dt)

              # put data into trace
              prel = 0
              prel = int(eval(args.pre) / eval(args.delta))
              length=numpy.arange(len(data)+prel)*0.0
              t=Trace(length)
              for i in range(len(data)):
                  t.data[i+prel]=data[i]
              t.stats['delta']   = dt
              t.stats['dist']    = Range[j]
              name = 'GREEN_' + str(t.stats['dist'])
              t.stats['station'] = name


              # update stats and apply -1 for tss,xds,zss,zdd
              if k==7:
                 t.stats['channel'] = 'tss'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(-1)
              if k==4:
                 t.stats['channel'] = 'tds'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(+1)
              if k==6:
                 t.stats['channel'] = 'xss'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(+1)
              if k==3:
                 t.stats['channel'] = 'xds'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(-1)
              if k==1:
                 t.stats['channel'] = 'xdd'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(+1)
              if k==5:
                 t.stats['channel'] = 'zss'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(-1)
              if k==2:
                 t.stats['channel'] = 'zds'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(+1)
              if k==0:
                 t.stats['channel'] = 'zdd'  
                 for i in range(len(data)):
                     t.data[i]             = t.data[i]*(-1)
              if k==8:
                 t.stats['channel'] = 'ex1'  
#                t.data[i]                 = t.data[i]
#                if(args.iso=='1'):
#                  for i in range(len(data)):
#                    t.data[i]             = 0.0
#                else:
#                  for i in range(len(data)):
#                    t.data[i]             = 0.0
              if k==9:
                 t.stats['channel'] = 'ex2'  
#                t.data[i]                 = t.data[i]
#                if(args.iso=='1'):
#                  for i in range(len(data)):
#                    t.data[i]             = 0.0 
#                else:
#                  for i in range(len(data)):
#                    t.data[i]             = 0.0



              green.append(t)


   return green

def velTodisData(data,dt):

    totint = 0.0
    for i in range(len(data)-1):
        prtint = 0.5*dt*(data[i]+data[i+1])
        totint = totint + prtint
        data[i]=totint
    ref  = data[0]
    for i in range(len(data)-1):
        data[i]=data[i]-ref 

    return data

def swap(t1,t2):
    return t2,t1

def four1(spec, nn, isign, dt, df):

    data=[-1]
    for i in range(len(spec)):
        data.append(spec[i].real)
        data.append(spec[i].imag)

    n=2 * nn
    j=1
    for i in range(1,n-1,2):
        if (j>i):
           data[j],data[i]     = swap(data[j],data[i])
           data[j+1],data[i+1] = swap(data[j+1],data[i+1])
        m=n/2
        while(m>=2 and j>m):
           j = j-m
           m = m/2
        j += m

    mmax=2

    while (n>mmax):
        istep=mmax * 2
        theta=isign*(6.28318530717959/mmax)
        wtemp=sin(0.5*theta)
        wpr = -2.0*wtemp*wtemp
        wpi=sin(theta)
        wr=1.0
        wi=0.0
        for m in range(1,mmax-0,2):
            for i in range(m,n-1,istep):

                j=i+mmax
                tempr     = wr*data[j]  -wi*data[j+1]
                tempi     = wr*data[j+1]+wi*data[j]
                data[j]   = data[i]  -tempr
                data[j+1] = data[i+1]-tempi
                data[i]   = data[i]  +tempr
                data[i+1] = data[i+1]+tempi

            wtemp=wr
            wr = wr*wpr-wi*wpi+wr
            wi = wi*wpr+wtemp*wpi+wi

        mmax = istep

    for i in range(n):
        data[i] = data[i] * df

    # remove first elemnt of dtata array
    del data[0]
    out = []
    for i in range(0,n,2):
        out.append(data[i])

    
    return out

