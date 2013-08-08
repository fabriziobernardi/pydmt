import numpy as np
from processData import getSTationslist
from obspy.core import Stream, read
from scipy.optimize import curve_fit
from pylab import plot,show
from obspy.signal.cross_correlation import xcorr
import sys

def purgeStream(st,l):

    # Parameters: st: Stream, l: linst with station number [0 to nr stations]

    new=Stream()

    ls_Sta = getSTationslist(st)
    for i in range(len(l)):
        a = l[i]-i
        ls_Sta.pop(a)

    for i in range(len(st)):
        sta = st[i].stats.station
        rem = [y for y in ls_Sta if sta == y]
        if (len(rem)==1):
           new.append(st[i])
    
    return new


def cleanNoise(ob,args):

    List = []
    Max_noise = float(args.noise)
    if Max_noise >= 0:

      for i in range(len(ob)/3):

        # set temporary data vectors
        r = np.zeros(len(ob[i*3+0])/3)
        t = np.zeros(len(ob[i*3+1])/3)
        z = np.zeros(len(ob[i*3+2])/3)

        # set noise (random uniform) for reference
        n = np.random.uniform(-1,1,len(ob[i*3+0]))
        w = len(n)/4

        # Allocate data values
        r = ob[i*3+0].data/max(ob[i*3+0].data)
        t = ob[i*3+1].data/max(ob[i*3+1].data)
        z = ob[i*3+2].data/max(ob[i*3+2].data)

        # 3 comp reference values
        N_r = xcorr(n,r,w) 
        N_t = xcorr(n,t,w) 
        N_z = xcorr(n,z,w) 
        # find for
        R_t = xcorr(r,t,w)
        R_z = xcorr(r,z,w)
        T_z = xcorr(t,z,w)

        M_Rp = (R_t[0] + R_z[0] + T_z[0])/3.
        M_Np = (N_t[0] + N_z[0] + N_z[0])/3.
        M_Rv = (R_t[1] + R_z[1] + T_z[1])/3.
        M_Nv = (N_t[1] + N_z[1] + N_z[1])/3.

        
        if (abs(M_Rv) < Max_noise):
           List.append(i)

      if (len(List) > 0):

         listToRemove = ' '.join(['%d:%s' % (List[n],ob[List[n]*3].stats.station) for n in xrange(len(List))])
         print "\n\nStation to clean because of noise: ",listToRemove,"\n\n"

         # Purge station from observed
         ob  = purgeStream(ob,List)
         if(len(ob)==0):
            print "No data left. No solution found. Exit"
            sys.exit()
 
    return (ob, List)

def cleanSpike(ob, sy, args):

    # set empty list of stations to remove because of spikes
    cleanList = []

    # set the length of the x and y array
    # X is distance in km from epicenter
    # Y is the maximal amplitude, for obs and for syn
    # M is the mean of the Y value
    X       = np.zeros(len(ob)/3)
    Y_obs_r = np.zeros(len(ob)/3)
    Y_obs_t = np.zeros(len(ob)/3)
    Y_obs_z = np.zeros(len(ob)/3)
    Y_syn_r = np.zeros(len(ob)/3)
    Y_syn_t = np.zeros(len(ob)/3)
    Y_syn_z = np.zeros(len(ob)/3)
    M_obs   = np.zeros(len(ob)/3)
    M_syn   = np.zeros(len(ob)/3)

    # now put max amplitudes into np.arrays
    for i in range(len(ob)/3):
        u = i*2
        X[i]       = float(ob[i*3+0].stats.dist) 
        Y_obs_r[i] = max(np.absolute(ob[i*3+0].data))
        Y_obs_t[i] = max(np.absolute(ob[i*3+1].data))
        Y_obs_z[i] = max(np.absolute(ob[i*3+2].data))
        Y_syn_r[i] = max(np.absolute(sy[i*3+0].data))
        Y_syn_t[i] = max(np.absolute(sy[i*3+1].data))
        Y_syn_z[i] = max(np.absolute(sy[i*3+2].data))
        M_obs[i]   = (Y_obs_r[i] + Y_obs_t[i] + Y_obs_z[i]) / 3
        M_syn[i]   = (Y_syn_r[i] + Y_syn_t[i] + Y_syn_z[i]) / 3

    # now normalize Y and M avlues to 1
    maxM_obs = max(M_obs)
    maxM_syn = max(M_syn)
    M_obs    = M_obs/maxM_obs 
    M_syn    = M_syn/maxM_syn 
  
 
    # now try fit 
    print "\n\n"
    obs_opt, obs_cov = curve_fit(func, X, M_obs)
    print "OBS ",obs_opt[0], obs_opt[1]
#   print "OBS ",obs_cov
    syn_opt, syn_cov = curve_fit(func, X,M_syn)
    print "SYN ",syn_opt[0],syn_opt[1]
#   print "SYN ",syn_cov
 
    # Now put condition of b: if b(obs)<0 and b(syn)>0, then search for station to clean

#   print "ZZZZZZZZ A Obs Syn - B Obs Syn",args.title,obs_opt[0],syn_opt[0],"  -  ",obs_opt[1] ,syn_opt[1] 
#   if   (obs_opt[1] >= 0):
#      pass

#   else:
    if (syn_opt[0] > 0 and syn_opt[1] > 0):
       if (obs_opt[1] < 0 or obs_opt[0] <0):
         # compute residuals of max amplitudes between observed and theoretical 
         res = getRes(X,M_obs,syn_opt[0],syn_opt[1],ob)
         # compute mead and std to remove outliers

         # to reject outlayers m*std thershold
         cleanList = reject_outliers(res,2)

#     yo = func(X,obs_opt[0],obs_opt[1])
#     ys = func(X,syn_opt[0],syn_opt[1])
#     plot(X,yo,color='k')
#     plot(X,ys,color='r')
#     show()
#     sys.exit()
 
    return cleanList
    
    
def getRes(x,m,a,b,ob):

    # Compute residuals of max amplitudes between observed and theoretical 
    # y: vector of residual
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = m[i] - func(x[i],a,b)
#       print("%d  %10.6f %10.6f %10.6f %s" % (i,y[i],m[i],func(x[i],a,b),ob[i*3+0].stats.station))
        
    return y

def reject_outliers(data, m):

    y = np.zeros(len(data))
    y = []

    mean_data = np.mean(data) 
    std_data  = np.std(data) 

    for i in range(len(data)):
       if (abs(data[i] - mean_data) > m * std_data):
          y.append(i)
#         y[i] = 1

    return y


def reject_outliers2(data, m=2):
    return data[abs(data - np.mean(data)) > m * np.std(data)]

def func(x,a,b):
    return a/(x)+b
