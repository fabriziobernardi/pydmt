#!/usr/bin/env python

import unittest
from obspy.imaging.beachball import PrincipalAxis

def stressRegime(P, N, T):
    """
    Stress regime assigment and maximum horizontal stress SH orientation

    From Zoback (1992), stress regime is defined from plunge (pl) of 
    the stress principal axes.  

    Table (from Zoback (1992))
    ------------------------------------------------------------------------------
    P/S1-axis           B/S2-axis       T/S3-axis       Regime  SH-azimuth
    ------------------------------------------------------------------------------
    pl >= 52                            pl <= 35        NF      azimuth of N
    40 <= pl <52                        pl <= 20        NS      azimuth of T+90
    pl < 40             pl >= 45        pl <= 20        SS      azimuth of T+90
    pl <= 20            pl >= 45        pl < 40         SS      azimuth of P
    pl <= 20                            40 <= pl < 52   TS      azimuth of P
    pl <= 35                            pl >= 52        TF      azimuth of P

    Orientations that do not satisfy the table above, resunts into "undefined" 
    stress regimes

    Stress Regimes:
    NF: Normal Fault
    NS: Normal Slip
    SS: Strike Slip
    TS: Thrust Slip
    TF: Thrust Fault
    
    References:
    Zoback, M. L., 1992, First and second order patterns of stress in the lithosphere:  
    the World Stress Map project:  Journal Geophysical Research, v. 97, p. 11703-11728

    :param P: PrincipalAxis(val=0, strike=0, dip=0)
    :param N: PrincipalAxis(val=0, strike=0, dip=0)
    :param T: PrincipalAxis(val=0, strike=0, dip=0)
    :return: stress regime and SH and azimuth of the maximum horizontal stress SH
    """

    # Sort axis to be sure to be consistent with T-N-P 
    [P, N, T] = sorted([P, N, T], key=lambda x: x.val)

    # Define regimes
    regimes = {1: "Normal Fault", 
               2: "Normal Slip", 
               3: "Strike Slip", 
               4: "Thrust Slip", 
               5: "Thrust Fault", 
               6 : "Undefined"} 

    # Assign regimes
    if(P.dip >= 52 and T.dip <=35):                      # Normal Fault
      return (regimes[1], N.strike)
    elif(P.dip >= 40 and P.dip < 52 and T.dip <= 20):     # Normal Slip
      return (regimes[2], T.strike+90)
    elif(P.dip < 40 and N.dip >= 45 and T.dip <= 20):     # Strike Slip
      return (regimes[3], T.strike+90)
    elif(P.dip <= 20 and N.dip >= 45 and T.dip < 40):     # Strike Slip
      return (regimes[3], P.strike)
    elif(P.dip <= 20 and T.dip <= 40 and T.dip < 52):     # Thrust Slip 
      return (regimes[4], P.strike)
    elif(P.dip <= 35 and T.dip >= 52):                   # Thrust Fault 
      return (regimes[5], P.strike)
    else:                                              # Undefined
      return (regimes[6], 999)                        

class TestStressRegime(unittest.TestCase):
    def test_stressRegime(self):
        """
        """
        # Test Normal Faut
        A1 = PrincipalAxis(5.196e21, 241, 2)
        A2 = PrincipalAxis(1.696e20, 331, 2)
        A3 = PrincipalAxis(-5.366e21, 105, 87)
        (regime, sh) = stressRegime(A1, A2, A3)
        self.assertEqual(regime, "Normal Fault")
        self.assertEqual(sh, 331)

        # Test Strike Slip
        A1 = PrincipalAxis(2.564, 267, 10)
        A2 = PrincipalAxis(-0.171, 56, 78)
        A3 = PrincipalAxis(-2.393, 176, 6)
        (regime, sh) = stressRegime(A1, A2, A3)
        self.assertEqual(regime, "Strike Slip")
        self.assertEqual(sh, 357)

        # Test Thrust Fault
        A1 = PrincipalAxis(val=8.94, dip=75, strike=283)
        A2 = PrincipalAxis(val=1.26, dip=2, strike=19)
        A3 = PrincipalAxis(val=-10.19, dip=15, strike=110)
        (regime, sh) = stressRegime(A1, A2, A3)
        self.assertEqual(regime, "Thrust Fault")
        self.assertEqual(sh, 110)

def suite():
    return unittest.makeSuite(TestStressRegime, 'test')

if __name__ == '__main__':
    unittest.main()
#   import doctest
#   doctest.testmod()

