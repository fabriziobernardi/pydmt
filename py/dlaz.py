from math import fabs,sin,cos,sqrt
import numpy as np


def dlaz(lat_a,lon_a,lat_b,lon_b):

  d2r =  0.017453293    #degree-to-radians
  r2d = 57.295779515    #radians-to-degree
  deg2km = 6371 * d2r  #degree-to-kilometers
                        #this assumes the average Earth radius to be 6371 km

  #use geocentric latitude instead of geographic latitude
  #to calculate distance-aximuth-backazimuth

  #eccentricity e2:
  e2 = .0066943800      #Geodetic Reference System'80 GRS-80

  if lon_a > 360:
     lon_a -= 360
  if lon_b > 360:
     lon_b -= 360
  if lon_a < 0:
     lon_a += 360
  if lon_b < 0:
     lon_b += 360


  #convert from geograpic to geocentric latitude:
  lat_a_r = lat_a * d2r
  lat_b_r = lat_b * d2r
  lat_a_r_geoc = np.arctan2 ( (1-e2)*sin(lat_a_r)/cos(lat_a_r),1 )
  lat_b_r_geoc = np.arctan2 ( (1-e2)*sin(lat_b_r)/cos(lat_b_r),1 )
  lat_a = lat_a_r_geoc * r2d  #now geocentric latitude!
  lat_b = lat_b_r_geoc * r2d  #now geocentric latitude!

  #colatitude = 90 - latitude!
  colat_a = 90 - lat_a
  colat_b = 90 - lat_b

  colat_a_r = colat_a * d2r
  lon_a_r = lon_a * d2r
  colat_b_r = colat_b * d2r
  lon_b_r = lon_b * d2r

  #calculate the distance:
  #
  tmp = cos(colat_a_r)*cos(colat_b_r) + \
        sin(colat_a_r)*sin(colat_b_r)*cos(lon_b_r-lon_a_r)
  #acos_rad = np.arctan2(y, x)
  acos_rad = np.arctan2(sqrt(1-tmp*tmp),tmp)
  distance_r = acos_rad
  distance = distance_r * r2d
  distance_km = distance * deg2km

  # azimuth
  tmp = cos(colat_b_r)*sin(colat_a_r) - sin(colat_b_r)*cos(colat_a_r)*cos(lon_b_r-lon_a_r);
  tmp = tmp/sin(distance_r);
  acos_rad = np.arctan2(sqrt(1-tmp*tmp),tmp)
  az  = acos_rad * r2d
  tmp_2 = sin(colat_b_r)*sin(lon_b_r-lon_a_r)/sin(distance_r)
  if tmp_2 < 0:
      az = 360 - az

  # backazimut
  tmp = cos(colat_a_r)*sin(colat_b_r) - sin(colat_a_r)*cos(colat_b_r)*cos(lon_a_r-lon_b_r)
  tmp = tmp/sin(distance_r);
  acos_rad = np.arctan2(sqrt(1-tmp*tmp),tmp)
  baz = acos_rad * r2d
  tmp_2 = sin(colat_a_r)*sin(lon_a_r-lon_b_r)/sin(distance_r)
  if tmp_2 < 0:
      baz = 360 - baz


  return (distance,az,baz,distance_km)
