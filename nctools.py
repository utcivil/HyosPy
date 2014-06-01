# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 19:53:18 2013

@author: xianl_000
"""

#!/usr/bin/env python


from netCDF4 import num2date, date2num, date2index
import glob, os
import math


def show_tbounds(t):
    
    print 'Start date: ', num2date(t[0],t.units)
    try:
        print 'End date: ', num2date(t[-1],t.units)
    except IndexError:
        print num2date(t[0],t.units)
        
def get_tindex(t,start_date,end_date,stride=None):
        
    tindex = []
    tindex.append(date2index(start_date,t,select='before'))
    tindex.append(date2index(end_date,t,select='after') + 1)
    if stride is None:
        tindex.append(1)
    else:
        tindex.append(stride)
    return tindex


def adjust_time(t,t_units):
    #GNOME can't handle pre-1970 date units
    
    dtime = num2date(t,t_units)
    new_units = 'days since 1980-1-1 00:00:00'
    new_time = date2num(dtime,units=new_units)
    
    return new_time,new_units


def make_filelist_for_GNOME(file_dir,file_match='*.*',outfilename='filelist.txt'):
    #used to load multiple files as one mover in GNOME
    flist = glob.glob(os.path.join(file_dir,file_match))
    f = open(os.path.join(file_dir,outfilename),'w')
    f.write('NetCDF Files\n')
    f.write('\n'.join(['[FILE] ' + os.path.split(file)[-1] for file in flist]))
    f.close()


def utmToLatLng(zone, easting, northing, northernHemisphere=True):
    # Convert UTM coordinates to lat/lon
    if not northernHemisphere:
        northing = 10000000 - northing


    a = 6378137
    e = 0.081819191
    e1sq = 0.006739497
    k0 = 0.9996


    arc = northing / k0
    mu = arc / (a * (1 - math.pow(e, 2) / 4.0 - 3 * math.pow(e, 4) / 64.0 - 5 * math.pow(e, 6) / 256.0))


    ei = (1 - math.pow((1 - e * e), (1 / 2.0))) / (1 + math.pow((1 - e * e), (1 / 2.0)))


    ca = 3 * ei / 2 - 27 * math.pow(ei, 3) / 32.0


    cb = 21 * math.pow(ei, 2) / 16 - 55 * math.pow(ei, 4) / 32
    cc = 151 * math.pow(ei, 3) / 96
    cd = 1097 * math.pow(ei, 4) / 512
    phi1 = mu + ca * math.sin(2 * mu) + cb * math.sin(4 * mu) + cc * math.sin(6 * mu) + cd * math.sin(8 * mu)


    n0 = a / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (1 / 2.0))


    r0 = a * (1 - e * e) / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (3 / 2.0))
    fact1 = n0 * math.tan(phi1) / r0


    _a1 = 500000 - easting
    dd0 = _a1 / (n0 * k0)
    fact2 = dd0 * dd0 / 2


    t0 = math.pow(math.tan(phi1), 2)
    Q0 = e1sq * math.pow(math.cos(phi1), 2)
    fact3 = (5 + 3 * t0 + 10 * Q0 - 4 * Q0 * Q0 - 9 * e1sq) * math.pow(dd0, 4) / 24


    fact4 = (61 + 90 * t0 + 298 * Q0 + 45 * t0 * t0 - 252 * e1sq - 3 * Q0 * Q0) * math.pow(dd0, 6) / 720


    lof1 = _a1 / (n0 * k0)
    lof2 = (1 + 2 * t0 + Q0) * math.pow(dd0, 3) / 6.0
    lof3 = (5 - 2 * Q0 + 28 * t0 - 3 * math.pow(Q0, 2) + 8 * e1sq + 24 * math.pow(t0, 2)) * math.pow(dd0, 5) / 120
    _a2 = (lof1 - lof2 + lof3) / math.cos(phi1)
    _a3 = _a2 * 180 / math.pi


    latitude = 180 * (phi1 - fact1 * (fact2 + fact3 + fact4)) / math.pi


    if not northernHemisphere:
        latitude = -latitude


    longitude = ((zone > 0) and (6 * zone - 183.0) or 3.0) - _a3


    return latitude, longitude
