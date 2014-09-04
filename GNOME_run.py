# -*- coding: utf-8 -*-
"""
Created on Tue Nov 05 19:48:40 2013

@author: xianl_000
"""
import os 
utm_x=666393 ; utm_y=3076780
coor=[[utm_x,utm_y],[utm_x-10,utm_y],[utm_x-20,utm_y],[utm_x+10,utm_y],[utm_x+20,utm_y],[utm_x,utm_y-10],
      [utm_x,utm_y-20],[utm_x,utm_y+10],[utm_x,utm_y+20],[utm_x-10,utm_y+10],[utm_x-10,utm_y-10],[utm_x+10,utm_y+10],[utm_x+10,utm_y-10]]

import math

def utmToLatLng(zone, easting, northing, northernHemisphere=True):
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

    return [latitude,longitude]


from datetime import datetime, timedelta
from urllib2 import HTTPError
import numpy
np = numpy
from gnome import scripting
from gnome.basic_types import datetime_value_2d
from gnome.utilities.remote_data import get_datafile
from gnome.model import Model
from gnome.map import MapFromBNA
from gnome.environment import Wind
from gnome.spill import point_line_release_spill
from gnome.movers import RandomMover, WindMover, GridCurrentMover
from gnome.outputters import Renderer
from gnome.outputters import NetCDFOutput


# define base directory

base_dir = os.path.dirname(__file__)

def make_model(images_dir=os.path.join(base_dir,"images")):
    print "initializing the model"

    start_time = datetime(2013, 7, 23, 0)
    model = Model(start_time = start_time,
                              duration = timedelta(hours=47),	# n+1 of data in file
                              time_step = 900, # 4 hr in seconds
                              uncertain = False,
                              )
    
    mapfile = os.path.join(base_dir, './coast.bna')
    print "adding the map"
    gnome_map = MapFromBNA(mapfile, refloat_halflife=6)  # hours
    
    print "adding renderer" 
    model.outputters += Renderer(mapfile, images_dir, size=(1800, 1600))

    print "adding a wind mover from a time-series"
    ## this is wind
    wind_file=get_datafile(os.path.join(base_dir, 'wind.WND'))
    wind = Wind(filename=wind_file)
    w_mover = WindMover(wind)
    model.movers += w_mover
    
    print "adding a current mover:"
    ## this is currents
    curr_file = get_datafile(os.path.join(base_dir, 'current.txt'))
    model.movers += GridCurrentMover(curr_file)

    ##
    ## Add some spills (sources of elements)
    ##
    print "adding 13 points in a cluster that has some small initial separation as the source of spill"
    
    for i in range(len(coor)):
        
        aaa=utmToLatLng(14,coor[i][0],coor[i][1],northernHemisphere=True)
        model.spills += point_line_release_spill(num_elements=1,
                                                start_position = (aaa[1],aaa[0], 0.0),
                                                release_time = start_time,
                                                )

    print "adding netcdf output"
    netcdf_output_file = os.path.join(base_dir,'GNOME_output.nc')
    scripting.remove_netcdf(netcdf_output_file)
    model.outputters += NetCDFOutput(netcdf_output_file, which_data='all')

    return model


if __name__ == "__main__":
    """ if called on its own -- run it """

    scripting.make_images_dir()
    model = make_model()
    model.full_run(log=True)
