# -*- coding: utf-8 -*-
"""
Oil spill wrapper tools
"""

import os
import subprocess
import shutil
import utools
import nctools 
import datetime as dt
import numpy as np
from netCDF4 import date2num
from netCDF4 import Dataset
import pyselfe_v1

base_dir = os.path.dirname(__file__)

def rk4(number,dt,h,hr,x,y):
    
    '''
    number--set up how many SELFE's outputs to use;
    dt,h--RK4 time step;
    hr--set up how many hours of hydrodynamic data to be used for oil spill modeling
    x,y--UTM coordinates for oil spill location
    '''
    
    a=open(base_dir+'/multiple_models_with_multiple_points_v1.py','r')
    row=[]
    for i in a.readlines():
        row.append(i)
        
    
    row[126]='number='+str(number)+'\n'
    row[131]='dt='+str(dt)+' ; h='+str(h)+'\n'
    row[138]='    '+'[t, t_iter, eta, dp, mdata] = selfe.read_time_series('+"'hvel.64'"+', nfiles='+str(hr)+', sfile=169, datadir=selfedir[nmodels])'+'\n'
    row[141]='    '+'utm_x='+str(x)+' ; utm_y='+str(y)+'\n'
    
    b=open(base_dir+'/rk4.py','w')
    for l in row:
        g=''.join([str(j) for j in l])
        b.write(g)
    b.close 
    
    subprocess.Popen(['python','rk4.py'])
    
    
def mul_GNOME_inputs(number,yr,month,day,hr,period):
    
    '''
    Transform multiple SELFE combined binary outputs into GNOME's inputs in NetCDF format;
    number is SELFE number;
    Oil spill starting time: yr,month,day,hr; month is based on 30 days
    Oil spill simulation time period: period (hours)
    
    '''
        
    data_file = os.path.join(base_dir,'169_hvel.nc')
    
    var_map = { 'longitude':'lon', \
                'latitude':'lat', \
                'time':'', \
                'u_velocity':'u', \
                'v_velocity':'v', \
                'nodes_surrounding_ele':'ele',\
                'eles_surrounding_ele':'',\
              }  


    txselfe = utools.ugrid(data_file)

    print 'Downloading data dimensions'

    x = txselfe.Dataset.variables['x'][:]
    y = txselfe.Dataset.variables['y'][:]
    lon = np.ones_like(x); lat = np.ones_like(x)
    for ii in range(len(x)):
        lat[ii], lon[ii] = nctools.utmToLatLng(14,x[ii],y[ii])
    txselfe.data['lon'] = lon
    txselfe.data['lat'] = lat
    txselfe.atts['lon'] = {'long_name': 'longitude'}
    txselfe.atts['lat'] = {'long_name': 'latitude'}
    
    
    # get grid topo variables (nbe, nv)
    print 'Downloading grid topo variables'
    try:
        txselfe.get_grid_topo(var_map)
    except KeyError: #model output on server doesn't have nbe
        txselfe.build_face_face_connectivity()
        
    # GNOME requires boundary info -- this file can be read form data_files directory
    # if saved or generated
    print 'Loading/generating boundary segments'
    bndry_file = os.path.join(base_dir, 'txselfe.bry')
    try:
        txselfe.read_bndry_file(bndry_file)
    except IOError:
        txselfe.write_bndry_file('txselfe',bndry_file)
        txselfe.read_bndry_file(bndry_file)
    txselfe.data['nbe'] = txselfe.data['nbe']
    txselfe.data['nv'] = txselfe.data['nv']
    # GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
    txselfe.atts['nbe']['order'] = 'ccw'
    
    
    # get the SELFE data
    print 'Loading u/v'

    # specify SELFE output path
    for n in range(number):

        shutil.copy(base_dir+'/txselfe.bry',base_dir+'/'+str(n+1)+'/GNOME')
    
        selfe = pyselfe_v1.Dataset(base_dir+'/'+str(n+1)+'/outputs/169_hvel.64')

        for j in range(period):
    
            model_time = dt.datetime(yr,month+j/720,day+(j/24)%30,hr+j%24,0,0) 
            t_units = 'hours since 2012-01-01 00:00:00'
            txselfe.data['time'] = [date2num(model_time,t_units),]
            txselfe.atts['time'] = {'units':t_units}
        
            num_nodes = len(txselfe.data['lon'])
            txselfe.data['u'] = np.ones([1,num_nodes],)
            txselfe.data['v'] = np.ones([1,num_nodes],)
            txselfe.atts['u'] = {'long_name':'eastward_velocity','units':'m/s'}
            txselfe.atts['v'] = {'long_name':'northward_velocity','units':'m/s'}
    
            [t, t_iter, eta, dp, mdata] = selfe.read_time_series('hvel.64', nfiles=1, sfile=169+j, datadir=base_dir+'/'+str(n+1)+'/outputs/')    
    
            for i in range(num_nodes):
                txselfe.data['u'][0][i]=mdata[0,i,5,0]
                txselfe.data['v'][0][i]=mdata[0,i,5,1]
 
            print 'Writing to GNOME file'
            txselfe.write_unstruc_grid(os.path.join(base_dir, str(n+1)+'/GNOME/'+str(169+j)+'_hvel.nc'))


def run_mul_GNOME(number,x,y,yr,month,day,period,dt):
    
    '''
    number--set up how many GNOME to run
    Run multiple GNOME;
    x,y--UTM coordinates for oil spill location;
    yr,month,day---oil spill start time;
    period---oil spill simulation duration;
    dt--oil spill time step in second
    '''
    
    a=open(base_dir+'/GNOME_run.py','r')
    row=[]
    for i in a.readlines():
        row.append(i)
    
    row[7]='utm_x='+str(x)+' ; utm_y='+str(y)+'\n'
    row[88]='    '+'start_time = datetime('+str(yr)+','+str(month)+','+str(day)+',0)'+'\n'
    row[90]=' '*30+'duration = timedelta(hours='+str(period-1)+'),'+'\n'
    row[91]=' '*30+'time_step ='+str(dt)+','+'\n'
    
    for j in range(number):
        
        h=open(base_dir+'/'+str(j+1)+'/GNOME/GNOME_run.py','w')
        for x in row:
            g=''.join([str(n) for n in x])
            h.write(g)
        h.close 
        
        subprocess.Popen(['python','GNOME_run.py'], cwd=base_dir+'/'+str(j+1)+'/GNOME')
        

def GNOME_GM_visualization(number):
    
    '''
    Visualize different GNOME oil spill trajectories on 2D Google Map GIS;
    number --- set up how many set of GNOME outputs to visualize
    '''
        
    a=[None]*number ; locations=[]

    print "Generating different oil spill tracks (GNOME) on Google Map"
    
    for i in range(number):
    
        a[i]=Dataset(base_dir+'/'+str(i+1)+'/GNOME/GNOME_output.nc','r')
    
        for j in range(len(a[i].variables[u'longitude'][:])):
        
            locations.append([a[i].variables[u'latitude'][j],a[i].variables[u'longitude'][j]])
        
    file1 = open('./javascript.txt','r')     # open the javascript
    row=[]
    for s in file1.readlines():        
        row.append(s)

    row[11]='    '+'var'+' '+'locations'+'='+str(locations)   
    
    h=open('./GNOME_Google_map.html','w')
    for l in row:
        g=''.join([str(j) for j in l])
        h.write(g+'\n')
    h.close        
        

def GNOME_GE_animation(number,np,yr,month,day):
    
    '''
    Animate different GNOME oil spill trajectories on 3D Google Earth GIS;
    number --- set up how many set of GNOME outputs to visualize;
    np is the total particle number; 
    yr,month,day---oil spill start time
    '''
    
    a=[None]*number
    
    for i in range(number):
        a[i]=Dataset(base_dir+'/'+str(i+1)+'/GNOME/GNOME_output.nc','r')
        
    file1 = open('./GE_animation.txt','r')     
    row=[]
    for s in file1.readlines():     
        row.append(s)
        
    nt=len(a[0].variables[u'time'][:])   #  nt is the total time steps

    print "Generating multiple oil spill tracks (GNOME) on 3D Google Earth"
    
    for n in range(nt):   
    
        for m in range(number):                            
         
            for j in range(np*n,np*(n+1)):                       
            
                row.append('  <Placemark>'+'\n'+'    <TimeStamp>'+'\n'+'      <when>'+str(yr)+'-'+str('%02d' % month)+'-'+str('%02d' % (day+n/96))+'T'+ str("%02d" % ((n/4)%24))+':' +str("%02d" %((n*15)%60))+':00Z</when>'+'\n'+'    </TimeStamp>'+'\n'+'    <styleUrl>#'+str(m+1)+'</styleUrl>'+'\n'+'    <Point>'+'\n'+'      <coordinates>'+str(a[m].variables[u'longitude'][j])+','+str(a[m].variables[u'latitude'][j])+'</coordinates>'+'\n'+'    </Point>'+'\n'+'  </Placemark>'+'\n')
                
    row.append('</Document>'+'\n'+'</kml>')    
      
    h=open('./GNOME_GE.kml','w')
    for l in row:
        g=''.join([str(j) for j in l])
        h.write(g+'\n')
    h.close 
     
        
        
        
    
    
