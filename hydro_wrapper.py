# -*- coding: utf-8 -*-
"""
Hydrodynamic wrapper tools

"""

import os
import subprocess
import time
import shutil

base_dir = os.path.dirname(__file__)

def timer(seconds):
    
    time.sleep(seconds)
    
def runSELFE(number,interval):
    """
    run "number" of SELFE with "interval" time interval;
    clean up unnecessary binary outputs;
    copy useful data
    """
    
    for i in range(number):
        path=base_dir+'/'+str(i+1)
        subprocess.Popen(['python','data.py'],cwd=path)
        time.sleep(interval)
        
        alldata = os.listdir(base_dir+'/'+str(i+1)+'/outputs/')
        
        for data in alldata:
            if data.endswith('_hotstart'):
                os.remove(base_dir+'/'+str(i+1)+'/outputs/'+data)
            if data.endswith('_0000_elev.61'):
                os.remove(base_dir+'/'+str(i+1)+'/outputs/'+data)
            if data.endswith('_0001_elev.61'):
                os.remove(base_dir+'/'+str(i+1)+'/outputs/'+data)
            if data.endswith('_0000_wind.62'):
                os.remove(base_dir+'/'+str(i+1)+'/outputs/'+data)
            if data.endswith('_0001_wind.62'):
                os.remove(base_dir+'/'+str(i+1)+'/outputs/'+data)
            if data.endswith('_0000_hvel.64'):
                os.remove(base_dir+'/'+str(i+1)+'/outputs/'+data)
            if data.endswith('_0001_hvel.64'):
                os.remove(base_dir+'/'+str(i+1)+'/outputs/'+data)

        alldata1=os.listdir(base_dir+'/'+str(i+1)+'/outputs/')

        for data in alldata1:         
            if data.endswith('.64'):
                for j in range(i+1,number):
                    shutil.copy(base_dir+'/'+str(i+1)+'/outputs/'+data,base_dir+'/'+str(j+1)+'/outputs')
        
        
def change_param(number,run_time, time_step, elev, wind, hvel):
    
    '''
    number is SELFE #;
    run_time in days; time_step in seconds;
    elev,wind,hvel has value 1(switch on) or 0(switch off)
    '''    
    
    param=open(base_dir+'/param.in','r')
    row=[]
    for i in param.readlines():
        row.append(i)
        
    row[164]='  '+'rnday = '+str(run_time)+' '+'!total run time in days'+'\n'
    row[167]='  '+'dt = '+str(time_step)+'.'+' '+'!Time step in sec'+'\n'
    row[277]='  '+'elev.61 = '+str(elev)+' '+'!0: off; 1: on'+'\n'
    row[289]='  '+'wind.62 = '+str(wind)+'\n'
    row[302]='  '+'hvel.64 = '+str(hvel)+'\n'
    
    for j in range(number):
        h=open(base_dir+'/'+str(j+1)+'/param.in','w')
        for x in row:
            g=''.join([str(n) for n in x])
            h.write(g)
        h.close 
        
        
def change_data_time(number,start_mon,start_day,start_yr,end_mon,end_day,end_yr):
    
    '''
    number is SELFE #;
    customize simulation time period
    '''
    
    autodata=open(base_dir+'/data.py','r')
    row=[]
    for i in autodata.readlines():
        row.append(i)
    
    # note: the start time = simulation start time-7)Day due to a 7 days model spin up
    row[15]='timeperiod=['+"'start:','"+str('%02d' % start_mon)+"','"+str('%02d' % start_day)+"','"+str(start_yr)+"','end:','"+str('%02d' % end_mon)+"','"+str('%02d' % end_day)+"','"+str(end_yr)+"'"+']'
    
    for j in range(number):
        
        row[140]='f=open(dire+'+"'wind.th'"+','+"'w')"+'\n'
        row[264]='h=open(dire+'+"'elev.th'"+','+"'w')"+'\n'
        row[369]='flux_th=open(dire+'+"'flux.th'"+','+"'w')"+'\n'
        
        h=open(base_dir+'/'+str(j+1)+'/data.py','w')
        for x in row:
            g=''.join([str(n) for n in x])
            h.write(g)
        h.close 
    
    
