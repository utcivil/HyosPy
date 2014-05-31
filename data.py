# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 20:24:45 2013

@author: xianl_000
"""

######## This is a combined module of wind, tide and river flow hindcast/forecast automation
######## Separate modules for wind, tide and inflow hindcast/forecast are also available at (our github or our research website)
############################################################################################



## simulation time period: 7/21 0:00-7/23 0:00
## because of model a 7 days model spin-up (7 more days hindcasts), the data start time should be 7 days less than the model simulation start time
timeperiod=['start:','07','14','2013','end:','07','23','2013']

####### Wind ##########
# Download wind forecast data from "http://seawater.tamu.edu/tglopu/twdb_lc.tar" and save it at "D:/twdb_1c.tar".
import os
dire=os.path.dirname(__file__)
import urllib
os.mkdir(dire+'wind')
url = 'http://seawater.tamu.edu/tglopu/twdb_lc.tar'
path = dire+'wind/twdb_1c.tar'
data = urllib.urlopen(url).read()
f = file(path,'wb')
f.write(data)
f.close()

# decide whether the download process is over and unzip the tar file to 'D:/wind' when the download is over.

import tarfile
from contextlib import closing
path = dire+'wind/twdb_1c.tar'
if os.path.isfile(dire+'wind/twdb_1c.tar'):
   os.mkdir(dire+'wind/wind_data')
   with closing(tarfile.open(dire+'wind/twdb_1c.tar','r')) as t:
        t.extractall(dire+'wind/wind_data')

# Step one to create wind.th--delet the first few lines of twdb051.wndq file and save it as wind1.th

string1='*'
a=open(dire+'wind/wind_data/twdb051.wndq', 'r').readlines()   # read file in lines
open(dire+'wind/wind.txt','w').write('')   # empty wind.txt file
for x in a:
    if string1 in x:     # if c in x line, skip
         continue
    open(dire+'wind/wind.txt','a').write(x)

string2='days'
d=open(dire+'wind/wind.txt','r').readlines()
open(dire+'wind/wind1.th','w').write('')
for y in d:
    if string2 in y:
        continue
    open(dire+'wind/wind1.th','a').write(y)

### UPDATE_5_29_2014 ###
## create wind.WND for PyGNOME ##

file000 = open(dire+'wind/wind1.th','r')
row000=[]
for s in file000.readlines():    
    column=[]
    line=s.split()
    for field in line: column.append(field)
    row000.append(column)
file000.close

for i in range(len(row000)):
    row000[i].append(i)

for i in range(len(row000)):
    ## wind direction
    row000[i][6]=row000[i][5]
    
    ## wind magnitude
    row000[i][5]=str(float(row000[i][4])*0.447)+','
    
    ## time
    row000[i][4]='00,'
    row000[i][3]=row000[i][3]+','
    row000[i][1]=str('%01d' % int(row000[i][1]))+','
    stock=row000[i][0]
    row000[i][0]=str('%01d' % int(row000[i][2]))+','
    row000[i][2]=stock[-2:]+','

ff=open(dire+'GNOME/wind.WND','w')

ff.write('twdb051'+'\n')
ff.write('27.798 -97.311'+'\n')
ff.write('mps'+'\n')
ff.write('LTime'+'\n')
ff.write('0,0,0,0,0,0,0,0'+'\n')
for i in (row000):
    k='   '.join([str(j) for j in i])
    ff.write(k+'\n')
ff.close
ff.close    

# Step two of create wind.th---delet the first few useless columns of wind1.th and save it in the list 'row'

file0 = open(dire+'wind/wind1.th','r')
row=[]
for s in file0.readlines():    # creat the 2D list to store the data
    column=[]
    line=s.split()
    for field in line: column.append(field)
    row.append(column)    
file0.close
  
# make record of what time period of wind data are required
# note that there is a 7 day's spin-up, so the the required time period should add a 7 days spin-up to form a hindcast

days=[]
for day in range(len(row)):
    days.append(row[day][1]+row[day][2])

start=days.index(timeperiod[1]+timeperiod[2])     # input start time in terms of "Mon+Day"
end=days.index(timeperiod[5]+timeperiod[6])     # input end time in unit "Mon+Day"  (note: the start time = simulation start time-7)Day 

for m in range(len(row)):
     row[m][0:4]=[]


# Step three of create wind.th---transfer data format and save it as 'D:/wind.th'

import string       
import math
    
for n in range(len(row)):
    aa=string.atof(row[n][0]);bb=string.atof(row[n][1])
    row[n][0]=0.447*aa*math.sin(bb*math.pi/180);row[n][1]=0.447*aa*math.cos(bb*math.pi/180)

winds=[]
for wind in range(start,end+8):
    winds.append(row[wind])
    

f=open(dire+'SELFE/1/wind.th','w')
for i in (winds):
    k='   '.join([str(j) for j in i])
    f.write(k+'\n')
f.close

import shutil
# shutil.rmtree(dire+'wind') # whether or not to remove the raw data

####### Tide #########
# Download pwl and harmwl data from TCOON and save it to the location as path
# change time section of url for different time in order to get different elev.th

os.mkdir(dire+'tide')

url = 'http://lighthouse.tamucc.edu/pd?stnlist=014&serlist=pwl%2Charmwl&when='+timeperiod[1]+'.'+timeperiod[2]+'.'+timeperiod[3]+'-'+timeperiod[5]+'.'+timeperiod[6]+'.'+timeperiod[7]+'&whentz=UTC0&-action=c&unit=metric&elev=msl'
path = dire+'tide/elevation.txt'
data=urllib.urlopen(url).read()
f = file(path,'wb')
f.write(data)
f.close

# Step one to create elev.th---delet the first few lines and NA lines of elevation.txt file and save it as elev1.th

string1='#'
a=open(dire+'tide/elevation.txt', 'r').readlines()   # read file in lines
open(dire+'tide/elev.txt','w').write('')   # empty elevation.txt file
for x in a:
    if string1 in x:     # if c in x line, skip
         continue
    open(dire+'tide/elev.txt','a').write(x)

string2='NA'
d=open(dire+'tide/elev.txt','r').readlines()
open(dire+'tide/elev1.th','w').write('')
for y in d:
    if string2 in y:
        continue
    open(dire+'tide/elev1.th','a').write(y)
    
    
# Step two of create elev.th---delet the first column of elev1.th and save it in the list 'row'

file1 = open(dire+'tide/elev1.th','r')
row=[]
for s in file1.readlines():    # creat the 2D list to store the data
    column=[]
    line=s.split()
    for field in line: column.append(field)
    row.append(column)    
file1.close

row.pop()    # delet the last useless element
for m in range(len(row)):
     del row[m][0] 

# calculate the mean difference between pwl and harmwl with N = 30 (3 hours)
     
sum=0
for n in range(len(row),len(row)-30,-1):
    aa=string.atof(row[n-1][0]);bb=string.atof(row[n-1][1])
    diff=aa-bb
    sum=sum+diff
x=sum/len(row)

# re-open elev.txt to edit it into selfe input

file1 = open(dire+'tide/elev.txt','r')
row1=[]

for q in file1.readlines():    # creat the 2D list to store the data
    column1=[]
    line1=q.split()
    for field1 in line1: column1.append(field1)
    row1.append(column1) 

file1.close
row1.pop() 


# find the index of the first 'NA'

for aa in range(len(row1)):
    del row1[aa][0];del row1[aa][1]  # because the first del delete the first column of the data, which make the data format change into [[],[]], so to delect the 'third' column of the 'raw' data, the second del's index becomes to be 1.
    
y=row1.index(['NA'])

file2 = open(dire+'tide/elev.txt','r')
row2=[]

for q in file2.readlines():    # creat the 2D list to store the data
    column2=[]
    line2=q.split()
    for field2 in line2: column2.append(field2)
    row2.append(column2) 

file2.close
row2.pop() 

for rr in range(len(row2)):  
    if row2[rr][1]== 'RM':                
        row2[rr][1] = row2[rr][2]

# add x with harmwl to creat the forecast part of pwl, set the return interval as 72 hours (72*60/6=720)

for i in range(y,len(row2)):
    row2[i][1]=string.atof(row2[i][2])+x*(1-i/720)

for k in range(y):
    row2[k][1]=string.atof(row2[k][1])

for m in range(len(row2)):     # change time step and delet the third column to suit the selfe input format
    row2[m][0]=360*m
    del row2[m][2]

# interpolate in order to suit the selfe input

for t in range(len(row2)-1):
    row2.insert(2*t+1,[(2*t+1)*180,(row2[2*t][1]+row2[2*t+1][1])/2])

del row2[0]   

# save elev.th at specific location 'F:/elev.th'

h=open(dire+'SELFE/1/elev.th','w')
for l in row2:
    g='   '.join([str(j) for j in l])
    h.write(g+'\n')
h.close 

#shutil.rmtree(dire+'tide') # whether or not to remove the raw data

####### River flow #########
# USGS websites for river flow data 5(available)/11(total)

os.mkdir(dire+'flux')
# time period 
time='http://waterdata.usgs.gov/tx/nwis/uv?cb_00060=on&format=rdb&period=&begin_date='+timeperiod[3]+'-'+timeperiod[1]+'-'+timeperiod[2]+'&end_date='+timeperiod[7]+'-'+timeperiod[5]+'-'+timeperiod[6]+'&site_no='
# site number
copano2='08189200'
mission3='08189500'
aransas4='08189700'
nueces5='08211200'
oso8='08211520'

sitelist=[]
sitelist.append(copano2);sitelist.append(mission3);sitelist.append(aransas4)
sitelist.append(nueces5);sitelist.append(oso8)

weblist=[]
weblist.append(time+copano2);weblist.append(time+mission3);weblist.append(time+aransas4)
weblist.append(time+nueces5);weblist.append(time+oso8)

# read in 10 days flux.th (historical)
j = open(dire+'prototype_flux.th','r')
flux0=[]
for ss in j.readlines():    # creat the 2D list to store the data
        columnn=[]
        lines=ss.split()
        for fields in lines: columnn.append(fields)
        flux0.append(columnn)    
j.close

del flux0[4800:]   # only use 10(7 days hindcast(spin up) and 3 days forecast) days data

string1='#'
fflux=[]

# main loop for creating flux.th
for n in range(len(weblist)):
    
    url = weblist[n]
    path = dire+'flux/'+sitelist[n]+'.txt'
    data=urllib.urlopen(url).read()
    f = file(path,'wb')
    f.write(data)
    f.close
    
    # Delete the first few line with #
    
    a=open(dire+'flux/'+sitelist[n]+'.txt', 'r').readlines()   # read file in lines
    open(dire+'flux/'+sitelist[n]+'a.txt','w').write('')   # empty elevation.txt file
    for x in a:
        if string1 in x:     # if c in x line, skip
             continue
        open(dire+'flux/'+sitelist[n]+'a.txt','a').write(x)

    # put data into list "row"
    b = open(dire+'flux/'+sitelist[n]+'a.txt','r')
    row=[]
    for s in b.readlines():    # creat the 2D list to store the data
        column=[]
        line=s.split()
        for field in line: column.append(field)
        row.append(column)    
    b.close
     
    del row[0]; del row[0]; row.pop() # delete useless stuff
    
    # store river discharge in the list "flux"
    flux=[]
    for i in range(len(row)):
        flux.append(row[i][5])
        flux[i]=string.atof(flux[i])*0.02831685   # convert string to interger and unit of m3/s
        
    # interpolate in order to suit the selfe input
    for m in range(len(flux)-1):
        aa=5*m;bb=5*m+1
        flux[bb:bb]=[(flux[aa]+flux[bb])/2,(flux[aa]+flux[bb])/2,(flux[aa]+flux[bb])/2,(flux[aa]+flux[bb])/2]

    # extend forecast data
    
    number=len(flux0)-len(flux)   # set how many data to extend
    extend=flux[len(flux)-1]  # use the last available data to extend
    for y in range(number):
        flux.append(extend)
    
    fflux.append(flux)

# replace part of the historical data with the available forecast data

for v in range(len(flux0)):
    flux0[v][2]=fflux[0][v]   # copano2
    flux0[v][3]=fflux[1][v]   # mission3
    flux0[v][4]=fflux[2][v]   # aransas4
    flux0[v][5]=fflux[3][v]   # nueces5
    flux0[v][8]=fflux[4][v]   # oso8

# put the final flux.th in the SELFE direction
flux_th=open(dire+'SELFE/1/flux.th','w')
for i in (flux0):
    k='   '.join([str(j) for j in i])
    flux_th.write(k+'\n')
flux_th.close

# shutil.rmtree(dire+'flux') # whether or not to remove the raw data

###### run SELFE  #####

import subprocess
subprocess.Popen('./run.sh',shell=False)



