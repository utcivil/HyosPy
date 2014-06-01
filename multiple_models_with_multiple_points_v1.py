# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 14:26:23 2013
"""
import os
import string
import heapq
import pyselfe_v1
import math
base_dir = os.path.dirname(__file__)

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


# Get all of the coordinate info and save them in the list of location

file = open(base_dir+'drag.gr3','r')
location=[]
for s in file.readlines():    # creat the 2D list to store the data
    column=[]
    line=s.split()
    for field in line: column.append(field)
    location.append(column)    
file.close

del location[0:2];del location[23286:]

for i in range(len(location)):
    del location[i][0]
for j in range(len(location)):
    del location[j][2]


# read velocity at point(x,y) at time t
def readVel(t,x,y):
    
    # find the 3 nearest nodes to (x,y)

    distance=[]
    for m in range(len(location)):
        aa=string.atof(location[m][0])-x ; bb=string.atof(location[m][1])-y
        distance.append(aa*aa+bb*bb)
    
    nearest=heapq.nsmallest (3, distance)   #find 3 smalleset number and save them at list 'nearest'

    n1=distance.index(nearest[0]); n2=distance.index(nearest[1]); n3=distance.index(nearest[2])

    x01=string.atof(location[n1][0]);y01=string.atof(location[n1][1])
    x02=string.atof(location[n2][0]);y02=string.atof(location[n2][1])
    x03=string.atof(location[n3][0]);y03=string.atof(location[n3][1])
    
    # save the velocity in x direction at time step n in ui, y direction at time step n in wi, i=1,2,3
    if isinstance(t,int):
        
        u01=mdata[t,n1,5,0];w01=mdata[t,n1,5,1]
        u02=mdata[t,n2,5,0];w02=mdata[t,n2,5,1]
        u03=mdata[t,n3,5,0];w03=mdata[t,n3,5,1]
        
    else :
        u01=(mdata[int(t),n1,5,0]+mdata[int(t)+1,n1,5,0])/2
        w01=(mdata[int(t),n1,5,1]+mdata[int(t)+1,n1,5,1])/2
        u02=(mdata[int(t),n2,5,0]+mdata[int(t)+1,n2,5,0])/2
        w02=(mdata[int(t),n2,5,1]+mdata[int(t)+1,n2,5,1])/2
        u03=(mdata[int(t),n3,5,0]+mdata[int(t)+1,n3,5,0])/2
        w03=(mdata[int(t),n3,5,1]+mdata[int(t)+1,n3,5,1])/2

    # save the distance between (x,y) and the three nearest points to Li,i=1,2,3
    L1=math.sqrt((x01-x)*(x01-x)+(y01-y)*(y01-y))
    L2=math.sqrt((x02-x)*(x02-x)+(y02-y)*(y02-y))
    L3=math.sqrt((x03-x)*(x03-x)+(y03-y)*(y03-y))
    Lsum=L1+L2+L3

    # interpolate to find velocity component of (x,y) at time step n and save them as u,w
    u=L1/Lsum*u01+L2/Lsum*u02+L3/Lsum*u03 ; w=L1/Lsum*w01+L2/Lsum*w02+L3/Lsum*w03
    return t,u,w

# set how many SELFE folders to read
number=12 
selfedir=[]
for i in range(1,number+1):
    selfedir.append(base_dir+str(i)+'/outputs/')

dt=900 ; h=900 
traceroutput=[] ; gradientoutput=[] 

# main loop(loop with different models)
for nmodels in range(number):
    
    selfe = pyselfe_v1.Dataset(selfedir[nmodels]+'169_hvel.64')
    [t, t_iter, eta, dp, mdata] = selfe.read_time_series('hvel.64', nfiles=48, sfile=169, datadir=selfedir[nmodels])
    
    length=len(mdata[:,7,5,0])
    utm_x=666393 ; utm_y=3076780
    for a in range(13):
    
        x=[None]*(length+1) ; y=[None]*(length+1) ; u=[None]*length ; w=[None]*length 
        tracer=[] ; element=[] ; velocity=[] ; uv=[] ; gradient=[] ; xy=[] ; vel0=[] ; vel1=[] ; vel2=[] ; vel3=[]
    
    # set N particles at different starting places
        if a ==0:
            x[0]=utm_x;y[0]=utm_y
        elif a ==1:
            x[0]=utm_x-10;y[0]=utm_y
        elif a ==2:
            x[0]=utm_x-20;y[0]=utm_y
        elif a ==3:
            x[0]=utm_x+10;y[0]=utm_y
        elif a ==4:
            x[0]=utm_x+20;y[0]=utm_y
        elif a ==5:
            x[0]=utm_x;y[0]=utm_y+10
        elif a ==6:
            x[0]=utm_x;y[0]=utm_y+20
        elif a ==7:
            x[0]=utm_x;y[0]=utm_y-10
        elif a ==8:
            x[0]=utm_x;y[0]=utm_y-20
        elif a ==9:
            x[0]=utm_x-10;y[0]=utm_y+10
        elif a ==10:
            x[0]=utm_x-10;y[0]=utm_y-10
        elif a ==11:
            x[0]=utm_x+10;y[0]=utm_y+10   
        else :
            x[0]=utm_x+10;y[0]=utm_y-10
    
        for n in range(length-1):    
          
            vel0=list(readVel(n,x[n],y[n])) 
            u[n]=vel0[1];w[n]=vel0[2]

            dx=dt*u[n] ; dy=dt*w[n]
            uv=[u[n],w[n],dx,dy]
            velocity.append(uv) 

    # transport the oil spill point by using forward Euler algorithm
    
        #x[n+1]=x[n]+dx ; y[n+1]=y[n]+dy

    # transport the oil spill point by using RK4 algorithm
            k1x=h*u[n] ; k1y=h*w[n]
            vel1=list(readVel(n+0.5,x[n]+k1x/2,y[n]+k1y/2))     #read velocity f(t+h/2,x+k1x/2),f(t+h/2,y+k1y/2) at time t+h/2 at point(x+k1x/2,y+k1y/2)
            k2x=h*vel1[1] ; k2y=h*vel1[2]
            vel2=list(readVel(n+0.5,x[n]+k2x/2,y[n]+k2y/2))
            k3x=h*vel2[1] ; k3y=h*vel2[2]
            vel3=list(readVel(n+1,x[n]+k3x,y[n]+k3y))
            k4x=h*vel3[1] ; k4y=h*vel3[2]
        
            x[n+1]=x[n]+(k1x+2*k2x+2*k3x+k4x)/6 ; y[n+1]=y[n]+(k1y+2*k2y+2*k3y+k4y)/6
    
            element=[(n+1)*900,x[n],y[n]] 
            tracer.append(element)
           
    # velocity gradient

        for t in range(len(velocity)-1):
            du=velocity[t+1][0]-velocity[t][0] ; dw=velocity[t+1][1]-velocity[t][1]
            xy=[du/velocity[t][2],du/velocity[t][3],dw/velocity[t][2],dw/velocity[t][3]]
            gradient.append(xy)

        traceroutput.append(tracer)
        gradientoutput.append(gradient)


lonlat=[]      
for nmodelpoint in range(number*13):  # nmodelpoint is the number of "# of model * # of point" you determine to set up 11
    for nposition in range(len(traceroutput[nmodelpoint])):
        aaa=utmToLatLng(14,traceroutput[nmodelpoint][nposition][1],traceroutput[nmodelpoint][nposition][2],northernHemisphere=True)
        
        lonlat.append([aaa[0],aaa[1]])  

coordinates=str(lonlat)

file1 = open('./javascript.txt','r')     # open the javascript
row=[]
for s in file1.readlines():        # read the javascript in lines and save them in row[]
    row.append(s)

row[11]='    '+'var'+' '+'locations'+'='+coordinates     #change the javascript code with the coordinates you want to mark. in this case, it is the tracers' coordinates.
   
h=open('./RK4_GM.html','w')
for l in row:
    g=''.join([str(j) for j in l])
    h.write(g+'\n')
h.close
