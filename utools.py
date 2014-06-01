# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 19:51:31 2013

@author: xianl_000
"""

#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset, MFDataset, num2date, date2num, date2index
import nctools


class ugrid:
    """
    A class for dealing with unstructured grid model output and converting to GNOME format
    Although I use variable names consistent with FVCOM, by passing in a var_map dict the 
    variable names can be customized for SELFE or ADCIRC
    
    Right now the attribute specifying whether the elements are orderd clockwise or counter
    clockwise needs to be manually added before writing to GNOME format (GNOME requres this, 
    but its not often specified in the model output)
        
    """
    
    def __init__(self,FileName=None):
            
        if FileName is not None:
            self.FileName = FileName
            if isinstance(FileName,list):
                self.Dataset = MFDataset(FileName)
            else:
                self.Dataset = Dataset(FileName)
            self.data = dict()
            self.atts = dict()
            
    def get_dimensions(self,var_map):
        
        lat = self.Dataset.variables[var_map['latitude']]
        self.atts['lat'] = dict()
        for an_att in lat.ncattrs():
            self.atts['lat'][an_att] = getattr(lat,an_att)
        self.data['lat'] = lat[:]
        
        lon = self.Dataset.variables[var_map['longitude']]
        self.atts['lon'] = dict()
        for an_att in lon.ncattrs():
            self.atts['lon'][an_att] = getattr(lon,an_att)
        lon = lon[:]
        self.data['lon'] = (lon > 180).choose(lon,lon-360)
        
        self.time = self.Dataset.variables[var_map['time']]
        self.atts['time'] = dict()
        for an_att in self.time.ncattrs():
            self.atts['time'][an_att] = getattr(self.time,an_att)
        self.data['time'] = self.time[:]
    
    
    def get_grid_topo(self,var_map):
        
        nv = self.Dataset.variables[var_map['nodes_surrounding_ele']]
        self.atts['nv'] = dict()
        for an_att in nv.ncattrs():
            self.atts['nv'][an_att] = getattr(nv,an_att)
        self.data['nv'] = nv[:]
        
        nbe = self.Dataset.variables[var_map['eles_surrounding_ele']]
        self.atts['nbe'] = dict()
        for an_att in nbe.ncattrs():
            self.atts['nbe'][an_att] = getattr(nbe,an_att)
        self.data['nbe'] = nbe[:]
        
    def get_data(self,var_map,tindex=None,nindex=None,zindex=0):
    
        ''' 
        var_map is a dict mapping model variable names to common names
        tindex can be used to subset in time --> tindex = [start,stop,step]
        nindex is for subsetting grid -- list of nodes/elements to get u/v on
        zindex is for z layer (FVCOM has surface at zindex = 0, SELFE zindex = -1)
        '''
    
        if tindex:
            self.data['time_ss'] = self.data['time'][tindex[0]:tindex[1]:tindex[2]]
        else:
            tindex = [0,len(self.data['time']),1]
                        
        u = self.Dataset.variables[var_map['u_velocity']]
        self.atts['u'] = dict()
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att)
        v = self.Dataset.variables[var_map['v_velocity']]
        self.atts['v'] = dict()
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att)
        
        if nindex is None:
            if len(u.shape)==3:
                self.data['u'] = u[tindex[0]:tindex[1]:tindex[2],zindex,:]
                self.data['v'] = v[tindex[0]:tindex[1]:tindex[2],zindex,:]    
            elif len(u.shape)==2:
                self.data['u'] = u[tindex[0]:tindex[1]:tindex[2],:]
                self.data['v'] = v[tindex[0]:tindex[1]:tindex[2],:]
            else:
                print "Error:velocity is not 2 or 3 dimensional"
                raise
        else: #Spatial subset -- under development but *mostly* working
            if u.shape[-1] == max(self.data['nbe'].shape):
                # velocities on elements
                id = np.where(np.diff(self.eles_in_ss) > 1)[0]
            else:
                # velocities on nodes
                id = np.where(np.diff(self.nodes_in_ss) > 1)[0]
            id2 = [-1]; id2.extend(id)
            firsttime = 1
            print 'Number of contiguous segments: ', len(id2)+1 
            firsttime = 1
            for ii in range(len(id2)):
                if not np.mod(ii,20):
                    print ii,'of ', len(id2)
                sid = id2[ii]+1
                try:
                    fid = id2[ii+1]
                except IndexError:
                    fid = -1
                if  len(u.shape) == 3:
                    this_u = u[tindex[0]:tindex[1]:tindex[2],0,self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                    this_v = v[tindex[0]:tindex[1]:tindex[2],0,self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                elif len(u.shape)==2:
                    this_u = u[tindex[0]:tindex[1]:tindex[2],self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                    this_v = v[tindex[0]:tindex[1]:tindex[2],self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                else:
                    print "Error:velocity is not 2 or 3 dimensional"
                    raise
                    
                if firsttime:
                    self.data['u'] = this_u.copy()
                    self.data['v'] = this_v.copy()
                    ax = len(this_u)
                    firsttime = 0
                else:  
                    self.data['u'] = np.concatenate((self.data['u'],this_u),axis=ax)
                    self.data['v'] = np.concatenate((self.data['v'],this_v),axis=ax)
        
        # sometimes these are numpy masked arrays -- GNOME can't deal
        if type(self.data['u']) is np.ma.core.MaskedArray:
            self.data['u'] = self.data['u'].data
        if type(self.data['v']) is np.ma.core.MaskedArray:
            self.data['v'] = self.data['v'].data
        
        #self.atts['v']['fill_value'] = self.atts['v']['missing_value']         
        #self.atts['u']['fill_value'] = self.atts['u']['missing_value']         


    def read_bndry_file(self,bnd_file): 
 
        bnd = []
        f = open(bnd_file,'r')
        for line in f:
            vals = [int(val) for val in line.split()]
            bnd.append(vals)
        
        self.data['bnd'] = np.array(bnd)
        self.atts['bnd'] = {'long_name':'Boundary segment information required for GNOME model'} 
    
    def write_bndry_file(self,grid,bndry_file):
          
        '''
        Determine boundary nodes, generate ordered list of segments
        from nv (3 nodes of each element) and nbe (3 elements surrounding
        each element)
        
        grid: 'GOM2','MassB', or, 'NGOFS'
        Output ordered segment list to text file
        
        !!!HARD-CODED to tell which nodes are open-water!! If not know -- all 
        boundaries are specified as open water
        
        '''
        
        #open water nodes -- determined by plotting grid
        if grid.lower() == 'gom2':
            ow1 = 1; ow2 = 60;
        elif grid.lower() == 'massb':
            ow1 = 1; ow2 = 124;
        elif grid.lower() == 'ngofs':
            ow1 = 1; ow2 = 180;
        elif grid.lower() == 'creofs':
            ow = [68408,68409,68410,68411,68412,68414,68604,68605,68606,68607,68608,68791,68792,68793,68962,68963,68964,68965,69130,69131,69132,69133,69303,69304,69305,69479,69481,69669,69670,69671,69672,69674,69675,69866,69867,69868,69869,69870,70062,70063,70064,70065,70271,70272,70489,70490,70704,70705,70927,70928,71144,71346,71520,71683,71844,72001,72154,72281,72377,72462,72532,72583,72631,72676,72720,72765,72810,72851,72897,72939,72981,73023,73061,73099,73138,73178,73215,73251,73283,73313,73346,73381,73417,73453,73454,73481,73502,73523]
        else:
            if grid.lower() != 'subset':
                print 'No grid match -- setting all to open water'
                ow1 = 1; ow2 = 1e6 #set all to open water
        
        if grid.lower() == 'subset':
            lon = self.data['lon_ss']
            lat = self.data['lat_ss']
            nv = self.data['nv_ss']
            nbe = self.data['nbe_ss']
        else:
            lon = self.data['lon']
            lat = self.data['lat']
            nv = self.data['nv']
            nbe = self.data['nbe']
            
        [three,ngl]=nv.shape
        mgl=nv.max()
        
        print 'Finding boundary pts'
        #find all boundary pts (unordered)
        isonb = np.zeros([mgl+1],int)
        for i in range(0,ngl):
            if(np.min(nbe[:,i])==0):
                if(nbe[0,i]==0):
                    isonb[nv[1,i]]=1;isonb[nv[2,i]]=1
                if(nbe[1,i]==0):
                    isonb[nv[0,i]]=1;isonb[nv[2,i]]=1
                if(nbe[2,i]==0):
                    isonb[nv[0,i]]=1;isonb[nv[0,i]]=1
        bnd = list(np.nonzero(isonb)[0])
    
        print 'Ordering boundary pts'
        #now find boundary pts in order starting from a pt in outer bndry which for this
        #case happens to be the first pt in the unorderd bndry list
        pnum=1
        polydict = {}
        obnd = [bnd.pop(1),]
        while len(bnd)>0:
            '''ALGORITHM: from starting node (point) IN OUTER BNDRY, first determine which elements
            contain this point. Make list of all nodes defining these elements. The next bndry node
            will be the one that only appears once in this list.
            '''
            elems = np.nonzero(nv == obnd[-1])[1]
            points = nv[:,elems].flatten() 
            upts = np.unique(points)
            candidates = []
            points = list(points)
            nextpt = []
            for p in upts:
                if points.count(p) == 1:
                    candidates.append(p)
            if len(candidates)>1:
                if len(obnd) == 1:
                    nextpt = candidates[-1]
                else:
                    for c in candidates:
                        if c not in obnd:
                            nextpt = c
                            break
            if nextpt == []:
                polydict[pnum] = obnd
                pnum = pnum+1
                obnd = [bnd.pop(1),]
            else:
                obnd.append(nextpt)
                bnd.remove(nextpt)


        #iterate thru boundary polygons, reverse order if necessary, write to file      
        print 'Checking topology and writing to file' 
        f = open(bndry_file,'w') 
        for key,val in polydict.iteritems():
            area = 0.0
            plon = lon[np.array(val)-1]
            plat = lat[np.array(val)-1]
            j = 0
            for i in range(0,len(plon)-1):
                j += 1
                if j == len(plon)-1: j = 0
                area += (plon[i] + plon[j]) * (plat[i] - plat[j])
            if area > 0:
                polydict[key].reverse
            for i in range(0,len(val)):
                p1 = val[i]
                try:
                    p2 = val[i+1]
                except IndexError:
                    p2 = val[0]
                      
                if grid != 'subset':
                    try:
                        if p1 >= ow1 and p1<=ow2 and p2>=ow1 and p2<=ow2:
                            lw = 1
                        else:
                            lw = 0
                    except NameError: #list of nodes
                        if ow.count(p1) + ow.count(p2) == 2:
                            lw = 1
                        else:
                            lw = 0
                else:
                    lw = 1
                    for sid,seg in enumerate(self.ss_land_bry_segs):
                        if seg.count(p1) + seg.count(p2) == 2: #matching seg
                            lw = 0
                            break
                          
                line = ' '.join(map(str,[p1,p2,key-1,lw]))
                f.write(line + '\n')      
        f.close()    
    
    def write_unstruc_grid(self,ofn):
        
        """
        
        Write GNOME compatible netCDF file (netCDF3) from unstructured (triangular) grid data
        
        """  
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','Triangular')
        
        # test u/v dimensions
        if self.data['u'].shape != self.data['v'].shape:
            print 'u/v dimensions differ'
            raise
        
        # determine if its a subset in time
        t_key = 'time'
        try:
            if self.data['u'].shape[0] == len(self.data['time_ss']):
                t_key = 'time_ss'
        except KeyError:
            if self.data['u'].shape[0] != len(self.data['time']):
                print 'Dimensions of u/v do not match time variable'
                raise
                
        lon_key = 'lon'; lat_key = 'lat'
        nv_key = 'nv'; nbe_key = 'nbe'
        # determine if its a subset of the grid
        try:
            if self.data['u'].shape[-1] == len(self.data['lon_ss']) or \
                self.data['u'].shape[-1] == self.data['nbe_ss'].shape[-1]:
                lon_key = 'lon_ss'; lat_key = 'lat_ss'
                nv_key = 'nv_ss'; nbe_key = 'nbe_ss'
        except KeyError:
            if self.data['u'].shape[-1] != len(self.data['lon']) and \
                self.data['u'].shape[-1] != self.data['nbe'].shape[-1]:
                print 'Dimensions of u/v do not match grid variables'
                raise 
                
        # add Dimensions
        nc.createDimension('time',None)
        nc.createDimension('node',len(self.data[lon_key]))
        nc.createDimension('nele',np.shape(self.data[nbe_key])[1])
        nc.createDimension('nbnd',len(self.data['bnd']))
        nc.createDimension('nbi',4)
        nc.createDimension('three',3)
        #nc.createDimension('sigma',1) #coming soon?
        
        # create variables
        nc_time = nc.createVariable('time','f4',('time',))
        nc_lon = nc.createVariable('lon','f4',('node'))
        nc_lat = nc.createVariable('lat','f4',('node'))
        nc_nbe = nc.createVariable('nbe','int32',('three','nele'))
        nc_nv = nc.createVariable('nv','int32',('three','nele'))
        nc_bnd = nc.createVariable('bnd','int32',('nbnd','nbi'))
        
        if self.data['u'].shape[-1] == len(self.data[lon_key]): #velocities on nodes
            nc_u = nc.createVariable('u','f4',('time','node'))
            nc_v = nc.createVariable('v','f4',('time','node'))
        else: #velocities on elements
            nc_u = nc.createVariable('u','f4',('time','nele'))
            nc_v = nc.createVariable('v','f4',('time','nele'))
        
        #adjust time if necessary
        ref_time = self.atts['time']['units'].split(' since ')[1]
        ref_year = int(ref_time[0:4])
        if ref_year < 1970:
            print 'Adjusting reference time'
            self.data[t_key],self.atts['time']['units'] = \
                nctools.adjust_time(self.data[t_key],self.atts['time']['units'])
        
        #add data to netcdf file
        nc_time[:] = self.data[t_key]
        nc_lon[:] = self.data[lon_key]
        nc_lat[:] = self.data[lat_key]
        nc_u[:] = self.data['u']
        nc_v[:] = self.data['v']
        nc_bnd[:] = self.data['bnd']
        nc_nbe[:] = self.data[nbe_key]
        nc_nv[:] = self.data[nv_key]
        
        #add variable attributes to netcdf file
        for an_att in self.atts['time'].iteritems():
           setattr(nc_time,an_att[0],an_att[1])
        
        for an_att in self.atts['lon'].iteritems():
            setattr(nc_lon,an_att[0],an_att[1])
        
        for an_att in self.atts['lat'].iteritems():
            setattr(nc_lat,an_att[0],an_att[1])
        
        for an_att in self.atts['bnd'].iteritems():
            setattr(nc_bnd,an_att[0],an_att[1])


        for an_att in self.atts['nbe'].iteritems():
            setattr(nc_nbe,an_att[0],an_att[1])


        for an_att in self.atts['nv'].iteritems():
            setattr(nc_nv,an_att[0],an_att[1])
        
        for an_att in self.atts['u'].iteritems():
           setattr(nc_u,an_att[0],an_att[1])
        
        for an_att in self.atts['v'].iteritems():
           setattr(nc_v,an_att[0],an_att[1])
        
        nc.close()
    
    def find_nodes_eles_in_ss(self,nl,sl,wl,el):
       
        print 'Total number of eles: ', self.data['nbe'].shape[1]
        print 'Total number of nodes: ', self.data['nv'].shape[1]
        
        #returns lists of eles and nodes, plus truncated and edited topology arrays (nbe, nv)
        subset_lat = np.nonzero(np.logical_and(self.data['lat']>=sl,self.data['lat']<=nl))[0]
        subset_lon = np.nonzero(np.logical_and(self.data['lon']>=wl,self.data['lon']<=el))[0]
        self.nodes_in_ss = np.intersect1d(subset_lat,subset_lon) + 1 #node numbering starts at 1 (subtract one for indexing lon/lat)
        
        self.eles_in_ss = []
        nv_ss = []; nbe_ss = []; 
        
        #determine which nodes are in subset boundary and elements with all nodes in ss
        print 'Finding nodes and entire elements in ss'
        for ii, ele in enumerate(self.data['nv'].transpose()):
            #if all of the nodes are in subset domain keep -- otherwise get rid of it
            if (self.nodes_in_ss == ele[0]).any() and (self.nodes_in_ss == ele[1]).any() \
             and (self.nodes_in_ss == ele[2]).any():
                nv_ss.append(ele)
                nbe_ss.append(self.data['nbe'][:,ii])
                self.eles_in_ss.append(ii+1) #ele numbering starts at 1
            else:
                pass
        print 'Number of eles in ss: ', len(self.eles_in_ss)
        print 'Number of nodes in ss: ', len(self.nodes_in_ss)
              
        nbe_ss = np.array(nbe_ss).transpose()
        nv_ss = np.array(nv_ss).transpose()
        self.eles_in_ss = np.array(self.eles_in_ss)


        print 'Remapping nodes and elements'
        #now remap nbe2, nv2 to number of remaining nodes, and elements
        nv_ssr = nv_ss.copy()
        nbe_ssr = nbe_ss.copy()
        for ii in range(len(self.eles_in_ss)):
            for jj in range(3):
                nid = np.searchsorted(self.nodes_in_ss, nv_ss[jj,ii], side='left')
                nv_ssr[jj,ii] = nid+1
                if nbe_ss[jj,ii] != 0:
                    eid = np.searchsorted(self.eles_in_ss, nbe_ss[jj,ii], side='left')
                    if self.eles_in_ss[eid] != nbe_ss[jj,ii]:
                        nbe_ssr[jj,ii] = 0
                    else:
                        nbe_ssr[jj,ii] = eid+1
                                   
        self.data['nbe_ss'] = nbe_ssr
        self.data['nv_ss'] = nv_ssr
        self.data['lon_ss'] = self.data['lon'][self.nodes_in_ss-1]
        self.data['lat_ss'] = self.data['lat'][self.nodes_in_ss-1]
        
    def remap_bry_nodes(self,bndry_file):
      
        print 'Remapping boundary segs to new subset node numbers'
  
        f = open(bndry_file)
        ss_bry_land_segs = []
        for line in f:
            node1,node2,bnumber,flag = map(int,line.split())
            node1_id = np.where(self.nodes_in_ss == node1)
            node2_id = np.where(self.nodes_in_ss == node2)
            if len(node1_id) > 0 and len(node2_id) > 0 and flag == 0:
                ss_bry_land_segs.append([node1_id[0]+1,node2_id[0]+1,flag])
    
        f.close()
    
        return ss_bry_land_segs
    
    def build_face_face_connectivity(self):
        """
        builds the triangular connectivity array (nbe)
        essentially giving the neighbors of each triangle (face)
        """        
        num_vertices = 3
        num_faces = len(self.data['nv'].transpose())
        face_face = np.zeros( (num_faces, num_vertices), dtype=np.int32  )
        face_face += -1 # fill with -1


        # loop through all the triangles to find the matching edges:
        edges = {} # dict to store the edges in 
        for i, face in enumerate(self.data['nv'].transpose()):
            # loop through edges:
            for j in range(num_vertices):
                edge = (face[j-1], face[j])
                if edge[0] > edge[1]: # sort the node numbers
                    edge = (edge[1], edge[0]) 
                # see if it is already in there
                prev_edge = edges.pop(edge, None)
                if prev_edge is not None:
                    face_num, edge_num = prev_edge
                    face_face[i,j] = face_num
                    face_face[face_num, edge_num] = i
                else:
                    edges[edge] = (i, j)
        nbe = (face_face + 1)
        nbe = nbe[:,[2,0,1]].transpose()
        self.data['nbe'] = nbe
        self.atts['nbe'] = {'long_name': 'elements surrounding element'}

