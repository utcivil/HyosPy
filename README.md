HyosPy
======

HyosPy (Hydrodynamic and oil spill Python) is a model coupling system developed by Dr. Ben R. Hodges' research team in CRWR at UT-Austin. The HyosPy is able to automate multiple simulations of a hydrodynamic model (SELFE) using real-time downloaded hindcast/foreast data, then couples the results to an oil spill transport model (PyGNOME or RK4) and provides visualization/animation in the 2D Google Map and 3D Google Earth GIS.

The current version of HyosPy has the following new features:

1. The structure of HyosPy has updated by adding a upper wrapper that controls all the necessary applications, such as changing parameters, executing models, and etc., so that people don't need to understand the logic and dig into the code but just use the upper wrapper to do everything they need.

2. The original RK4 has been replaced by the PyGNOME (the Linux version of GNOME) and it could run much more faster (more than 100X times faster) than the RK4.

3. The PyGNOME was initiated based on the forecast current fields by SELFE and also the wind force.

4. The original 2D Google Map visualization tool has been updated to the 3D Google Earth interface, with which multiple oil spill forecast trajectories could be animated simultaneously. You may want to check it out the sample animation we have done at: https://www.youtube.com/watch?v=-5yNl0K6wJ4

Requirements
-----
1. Parallel SELFE with MPI protocol (v3.1dc) : http://www.stccmop.org/knowledge_transfer/software/selfe/sourcecode
2. PyGNOME (Linux version of GNOME) : https://github.com/NOAA-ORR-ERD/GNOME2
3. NetCDF library (4.3.0 or higher) : http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
4. Python (2.7.2 or higher) : https://www.python.org/download/

Tutorial
--------
Below is an example of HyosPy structure with 12 SELFE embedded:


HyosPy/

      /1/GNOME/
              /current.txt
              /coast.bna
        /hotstart.in
        /run.sh
        /autocombine_MPI_elfe.pl
        /machines
        /all the other mandatory SELFE inputs except 'param.in', 'wind.th', 'elev.th', and 'flux.th'
        
      /2/GNOME/
              /current.txt
              /coast.bna
        /hotstart.in
        /run.sh
        /autocombine_MPI_elfe.pl
        /machines
        /all the other mandatory SELFE inputs except 'param.in', 'wind.th', 'elev.th', and 'flux.th'
        
      ...
      
      /12/GNOME/
               /current.txt
               /coast.bna
         /hotstart.in
         /run.sh
         /autocombine_MPI_elfe.pl
         /machines
         /all the other mandatory SELFE inputs except 'param.in', 'wind.th', 'elev.th', and 'flux.th'
      
      /hydro_wrapper.py
      /data.py
      /GNOME_run.py
      /oilspill_wrapper.py
      /multiple_models_with_multiple_points_v1.py
      /upper_hydro_wrapper.py
      /upper_oilspill_wrapper.py
      /GE_animation.txt
      /utools.py
      /nctools.py
      /169_hvel.nc
      /param.in
      /pyselfe_v1.py
      /javascript.txt
      /numpyIO.py
      /drag.gr3
      
        
       
        
     
