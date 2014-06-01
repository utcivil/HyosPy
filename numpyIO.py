'''
Author: Benyang Tang, btang(noSpam)@pacific.jpl.nasa.gov
Modified by: Nickolas Fotopoulos, nvf(noSpam)@mit.edu

========
Document:
========

This is an implementation in python of the c version of numpyio by Travis
Oliphant. The advantage of using numpyIO instead of numpyio is the ease of
installation: All you have to do it to put this file into a directory of your
python path.

I have not benchmarked this module. I expect numpyIO to be slower than numpyio,
but not by much.

Only 2 functions were implemented: fread and fwrite.

The interfaces of the 2 functions here are exactly the same as those of numpyio.
If you have codes using numpyio, you don't have to change anything to call
fread and fwrite, except changing numpyio to numpyIO. 

Here is the document of fread and fwrite from Travis Oliphant:

********************************************************************

g = numpyio.fread( fid, Num, read_type { mem_type, byteswap})

     fid =       open file pointer object (i.e. from fid = open("filename") )
     Num =       number of elements to read of type read_type
     read_type = a character in 'cb1silfdFD' (PyArray types)
                 describing how to interpret bytes on disk.
OPTIONAL
     mem_type =  a character (PyArray type) describing what kind of
                 PyArray to return in g.   Default = read_type
     byteswap =  0 for no byteswapping or a 1 to byteswap (to handle
                 different endianness).    Default = 0

************************************************************************

numpyio.fwrite( fid, Num, myarray { write_type, byteswap} )
 
     fid =       open file stream
     Num =       number of elements to write
     myarray =   NumPy array holding the data to write (will be
                 written as if ravel(myarray) was passed)
OPTIONAL
     write_type = character ('cb1silfdFD') describing how to write the 
                  data (what datatype to use)  Default = type of
                  myarray.
     byteswap =   0 or 1 to determine if byteswapping occurs on write.
                  Default = 0.

************************************************************************

One convenience I have coded into fread: If you want to read all bytes from a
file, you can just specify the second argument Num as None.

Similarly for fwrite: If you want write all elements of myArray to a file, you
can just specify the second argument Num as None.
  


====
test:
====

import numpy
import numpyIO
a = numpy.array([1.0, 1.1, 1.2]).astype('d')

fn = '/home/btang/tmp/temp99'
fid = open(fn,'w')
numpyIO.fwrite(fid,None,a,write_type='d', byteswap=1)

fid.close()

fid = open(fn,'r')
a1 = numpyIO.fread(fid,None,'d','f',byteswap=1)
fid.close()

a1

=======
History:
=======
2005-09-02: Added a few new readTypes to support Matlab R14 mat-files.  Not
            comprehensive.
2003-03-31: Coded and tested.
'''

import numpy

#=======================
def fread(fid, num=None, readType='f', mem_type=None, byteswap=0):
#=======================

#=== how many bytes per number
  if readType=='1' or readType=='b' or readType=='c' or readType=='B':
    byteSize = 1
  elif readType=='s' or readType=='w' or readType=='h' or readType=='H':
    byteSize = 2
  elif readType=='f' or readType=='i' or readType=='I'\
    or readType=='l' or readType=='u':
    byteSize = 4
  elif readType=='d' or readType=='F':
    byteSize = 8
  elif readType=='D':
    byteSize = 16
  else:
    print readType
    raise TypeError

#=== figure out num
  if not num:
    import os
    fName = fid.name
    fileSize = os.path.getsize(fName)
    num = fileSize / byteSize

#=== figure out mem_type
  if not mem_type:
    mem_type = readType

#=== read in
  nByte = num * byteSize

  if byteswap==1:
    a1 = numpy.fromstring(fid.read(nByte),readType).byteswapped()
  else:
    a1 = numpy.fromstring(fid.read(nByte),readType)

#=== mem_type
  if readType!=mem_type:
    a1 = a1.astype(mem_type)

#=== return
  return a1

#=======================
def fwrite(fid, num, myArray, write_type=None, byteswap=0):
#=======================

#=== figure out write_type
  mem_type = myArray.typecode()
  if write_type==0 or write_type==None:
    write_type = mem_type

#=== figure out total size of myArray
  temp1 = myArray.shape
  totalSize = 1
  for i in temp1: 
    totalSize = totalSize * i

  if num==None:
    num = totalSize

#=== if not all elements of myArray are written
  if num!=totalSize:
    myArray1 = numpy.ravel(myArray)[:num]
  else:
    myArray1 = myArray

#=== write
  if mem_type!=write_type:
    if byteswap: 
      fid.write( myArray1.astype(write_type).byteswapped().tostring() )
    else:
      fid.write( myArray1.astype(write_type).tostring() )
  else:
    if byteswap: 
      fid.write( myArray1.byteswapped().tostring() )
    else:
      fid.write( myArray1.tostring() )
