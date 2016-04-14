#!/usr/bin/env python

""" Fortran unformatted binary read and write module

Fortran unformatted binary data format is not specified by the standard, i.e.
it depends on the platforms and compilers.

Since USUAL compilers on x86 platforms (such as g77) use 4byte header and
footer (record marker) for unformatted I/O, fort_read/fort_write assume 4byte
by default. However, you can use the differnct size of record marker by
properly specifying the `length' argument of both functions. If the size is 0,
fort_read and fort_write are identical to fread and fwrite of scipy.io module,
respecively.

Note that this module use scipy package.

Copyright (c) 2007, Takanobu Amano <amano@eps.s.u-tokyo.ac.jp>
All rights reserved.

This software is provided under the new BSD license.
See for detail : http://www.opensource.org/licenses/bsd-license.php

$Id: fio.py 228 2007-03-15 10:02:39Z amano $
"""
import os
import sys
import struct
import scipy
from scipy import io

# byteorder
if sys.byteorder == 'little':
    _swapchar = '>' # big endian
elif sys.byteorder == 'big':
    _swapchar = '<' # little endian
else:
    raise "invalied sys.byteorder : %s" % (sys.byteorder)

# for compiler checking
_FC = 'g77 -m32'
_LENGTH = 4
_FRSRC = 'fior.f'
_FWSRC = 'fiow.f'
_read_code = """\
      program read_binary
      implicit none
      integer i
      real*8 x(4)

      open(unit=9, status='old', form='unformatted', file='fio.dat')
      read(9) x
      close(9)

      ! show contents
      do i = 1, 4
        write(*,*) 'read from fortran : ', i, x(i)
      end do

      end
"""
_write_code = """\
      program write_binary
      implicit none
      integer i
      real*8 x(4)

      ! show contents
      do i = 1, 4
        x(i) = 0.1 * i
        write(*,*) 'write from fortran : ', i, x(i)
      end do

      open(unit=9, status='replace', form='unformatted', file='fio.dat')
      write(9) x
      close(9)

      end
"""

def fort_read(fid, num, read_type, mem_type=None, byteswap=0, length=4):
    """read fortran unformatted binary data

    the meaning of argumetn is the same as scipy.io.fread except for `length',
    which is is the size (in bytes) of header/footer.
    """
    if mem_type == None:
        mem_type = read_type
    if length == 0: # without header/footer
        return io.fread(fid, num, read_type, mem_type, byteswap)
    # with header/footer
    fid.read(length) # header
    result = io.fread(fid, num, read_type, mem_type, byteswap)
    fid.read(length) # footer
    return result

def fort_write(fid, num, myarray, write_type=None, byteswap=0, length=4):
    """write fortran unformatted binary data

    the meaning of argumetn is the same as scipy.io.fread except for `length',
    which is is the size (in bytes) of header/footer.
    """
    if write_type == None:
        write_type = myarray.dtype.char
    if length == 0:
        return io.fwrite(fid, num, myarray, write_type, byteswap)
    elif length == 4:
        if byteswap:
            marker = _swapchar + 'i'
        else:
            marker = 'i'
    elif length == 8:
        if byteswap:
            marker = _swapchar + 'l'
        else:
            marker = 'l'
    else:
        raise ValueError, "length argument should be either 0, 4, or 8"

    # write
    nbyte = num * myarray.itemsize
    fid.write(struct.pack(marker, nbyte)) # header
    result = io.fwrite(fid, num, myarray, write_type, byteswap)
    fid.write(struct.pack(marker, nbyte)) # footer
    return result

def check_compiler(fc, length):
    """check compiler behaviour

    fc     : compiler command
    length : size of header/footer (reocrd marker)
    """
    # create read_binary program
    fr = open(_FRSRC, 'w')
    fr.write(_read_code)
    fr.close()
    robj = _FRSRC.split('.')[0]
    os.system("%s %s -o %s" % (fc, _FRSRC, robj))
    # create write_binary program
    fw = open(_FWSRC, 'w')
    fw.write(_write_code)
    fw.close()
    wobj = _FWSRC.split('.')[0]
    os.system("%s %s -o %s" % (fc, _FWSRC, wobj))

    ## write by fortran program
    os.system("./%s" % (wobj))
    ## read from python
    f = open('fio.dat', 'rb')
    x = fort_read(f, 4, 'd', length=length)

    ## show contents
    for i in range(4):
        print "show contents from python : %2d %5.2f" % (i, x[i])
    f.close()

    ## write by python
    f = open('fio.dat', 'wb')
    fort_write(f, 4, x, length=length)
    f.close()
    ## read from fortran
    os.system("./%s" % (robj))

    # remove
    [os.remove(f) for f in [_FRSRC, robj, _FWSRC, wobj, 'fio.dat']]

if __name__ == "__main__":
    "test fortran unformatted I/O"
    #
    # default compilers behaviour
    #
    # * x86 linux
    # most of compilers will use 4byte record marker
    #
    # * x86_64 linux
    # g77-3.4       : 8byte
    # gfortran-4.1  : 8byte
    # gfortran>4.2  : 4byte (?)
    # ifort-9.1     : 4byte (for large file, automatically use 8byte)
    #
    print """Fortran unformatted I/O test code

    Since the format depends on the platform and compilers,
    you should check it from the output of this program.
    Please check if the content of array is the same between fortran
    and python program.
    ***** information *****
    the compiler used here           : %-20s
    the size of record marker        : %-2d
    """ % (_FC, _LENGTH)
    # check
    check_compiler(_FC, _LENGTH)
