import urllib
from datetime import datetime,timedelta
import os
import itertools


#retrieve desired cadence from file list
def des_cad(self,start,end,delta):
    """Create an array from start to end with desired cadence"""
    curr = start
    while curr < end:
        yield curr

#wrapper for download file for par. processing
def wrap_download_file(args):
    return download_file(*args)

#download files from archive for each wavelength
def download_file(time,wavl):
   global w_fmt,f_dir,b_dir
   #format wavelength
   w_fil = w_fmt.format(wavl)
   s_dir = f_dir.format(time)
   urllib.urlretrieve(syn_arch+s_dir+w_fil,b_dir+s_dir.split('/')[-1]+w_fil) 



#base local SDO archive directory
b_dir = 'sdo_archive/'

#location of syntopics
syn_arch = 'http://jsoc.stanford.edu/data/aia/synoptic/'

#Wavelength download
d_wav = [94, 193, 211, 131]

#wavelength format
w_fmt = '{0:04d}.fits'

#create directory path minus the wavelength
f_dir = '{0:%Y/%m/%d/H%H00/AIA%Y%m%d_%H%M_}'

#get starttime for observations
stime = datetime(2010,5,13,0,0,0)

#set endtime for observations
etime = datetime.utcnow()

#cadence to get observations
caden = timedelta(days=365)

#desired cadence for the observations
real_cad = [result for result in des_cad(self.start,self.end,timedelta(minutes=caden))]

#create an iterable combination of dates and wavelengths
inpt_itr = intertools.product(read_cad,d_wav)
