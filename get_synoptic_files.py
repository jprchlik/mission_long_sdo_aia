import urllib
from datetime import datetime,timedelta
import os
import itertools
from multiprocessing import Pool


#retrieve desired cadence from file list
def des_cad(start,end,delta):
    """Create an array from start to end with desired cadence"""
    curr = start
    while curr < end:
        yield curr
        curr += delta

#wrapper for download file for par. processing
def wrap_download_file(args):
    return download_file(*args)

#download files from archive for each wavelength
def download_file(time,wavl):
   global w_fmt,f_dir,b_dir
   #format wavelength
   w_fil = w_fmt.format(wavl)
   #format input time
   s_dir = f_dir.format(time)
   #create output file
   o_fil = b_dir+s_dir.split('/')[-1]+w_fil
   #file to download from archive
   d_fil = syn_arch+s_dir+w_fil

   #check if output file exists
   if os.path.isfile(o_fil) == False:
       #try to download file if fails continue on
       try:
           urllib.urlretrieve(d_fil,o_fil) 
       except:
           print("Cound not Download {0} from archive".format(d_fil))



#base local SDO archive directory
b_dir = 'sdo_archive/'

#try to make the local base direcory
try:
    os.makedirs(b_dir)
except OSError:
    pass

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
caden = timedelta(days=30)

#desired cadence for the observations
real_cad = [result for result in des_cad(stime,etime,caden)]

#create a list of combination of dates and wavelengths
inpt_itr = list(itertools.product(real_cad,d_wav))

#Download the files locally
pool = Pool(processes=8)
outp = pool.map(wrap_download_file,inpt_itr)
pool.close()
pool.join()
