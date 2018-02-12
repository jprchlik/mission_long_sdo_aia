import urllib
from datetime import datetime,timedelta
import os


#download files from archive for each wavelength
for i in d_wav:
   #format wavelength
   w_fil = w_fmt.format(i)
   urllib.urlretrieve(syn_arch+s_dir+w_fil,self.cmsdir+self.basedir+s_dir.split('/')[-1]+w_fil) 
   urllib.urlretrieve(syn_arch+e_dir+w_fil,self.cmsdir+self.basedir+e_dir.split('/')[-1]+w_fil) 




#location of syntopics
syn_arch = 'http://jsoc.stanford.edu/data/aia/synoptic/'


d_wav = [94, 193, 211, 131]

#wavelength format
w_fmt = '{0:04d}.fits'