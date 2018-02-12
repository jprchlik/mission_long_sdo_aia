from PIL import Image, ImageFont, ImageDraw
import subprocess
from aia_mkmovie.make_movie import create_movie
from glob import glob
from sunpy.cm import cm
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
import sunpy.map
from datetime import datetime,timedelta as dt
from astropy import units as u
from astropy.coordinates import SkyCoord
import os,sys
import pywt
#get a median filter for wavelet transformation
from scipy.signal import medfilt



#create directories without erroring if they already exist c
def create_dir(dirs):
    """
    Create directories used by the program.
    """
    try:
        os.mkdir(dirs)
    except OSError:
        sys.stdout.write('{0} Already Exists'.format(dirs))


def make_images(f):
    global wav,img_scale,wx,wy,h0,w0,sdir

    #try to make the image. If it fails just move on
    try:
        #width of ind. image
        img_w = w0/len(wav)
        if wx > wy:
           img_wx = wx/len(wav)
           img_wy = wy
        else:
           img_wy = wy/len(wav)
           img_wx = wx


        wavelet = False


        #read all images into sunpy maps
        img = sunpy.map.Map(*f)
  
        #dictionary of images
        img_dict = {} 
        scale = {}
        scale_list = []
        #create new image
        new_img = Image.new('RGB',(w0,h0))

        #image size of subwindow
        sub_img_size = (img_w,h0)
        #put parameters in a series of dictionaries
        for j,i in enumerate(img):

            #create image position based on index
            if j == 0:
                px,py = 0,0
            elif j == 2:
                px,py = 0,1024
            elif j == 1:
                px,py = 1024,0
            elif j == 3:
                px,py = 1024,1024

            #color mapping
            icmap = img_scale[wav[j]][0]
            ivmin = img_scale[wav[j]][1]
            ivmax = img_scale[wav[j]][2]
            #keep list of scale
            scale = [i.scale[0].value,i.scale[1].value]

            #set up for wavelet analysis
            wav_img = i.data
            f_img = wav_img

            #do wavelet analysis
            if wavelet:
                d_size = 15*4+1
                #get median filter
                n_back = medfilt(wav_img,kernel_size=d_size)
                #subtract median filter
                img_sub = wav_img-n_back
                
                #Use Biorthogonal Wavelet
                wavelet = 'bior2.6'
                
                #use 6 levels
                n_lev = 6
                o_wav = pywt.swt2(img_sub, wavelet, level=n_lev )
                #only use the first 4
                f_img = pywt.iswt2(o_wav[0:4],wavelet)
                #Add  wavelet back into image
                f_img = f_img+wav_img
                
                #remove zero values
                f_img[f_img < 0.] = 0.

            #do normalization
            b_img = (np.arcsinh(f_img) - ivmin) / (ivmax - ivmin)
            #format image with color scale
            img_n = np.array(b_img)
            #if greater than 1 set to 0.99
            img_n[img_n > 0.99] = 0.99
            #print(wav[j],img_n.max(),f_img.max(),np.percentile(f_img,[5.,99]))
            #img_n[img_n < -.2] = 0.99
            img_n[img_n < 0.] = 0.
            img_n = icmap(img_n)
            img_n = np.uint8(img_n*255)

            if img_n.max() > 255:
                print(img_n.max())
                print(img[0].date.strftime('%Y%m%d_%H%M%S'))
               

            

            img_n = Image.fromarray(img_n)
            img_dict[wav[j]] = img_n

            #resize image to img0 scale
            img_n = img_n.resize((1024,1024))


            #default image size
            old_size = img_n.size
            horizontal_padding = (1024 - old_size[0]) / 2
            vertical_padding   = (1024 - old_size[1]) / 2
            temp_img = img_n.crop(
                (
                    -horizontal_padding,
                    -vertical_padding,
                    old_size[0] + horizontal_padding,
                    old_size[1] + vertical_padding
                )
            )




            #Add image to array of images
            new_img.paste(temp_img,(px,py))
        

        #output file
        outfi = sdir+'/working/panel_{0}'.format(img[0].date.strftime('%Y%m%d_%H%M%S'))+'.png'

        #set scale for plotting 
        #observed time 
        obs_time = img[0].date

        
        #write on text 
        w_text = '{0:%Y/%m/%d %H:%M:%S} '.format(obs_time)
        #add wavelengths
        for i in wav: w_text += str(int(i))+u'\u212B/'
        #remove final /
        w_text = w_text[:-1]
        
        
        #Add text of datetime to image
        draw = ImageDraw.Draw(new_img)
        draw.text((10,10),w_text,(255,255,255),font=font)
   
        #save image
        new_img.save(outfi) 
    except: 
        pass

#make sure the input wavelength matches the searched wavelength
def check_wavelength(fil,wav,archive,xrt=False):
    """
    Check wavelength and cadence of images
    """

    new_fil = []
    fil_dat = []
    if str(wav) == 'xrt':
        datefmt = 'L1_XRT%Y%m%d_%H%M%S'
    else:
        datefmt = 'AIA%Y%m%d_%H%M'

    #retrieve file wavelength and observation time
    for i in fil:
        #try assuming download format of jsoc files
       # try:
        if str(wav) != 'xrt':
            date = datetime.strptime(i.strip(archive).split('.')[0][:-5],datefmt)
            wave = i.strip(archive).split('_')[2].strip('.fits')
            if int(wave) == int(wav):
                new_fil.append(i)
                fil_dat.append(date)
        else:
            date = datetime.strptime(i.strip(archive+'xrt/').split('.')[0],datefmt)
            new_fil.append(i)
            fil_dat.append(date)

    #convert fil_dat to numpy time array 
    timelist = np.array(fil_dat)

    final_list = [] #list of index to keep for image creation
    for p in real_cad: #loop over all cadence values to find best array values
            k = np.abs(timelist-p)
            rindex, = np.where(k == k.min()) #get the nonzero array index value
            final_list.append(new_fil[rindex[0]])

    return final_list

#retrieve desired cadence from file list
def des_cad(start,end,delta):
    """Create an array from start to end with desired cadence"""
    curr = start
    while curr < end:
        yield curr
        curr += delta

#input wavelengths and cadence
wav = ['0094','0193','0211','0131']
cad = dt(days=30)

#Time range to observe
start = datetime(2010,5,22,1,0,0)
end   = datetime.utcnow()
#start = datetime(2017,9,10,21,50,0)
#end   = datetime(2017,9,10,21,59,0)
span = (end-start).total_seconds() #seconds to run the movie over
sdir = start.strftime('%Y%m%d') 
#get set cadence to observe 
real_cad = [result for result in des_cad(start,end,cad)]

img_scale = {'0094':[cm.sdoaia94  ,np.arcsinh(1.),np.arcsinh(150.)],
                  '0131':[cm.sdoaia131 ,np.arcsinh(5.),np.arcsinh(3500.)],
                  '0171':[cm.sdoaia171 ,np.arcsinh(30.),np.arcsinh(15000.)],
                  '0193':[cm.sdoaia193 ,np.arcsinh(100.),np.arcsinh(4500.)],
                  '0211':[cm.sdoaia211 ,np.arcsinh(10.),np.arcsinh(4000.)],
                  '0304':[cm.sdoaia304 ,np.arcsinh(2.),np.arcsinh(300.)],
                  '0335':[cm.sdoaia335 ,np.arcsinh(1.),np.arcsinh(100.)],
                  '1600':[cm.sdoaia1600,np.arcsinh(20.),np.arcsinh(500.)],
                  '1700':[cm.sdoaia1700,np.arcsinh(200.),np.arcsinh(4000.)],
                  'xrt':[cm.hinodexrt,np.arcsinh(50.),np.arcsinh(20000.)]}


#create a directory which will contain the raw png files
#sdir = stard+eday.date().strftime('%Y%m%d')
#creating a subdirectory to extra step is not need
dirlist = [sdir,sdir+'/working',sdir+'/working/symlinks',sdir+'/final']
for i in dirlist: create_dir(i)

#location of file archive
arch = 'sdo_archive/'


fits_files_temp = []
#loop over all wavelengths in array
for i in wav:
    try:
        fits_files = glob(arch+'*_'+i+'.fits')
    except:
        fits_files = glob(sdir+'/xrt/*fits')
    #make sure the wavelength header agrees with found value
    fits_files = check_wavelength(fits_files,i,arch)
    fits_files_temp.append(fits_files)
#transpose list array
fits_files = map(list,zip(*fits_files_temp))


#use font
#font = ImageFont.truetype("/Library/Fonts/Times New Roman Bold.ttf", 56)
font = ImageFont.truetype("/Library/Fonts/Arial Unicode.ttf", 56)
font = ImageFont.truetype("/Volumes/Pegasus/jprchlik/anaconda2/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/DejaVuSans-Bold.ttf", 56)


#number of processors
n_proc = 8

#image and movie height and width
wy = 2048
wx = 2048
#Use the longer side to automatically set width
if wy > wx:
   w0 = wy
   h0 = wx
else:
   h0 = wy
   w0 = wx

#testing
#fits_files = fits_files[170:190]
#Create fits images
pool = Pool(processes=n_proc)
outp = pool.map(make_images,fits_files)
pool.close()


#create movie
outf = sdir+'/final/'+sdir+'.mp4'


#frame rate
frate = 12
#bit rate
bitr = int(w0*h0*frate)/1000
#create movie object
mo =create_movie(odir = sdir+'/final/',pdir = sdir+'/working/', ext = 'png', w0 = int(w0), h0=int(h0),frate=frate,outmov=outf)
#run movie object
mo.gather_files()

command = "/usr/local/bin/ffmpeg -y -f image2 -r {1:1d} -i {4}/seq%{5}d.png -an -pix_fmt 'yuv420p' "
command = command+"-vcodec libx264 -level 41 -crf 18.0 -b {6:1d}k -r {1:1d} -bufsize {6:1d}k -maxrate {6:1d}k "
command = command+"-g 8 -coder 1 -profile main -preset faster -qdiff 4 -qcomp 0.7 -directpred 3 -flags "
command = command+"+loop+mv4 -cmp +chroma -partitions +parti4x4+partp8x8+partb8x8 -subq {1:1d} -me_range 16 "
command = command+"-keyint_min 1 -sc_threshold 40 -i_qfactor 0.71 -rc_eq 'blurCplx^(1-qComp)' -s '{2:4.0f}x{3:4.0f}' "
command = command+"-b_strategy 1 -bidir_refine 1 -refs 6 -deblockalpha 0 "
command = command+"-deblockbeta 0 -trellis 1 -x264opts keyint={1:1d}:min-keyint=1:bframes=1 -threads 8 {0}"
command = command.format(outf,frate,w0,h0,sdir+'/working/symlinks',mo.lengs,bitr)

run = subprocess.call(['/bin/tcsh','-c',command])