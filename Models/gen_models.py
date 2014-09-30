# Generate calibration models from SUMMS and VLA
# Sources we need: 3C444, HerA, 3C353, VirA, PicA, PKS2356-61, PKS2153-69

# Natasha Hurley-Walker 10/07/2013
# Todos: get more info from image headers instead of hardcoding

# Note that this requires CASA 4+ to run, since CASA3 can't open images without beams in the headers

import re, math, os, glob, pyfits

rmtables('*.im')

filelist = glob.glob("*.fits")
for f in filelist:
    os.remove(f)

VLSSr_sources=['3C444', 'HerA', '3C353', 'VirA','HydA','3C33', '3C161', 'PKS0442-28', 'Crab']
SUMSS_sources=['PKS2356-61', 'PKS2153-69', 'PKS0408-65', 'PKS0410-75']
VLA333_sources=['PicA']
SUMSS_deconv_sources=['CenA']

sources=VLSSr_sources+SUMSS_sources+VLA333_sources+SUMSS_deconv_sources

spec={}
pos={}

# Spectra from bright_sources.vot because I'm too tired to do otherwise!
spec['3C444']=-0.95
spec['HerA']=-1.0
spec['3C353']=-0.85
spec['3C33']=-0.85
spec['VirA']=-0.86
spec['HydA']=-0.83
spec['PicA']=-0.97
spec['PKS2356-61']=-0.85
spec['PKS2153-69']=-0.85
spec['PKS0408-65']=-0.85
spec['PKS0410-75']=-0.85
spec['3C161']=-0.85
spec['PKS0442-28']=-0.85
spec['CenA']=-0.497 
# CenA: in http://ned.ipac.caltech.edu/cgi-bin/ex_refcode?refcode=1981A%26AS...45..367K
# this basically refers to what it is in NED, so 80 MHz: 1544, 160 MHz: 1094
# derived spectral index: -0.497
spec['Crab']=-0.30 # Hurley-Walker 2009

# Not used in the script, but could be scripted to grab from postage stamp server
pos['3C444']='22 14 25.752 -17 01 36.29'
pos['HerA']='16 51 08.2 +04 59 33'
pos['3C353']='17 20 28.147 -00 58 47.12'
pos['VirA']='12 30 49.42338 +12 23 28.0439'
pos['HydA']='09 18 05.651 -12 05 43.99'
pos['PKS2356-61']='23 59 04.365 -60 54 59.41'
pos['PKS2153-69']='21 57 05.98061 -69 41 23.6855'
pos['PicA']='05 19 49.735 -45 46 43.70'
pos['PKS0408-65']=''
pos['PKS0410-75']=''
pos['3C33']=['01 08 52.854 +13 20 13.75']
pos['3C161']=['06 27 10.0960 -05 53 04.768']
pos['PKS0442-28']=['04 44 37.707 -28 09 54.41']

# RW's deconvolved CenA image
f0=843 #MHz

for source in SUMSS_deconv_sources:
  print source
  model='templates/'+source+'_SUMSSd.fits'

# Only include real sources
  exp='iif(IM0>=0.033,IM0/((150/'+str(f0)+')^('+str(spec[source])+')),0.0)'
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile='new.im')
  exp='iif(IM0>=0.033,'+(str(spec[source]))+'*(IM0/IM0),0.0)'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile='newspec.im')

  ia.open('new.im')
  ia.adddegaxes(outfile=outname,spectral=True,stokes='I')
  ia.close()
  ia.open('newspec.im')
  ia.adddegaxes(outfile=outspec,spectral=True,stokes='I')
  ia.close()
  rmtables('new*')

# VLSSr
f0=74 #MHz
psf_rad='75arcsec' # beam in stamp headers
pix_size='5arcsec' # chosen on the postage stamp server
pix_area=qa.convert(pix_size,radians)['value']*qa.convert(pix_size,radians)['value']
psf_vol=pix_area/(1.1331*qa.convert(psf_rad,radians)['value']*qa.convert(psf_rad,radians)['value'])

for source in VLSSr_sources:
  print source
# Only include real sources
  exp='iif(IM0>=2.0,'+str(psf_vol)+'*IM0*((150/'+str(f0)+')^('+str(spec[source])+')),0.0)'
  model='templates/'+source+'_VLSS.fits'
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outname)
  exp=(str(spec[source]))+'*(IM0/IM0)'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outspec)

#SUMSS
f0=843 #MHz

for source in SUMSS_sources:
  print source
  model='templates/'+source+'_SUMSS.fits'
  importfits(fitsimage=model,imagename='test.im')
  dec=imhead('test.im',mode='get',hdkey='crval2')['value']
# psf = 43"x43" csc|delta| according to 1999AJ....117.1578B
  cosecdec=abs(43.0/math.sin(dec))
  psf_bmaj='43arcsec'
  psf_bmin=str(cosecdec)+'arcsec'
  pix_size='5arcsec' # chosen on the postage stamp server
  pix_area=qa.convert(pix_size,radians)['value']*qa.convert(pix_size,radians)['value']
  psf_vol=pix_area/(1.1331*qa.convert(psf_bmaj,radians)['value']*qa.convert(psf_bmin,radians)['value'])
  rmtables('test.im')

# Only include real sources
  exp='iif(IM0>=0.25,'+str(psf_vol)+'*IM0/((150/'+str(f0)+')^('+str(spec[source])+')),0.0)'
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile='new.im')
  exp='iif(IM0>=0.25,'+(str(spec[source]))+'*(IM0/IM0),0.0)'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile='newspec.im')

  ia.open('new.im')
  ia.adddegaxes(outfile=outname,spectral=True,stokes='I')
  ia.close()
  ia.open('newspec.im')
  ia.adddegaxes(outfile=outspec,spectral=True,stokes='I')
  ia.close()
  rmtables('new*')

# VLA 333 MHz proper observation
f0=333 #MHz
psf_rad='30arcsec' # convolving beam in fits history
pix_size='1.25arcsec' # from the fits header
pix_area=qa.convert(pix_size,radians)['value']*qa.convert(pix_size,radians)['value']
psf_vol=pix_area/(1.1331*qa.convert(psf_rad,radians)['value']*qa.convert(psf_rad,radians)['value'])

for source in VLA333_sources:
  print source
#  exp=str(psf_vol)+'*IM0/((150/'+str(f0)+')^('+str(spec[source])+'))'
  exp=str(psf_vol)+'*IM0*((150/'+str(f0)+')^('+str(spec[source])+'))'
  model='templates/'+source+'_VLA333.fits'
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outname)
  exp=(str(spec[source]))+'*(IM0/IM0)'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outspec)

# Delete excess keywords that cause problems in later versions of pywcs
dead_keywords=['PC001001', 'PC002001', 'PC003001', 'PC004001', 'PC001002', 'PC002002', 'PC003002', 'PC004002', 'PC001003', 'PC002003', 'PC003003', 'PC004003', 'PC001004', 'PC002004', 'PC003004', 'PC004004']

# Update headers and export all images
for source in sources:
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  imhead(imagename=outname,mode='put',hdkey='crval3',hdvalue='150MHz')
  imhead(imagename=outname,mode='put',hdkey='cdelt3',hdvalue='30.72MHz')
  imhead(imagename=outname,mode='put',hdkey='bunit',hdvalue='Jy/pixel')
  imhead(imagename=outspec,mode='put',hdkey='crval3',hdvalue='150MHz')
  imhead(imagename=outspec,mode='put',hdkey='cdelt3',hdvalue='30.72MHz')
  imhead(imagename=outspec,mode='put',hdkey='bunit',hdvalue='Jy/pixel')
# 1s not Is in the latest version of CASA
  currentstokes=imhead(imagename=outname,mode='get',hdkey='crval4')
  if currentstokes['value']!=1.0:
    imhead(imagename=outname,mode='put',hdkey='crval4',hdvalue='I')
    imhead(imagename=outspec,mode='put',hdkey='crval4',hdvalue='I')
# Delete excess keywords that cause problems in later versions of pywcs
  fitsimage=re.sub('.im','.fits',outname)
  exportfits(imagename=outname,fitsimage=fitsimage)
  hdu_in=pyfits.open(fitsimage,mode='update')
  for fitskey in dead_keywords:
    del hdu_in[0].header[fitskey]
  hdu_in.flush()
  fitsimage=re.sub('.im','.fits',outspec)
  exportfits(imagename=outspec,fitsimage=fitsimage)
  hdu_in=pyfits.open(fitsimage,mode='update')
  for fitskey in dead_keywords:
    del hdu_in[0].header[fitskey]
