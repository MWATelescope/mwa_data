# Generate calibration models from SUMMS and VLA
# Sources we need: 3C444, HerA, 3C353, VirA, PicA, PKS2356-61, PKS2153-69

# Natasha Hurley-Walker 10/07/2013
# Todos: get more info from image headers instead of hardcoding

import re

rmtables('*.im')

VLSSr_sources=['3C444', 'HerA', '3C353', 'VirA','HydA']
SUMSS_sources=['PKS2356-61', 'PKS2153-69']
VLA333_sources=['PicA']

sources=VLSSr_sources+SUMSS_sources+VLA333_sources

spec={}
pos={}

# Spectra from bright_sources.vot because I'm too tired to do otherwise!
spec['3C444']=-0.95
spec['HerA']=-1.0
spec['3C353']=-0.85
spec['VirA']=-0.86
spec['HydA']=-0.83
spec['PicA']=-0.97
spec['PKS2356-61']=-0.85
spec['PKS2153-69']=-0.85

# Not used in the script, but could be scripted to grab from postage stamp server
pos['3C444']='22 14 25.752 -17 01 36.29'
pos['HerA']='16 51 08.2 +04 59 33'
pos['3C353']='17 20 28.147 -00 58 47.12'
pos['VirA']='12 30 49.42338 +12 23 28.0439'
pos['HydA']='09 18 05.651 -12 05 43.99'
pos['PKS2356-61']='23 59 04.365 -60 54 59.41'
pos['PKS2153-69']='21 57 05.98061 -69 41 23.6855'
pos['PicA']='05 19 49.735 -45 46 43.70'

# VLSSr
f0=74 #MHz
psf_rad='75arcsec' # beam in stamp headers
pix_size='5arcsec' # chosen on the postage stamp server
pix_area=qa.convert(pix_size,radians)['value']*qa.convert(pix_size,radians)['value']
psf_vol=pix_area/(1.1331*qa.convert(psf_rad,radians)['value']*qa.convert(psf_rad,radians)['value'])

for source in VLSSr_sources:
  print source
  exp=str(psf_vol)+'*IM0/((150/'+str(f0)+')^('+str(spec[source])+'))'
  model='templates/'+source+'_VLSS.fits'
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outname)
  exp=(str(spec[source]))+'*(IM0/IM0)'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outspec)

#SUMSS
f0=843 #MHz
#psf_rad='43arcsec' #1999AJ....117.1578B
#pix_size='5arcsec' # chosen on the postage stamp server
#pix_area=qa.convert(pix_size,radians)['value']*qa.convert(pix_size,radians)['value']
#psf_vol=pix_area/(1.1331*qa.convert(psf_rad,radians)['value']*qa.convert(psf_rad,radians)['value'])

# Seems SUMSS images are already in Jy/pixel?

for source in SUMSS_sources:
  print source
#  exp=str(psf_vol)+'*IM0/((150/'+str(f0)+')^('+str(spec[source])+'))'
  exp='IM0/((150/'+str(f0)+')^('+str(spec[source])+'))'
  model='templates/'+source+'_SUMSS.fits'
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile='new.im')
  exp=(str(spec[source]))+'*(IM0/IM0)'
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
  exp=str(psf_vol)+'*IM0/((150/'+str(f0)+')^('+str(spec[source])+'))'
  model='templates/'+source+'_VLA333.fits'
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outname)
  exp=(str(spec[source]))+'*(IM0/IM0)'
  immath(imagename=[model],mode='evalexpr',expr=exp,outfile=outspec)

# Update headers and export all images
for source in sources:
  outname=source+'.im'
  outspec=source+'_spec_index.im'
  imhead(imagename=outname,mode='put',hdkey='crval3',hdvalue='150MHz')
  imhead(imagename=outname,mode='put',hdkey='cdelt3',hdvalue='30.72MHz')
  imhead(imagename=outname,mode='put',hdkey='bunit',hdvalue='Jy/pixel')
  imhead(imagename=outname,mode='put',hdkey='crval4',hdvalue='I')
  imhead(imagename=outspec,mode='put',hdkey='crval3',hdvalue='150MHz')
  imhead(imagename=outspec,mode='put',hdkey='cdelt3',hdvalue='30.72MHz')
  imhead(imagename=outspec,mode='put',hdkey='bunit',hdvalue='Jy/pixel')
  imhead(imagename=outspec,mode='put',hdkey='crval4',hdvalue='I')
  fitsimage=re.sub('.im','.fits',outname)
  exportfits(imagename=outname,fitsimage=fitsimage)
  fitsimage=re.sub('.im','.fits',outspec)
  exportfits(imagename=outspec,fitsimage=fitsimage)
