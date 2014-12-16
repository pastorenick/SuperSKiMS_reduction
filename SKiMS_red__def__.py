import glob, pyraf, time, numpy
from scipy import interpolate

from fitsfitter_MOD import *
from flux_level import *
from flux_level2 import *
from plotFLUX import *

#Retrieve dictionary
lib_path = os.path.abspath('/Users/npastorello/Desktop/Galaxies/General_studies/')
sys.path.append(lib_path)
from galaxyParametersDictionary_v5 import *


###############
### CLASSES ###
###############


class mask():
#
  def __init__(self):
    self.name, self.RA_c, self.Dec_c = '', nan, nan # RA_c and Dec_c are the centre coordinates
    self.maskPA, self.PA = nan, nan #maskPA is the real PA of the mask in the sky, PA that in the header
    self.galaxy, self.gal_RA, self.gal_Dec = '', nan, nan
    self.gal_PA0, self.gal_Re, self.gal_ba = nan, nan, nan
    self.header = {}
    self.listSlits = []
  #
  def assignGalaxy(self, namegal):
    self.galaxy = namegal
    RA0, Dec0 = CentreCoordinates[namegal]
    self.gal_RA, self.gal_Dec = convAngCoord(RA0)[4]*15., convAngCoord(Dec0)[4] #in degrees
    self.gal_PA0, self.gal_Re, self.gal_ba = PA0[namegal], Reff[namegal]/3600., b_a[namegal]
  #
  def assignMaskPars(self, filename): #Takes the fits file name from which extract the header parameters
    inTmp = pyfits.open(filename)
    self.header = inTmp[0].header
    #
    self.name = self.header['OBJECT'] 
    self.RA_c, self.Dec_c = convAngCoord(self.header['RA'])[4]*15., convAngCoord(self.header['DEC'])[4] # in degrees
    #
    self.PA = (90.+self.header['SKYPA1'])
    deltaPA = mod(self.PA-self.gal_PA0, 360.) #Difference between mask alignment and galaxy PA0
    #
    angleNE = mod(90.-self.gal_PA0, 360)  #angle between galaxy major axis and East axis
    self.maskPA = -mod(angleNE+deltaPA, 360) #angle between the East-West axis and the mask alignment
  #
  def createSlits(self, listFiles):  #Creates a list of objects 'slit' with the proper coordinates and parameters
    #
    for ii in numpy.arange(len(listFiles)):
      #
      tmpName = listFiles[ii]
      hdu_tmp = pyfits.open(listFiles[ii])
      hdr_tmp = hdu_tmp[0].header
      #
      if (hdr_tmp['RA'] != '0.0') and (hdr_tmp['DEC'] != '0.0'):
      #
        dummyRA, dummyDec = hdr_tmp['RA'], hdr_tmp['DEC']
      #
        RAslit = convAngCoord(dummyRA)[4]*15.
        Decslit = convAngCoord(dummyDec)[4]
      else: # SERENDIPITY OBJECTS DON'T HAVE COORDINATES
        #looking for associated object for coordinates
        try:
          indexSlitAssociated = findPosAssociatedSlit(ii, listFiles)
          hdu_tmp2 = pyfits.open(listFiles[indexSlitAssociated])
          hdr_tmp2 = hdu_tmp2[0].header
          dummyRA, dummyDec = hdr_tmp2['RA'], hdr_tmp2['DEC']
          #
          RAslit = convAngCoord(dummyRA)[4]*15.
          Decslit = convAngCoord(dummyDec)[4]
        except:
          RAslit, Decslit = nan, nan
      #
      #
      distRA = RAslit-self.RA_c    #Distance from mask CentreCoordinates (in degrees)
      distDEC = Decslit-self.Dec_c #Distance from mask centre (in degrees)
      #
      angrot = self.maskPA*numpy.pi/180.
      realRA = (distRA*numpy.cos(angrot)-distDEC*numpy.sin(angrot))+self.RA_c   #Coordinates given the rotation of the mask
      realDEC = (distRA*numpy.sin(angrot)+distDEC*numpy.cos(angrot))+self.Dec_c #Coordinates given the rotation of the mask
      #
      slitTmp = slit()
      slitTmp.assignSlitPars([tmpName, realRA, realDEC, self.maskPA, 'tar'])
      self.listSlits.append(slitTmp)
  #
  def assignExtraParametersSlits(self, filename): #To assign length, centre position in the mask
    inputFits = pyfits.open(filename)
    slitParLayer = inputFits[11].data
    for ii in self.listSlits:
      for jj in slitParLayer:
      	if ii.ID == jj[9]:
          ii.slitlength = jj[5] #in arcsec
          ii.PArel = self.PA-jj[6]


class slit(mask):
#
  def __init__(self):
    self.name, self.ID, self.PA, self.Dec, self.RA, self.action = '', '', nan, 0, nan, ''
    self.wv, self.flux, self.ivar = numpy.array([]), numpy.array([]), numpy.array([])
    self.skyspec, self.totspec, = numpy.array([]), numpy.array([])
    self.skyindex, self.slitlength, self.PArel = nan, nan, nan #PArel is the PA relative to the mask
    #
    self.check, self.heliocentric_corr = False, nan
    #
    self.vel, self.sigma, self.h3, self.h4, self.chi2 = nan, nan, nan, nan, nan #from pPXF
    self.errvel, self.errsigma, self.errh3, self.errh4 = nan, nan, nan, nan #from pPXF
    self.ppxf_obj, self.wv_bestfit = None, numpy.array([])
#
  def assignSlitPars(self, listValues):
    self.name, self.RA, self.Dec, self.PA, self.action = listValues
    self.ID = findListNumber(self.name)
#
  def assignSpectrum(self, wv, totSpec, skySpec, subSpec, ivarSpec):
    self.wv, self.flux, self.ivar = wv, subSpec, ivarSpec
    self.skyspec, self.totspec = skySpec, totSpec


#################
### FUNCTIONS ###
#################


def plotSpec(fitsfile, ax=''):
  if ax == '':
    ax = subplot(111)
  tmpIn = pyfits.open(fitsfile)
  flux = tmpIn[0].data
  wv = numpy.arange(len(flux))*tmpIn[0].header['CDELT1'] + tmpIn[0].header['CRVAL1']
  ax.plot(wv, flux)

def findListNumber(stringName):
  for ii in numpy.arange(len(stringName)):
    if stringName[ii] == '.':
      return stringName[ii+1:ii+4]

def findPosAssociatedSlit(index,listFILES):
  # 
  # Look for slit number 
  #
  tmp = []
  for ii in listFILES:
    tmp.append(findListNumber(ii))
  #
  listSameNumber = numpy.nonzero(numpy.array(tmp) == findListNumber(listFILES[index]))
  for ii in listSameNumber[0]:
    if not('serendip' in listFILES[ii]):
      return ii


def weighted_std(values, weights):
  average = numpy.average(values, weights=weights)
  variance = numpy.average((values-average)**2, weights=weights)  # Fast and numerically precise
  return (numpy.sqrt(variance))

# Montecarlo run for pPXF
def runMC(pp, dv, velScale, start, nReal=100, quiet=True):
  import ppxf
  #
  ppMC = []
  # 
  for ii in numpy.arange(nReal):
    #
    if verbose:
      stdout.write("\r MC  %i percent " % (((ii)+1)*100./nReal))
      stdout.flush()
    #
    # Create fake spectrum
    numpy.random.seed(int(time.time()*1e6 % 1e5))
    #
    fakeGal = numpy.random.normal(pp.bestfit, pp.noise)
    #
    # Run pPXF over fake spectrum and  store pPXF object
    ppMC.append(ppxf.ppxf(pp.star, fakeGal, pp.noise, velScale, start, 
                    plot=False, moments=pp.moments, degree=pp.degree, vsyst=dv, quiet=quiet,
                    bias=pp.bias, oversample=pp.oversample)
                )
    #
  if verbose:
    stdout.write("\n")
  return ppMC

def extractErrorKin(ppxfObjList):
  #
  errorKin = [] 
  #
  for ii in numpy.arange(ppxfObjList[0].moments):
    tmpList, tmpListWeights = [], []
    for jj in ppxfObjList:
      tmpList.append(jj.sol[ii])
      tmpListWeights.append(jj.error[ii])  #The weights are the errors on the single kin measurements
    #
    errorKin.append(weighted_std(tmpList, 1./numpy.array(tmpListWeights)))
  #
  return numpy.array(errorKin)



##### SKiMS reduction ######

def runSKiMSreduction(Deimos_mask):
  ### EXTRACTING SPECTRA FROM spec1d FILES
  # The output files (should be) already dispersion corrected
  # Also the .asc files are created
  print "EXTRACTING SPECTRA FROM FITS FILES"
  for ii in glob.glob('*.?.*fits'): os.remove(ii)
  #
  if not(os.path.exists('./login.cl')):
    os.system('mkiraf')
  #
  runFitsfitter(verbose=True)
  print "\nDONE"
  #  
  ### EXTRACTING COORDINATES FROM fits FILES AND SAVE "DATA.DAT" file
  #listFILES=asciidata.open('../R_sky/D_listaSPEC.txt') 
  #
  print "EXTRACTING COORDINATES FROM FITS FILES"
  listFILES = glob.glob('*.r.sky.fits')
  Deimos_mask.createSlits(listFILES)
  #
  listNAM, listRA, listDEC, action = [], [], [], []
  for ii in Deimos_mask.listSlits:
  	listNAM.append(ii.name)
  	listRA.append(ii.RA)
  	listDEC.append(ii.Dec)
  	action.append(ii.action)
  #
  outputTab = numpy.transpose(numpy.array([listNAM, action, listRA, listDEC]))
  numpy.savetxt('./DATA.DAT',outputTab,fmt='%s',header='file\taction\tRAdeg\tDecdeg',newline='\n',delimiter='\t')
  print "\nDONE"
  #
  ### MEASURING SKY INDEX ON ALL THE SPECTRA
  print "MEASURING SKY INDEX IN ALL THE SPECTRA"
  fluxLevel(Deimos_mask, inputFile='DATA.DAT', outputFile='DATA2.DAT', setLines='standard')
  print "\nDONE"
  #  
  ### SELECTION SKY SPECTRA
  print "SELECTION SKY SPECTRA"
  plotFlux(Deimos_mask, inputFile='DATA2.DAT')
  print "\nDONE"
  #  
  ### SKY_SUBTRACT
  print "CREATING MASTER SKY SPECTRUM"
  #
  ### Cleaning dir
  if len(glob.glob('./sky*fits')) > 0:
    for ii in glob.glob('./sky*fits'):
       os.remove(ii)
  #  
  if len(glob.glob('./sky*asc')) > 0:
    for ii in glob.glob('./sky*asc'):
       os.remove(ii)
  #
  filename = './DATA2.DAT'
  flag = True
  while flag:  #Check at least 1 sky-dominated spectrum exists
    listFILES = asciidata.open(filename)
    listNAMsky, listNAM, action = [], [], []
    l1, l2, l3, l4 = [], [], [], []
    #
    for ii in range(0,len(listFILES[0])):
      if (listFILES[1][ii] == 'sky'):
        listNAMsky.append("./"+listFILES[0][ii])
      listNAM.append(listFILES[0][ii])
      l1.append(listFILES[4][ii])
      l2.append(listFILES[5][ii])
      l3.append(listFILES[6][ii])
      l4.append(listFILES[7][ii])
    #
    if (len(listNAMsky) > 0):
      flag = False
    else:
      raw_input("No SKY spectra selected in '"+filename+"'. Please fix it and press a button.\n")
  #  
  if (os.path.isfile("./temp_sky.lis")): os.system("rm ./temp_sky.lis")
  if (os.path.isfile("./temp_gal.lis")): os.system("rm ./temp_gal.lis")
  #
  numpy.savetxt("./temp_sky.lis",listNAMsky,fmt='%s')
  numpy.savetxt("./temp_gal.lis",listNAM,fmt='%s')
  #  
  ############### Deleting old files
  #
  if os.path.exists('./temp.lis'): os.remove('./temp.lis')
  #
  if os.path.exists('./temp_sky.fits'): os.remove('./temp_sky.fits')
  #
  if os.path.exists('./norm_sky.fits'): os.remove('./norm_sky.fits')
  #
  if os.path.exists('./temp_bands.lis'): os.remove('./temp_bands.lis')
  #
  if os.path.exists('./temp_sky_bands.lis'): os.remove('./temp_sky_bands.lis')
  #
  ############### Make sky estimates
  #  
  fiddle = 1.00
  #
  # Create sky master spec
  pyraf.iraf.onedspec.scombine("@./temp_sky.lis",
                                  "sky.fits", noutput="", logfile="STDOUT", apertures="", group="apertures",
                                  combine="sum", reject="none", first='no', w1='INDEF', w2='INDEF', dw='INDEF',
                                  nw='INDEF', scale="none", zero="none", weight="none", sample="",
                                  lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=1, mclip='yes',
                                  lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1,
                                  pclip=-0.5, grow=0, blank=0.)
  #
  #
  # Convert sky master spec into txt table
  pyraf.iraf.onedspec.wspec("./sky.fits","./sky.asc", header='no', wformat="")
  #
  # plt.ion()
  # fig = figure()
  # clf()
  # ax = subplot(111)
  # plotSpec('sky.fits',ax=ax)
  #
  #####################
  #
  ############   Measure fluxes in sky estimates
  #
  os.system("ls ./sky.fits > ./temp.lis")
  #
  run_flux_level2()
  #  
  pyraf.iraf.columns("./temp_levels.lis", 5, outroot="./col.")
  os.system('mv ./col.3 ./temp_sky_bands.lis')
  pyraf.iraf.delete("./col.*", 'yes', verify='no', default_acti='yes', allversions='yes', subfiles='yes')
  #
  #####################
  #
  #
  ############  Normalise sky estimate
  #
  flux_level = numpy.loadtxt("./temp_sky_bands.lis")
  pyraf.iraf.imarith("./sky.fits", "/", flux_level, "./norm_sky.fits",
          title="", divzero=0., hparams="", pixtype="", calctype="", verbose='no',
          noact='no')
  #
  os.remove('./temp.lis')
  # 
  #####################
  #
  ############ Measure fluxes in target frames and remove sky  
  ### 1st mask
  #
  os.system('cp ./temp_gal.lis ./temp.lis')      
  #        cat ./temp.lis
  #
  run_flux_level2()
  #
  pyraf.iraf.columns("./temp_levels.lis", 5, outroot="./col.")
  os.system("mv ./col.3 ./temp_bands.lis")
  pyraf.iraf.delete("./col.*", 'yes', verify='no', default_acti='yes', allversions='yes', subfiles='yes')
  #
  listFiles = numpy.loadtxt("./temp_gal.lis",dtype='S50')
  flux_level = numpy.loadtxt("./temp_bands.lis",dtype='f50')
  #
  if len(glob.glob('sub_*.txt'))>0: os.system('rm sub_*.txt')
  #  
  tmpInput = pyfits.open('norm_sky.fits')
  fluxSkyMaster = tmpInput[0].data
  #
  print "DONE\n"
  #  
  #### CUT SKY SPECTRUM TO MATCH WITH THE SCIENCE SPECTRUM
  #
  print "SUBTRACTING MASTER SKY FROM ALL THE SPECTRA"
  #
  # INTERPOLATE SKY SPECTRUM
  flux_mastersky_norm = tmpInput[0].data
  wv_mastersky_norm = numpy.arange(len(tmpInput[0].data))*tmpInput[0].header['CDELT1']+tmpInput[0].header['CRVAL1']
  f_mastersky_norm = interpolate.interp1d(wv_mastersky_norm, flux_mastersky_norm)
  max_Sky_wv, min_Sky_wv = numpy.max(wv_mastersky_norm), numpy.min(wv_mastersky_norm)
  #
  # EXTRAPOLATE SKY SPECTRUM ON SAME WAVELENGTH ARRAY AS THE SCIENCE SPECTRUM
  for ii in numpy.arange(len(Deimos_mask.listSlits)):
    #
    stdout.write("\r Subtracting spectrum  %i / %i " % ((ii)+1, len(Deimos_mask.listSlits)))
    stdout.flush()
    #
    inputSpecFile = pyfits.open(Deimos_mask.listSlits[ii].name)
    flux_spectrum = inputSpecFile[0].data
    inputIvarFile = pyfits.open(Deimos_mask.listSlits[ii].name[:-5]+'ivar.fits')
    ivar_spectrum = inputIvarFile[0].data
    wv_spectrum = numpy.arange(len(inputSpecFile[0].data))*inputSpecFile[0].header['CDELT1']+inputSpecFile[0].header['CRVAL1']
    #
    selWV = numpy.nonzero((wv_spectrum <= max_Sky_wv) & (wv_spectrum >= min_Sky_wv))
    flux_spectrum_cut, wv_spectrum_cut, ivar_spectrum_cut = flux_spectrum[selWV], wv_spectrum[selWV], ivar_spectrum[selWV]
    #
    flux_sky_norm_extrapolated = f_mastersky_norm(wv_spectrum_cut)
    #
    # Multiplying sky per flux level
    flux_sky_extrapolated = flux_sky_norm_extrapolated*fiddle*flux_level[ii]
    # Subtracting sky from spectrum
    sub_spectrum = flux_spectrum_cut-flux_sky_extrapolated
    #
    # save in objects
    Deimos_mask.listSlits[ii].assignSpectrum(wv_spectrum_cut, flux_spectrum_cut, flux_sky_extrapolated, 
    								sub_spectrum, ivar_spectrum_cut)
    #
    outputTab = numpy.transpose(numpy.array([wv_spectrum_cut, sub_spectrum, ivar_spectrum_cut]))
    numpy.savetxt('sub_'+listFiles[ii][:-4]+'txt', outputTab, fmt='%.10f', delimiter='\t', newline='\n',header='WV\tFLUX\tIVAR')
    #   
  print "\nDONE"




##### Check Spectra ######
def runCheckSpectra(Deimos_mask, galname='NGC1023'):
  #
  plt.ion()
  #
  RA, Dec = [], []
  for ii in Deimos_mask.listSlits:
  	RA.append(ii.RA)
  	Dec.append(ii.Dec)
  #
  RA, Dec = numpy.array(RA), numpy.array(Dec)
  selNotNan = numpy.nonzero(~(isnan(RA)))
  minRA, maxRA = numpy.min(RA[selNotNan]), numpy.max(RA[selNotNan])
  minDec, maxDec = numpy.min(Dec[selNotNan]), numpy.max(Dec[selNotNan])
  #
  fig = figure(0, figsize=(12.6,6))
  ax2 = subplot(122)
  ax2.scatter(Deimos_mask.gal_RA, Deimos_mask.gal_Dec, marker='x', color='r')
  radiuses = numpy.array([1,2,3])
  ells = [Ellipse(xy=[Deimos_mask.gal_RA, Deimos_mask.gal_Dec], width=(2.*ii*Deimos_mask.gal_Re/numpy.sqrt(Deimos_mask.gal_ba)), 
                  height=(2.*ii*Deimos_mask.gal_Re*numpy.sqrt(Deimos_mask.gal_ba)), angle=90-Deimos_mask.gal_PA0, 
                  edgecolor = 'k', facecolor = 'none', fill = False, linestyle = 'dashed') for ii in radiuses]
  #
  for ee in ells: ax2.add_artist(ee)
  #
  ax2.set_xlim([maxRA, minRA])
  ax2.set_ylim([minDec, maxDec])
  ax1 = subplot(121)
  # 
  listNAME, listRA, listDEC, listCheck = [], [], [], []
  for sel in Deimos_mask.listSlits:
    cla()
    ax1 = subplot(121)
    wv, flux, ivar = sel.wv, sel.flux, sel.ivar
    ax1.plot(wv, flux)
    ax1.set_xlim([8500,8800])
    ax1.set_ylim([numpy.median(flux)/3,numpy.median(flux)*2])
    #
    answer = raw_input()
    if answer == '':
      cc = 'b'
      listCheck.append(True)
      sel.check = True
    else:
      cc = 'k'
      listCheck.append(False)
      sel.check = False
    ax2.scatter(sel.RA, sel.Dec, color=cc, marker='.')
    draw()
    #
  return True




# To save plot of distro 
def drawDistro(Deimos_mask, galname='NGC1023', outputname='outputSlit.pdf'):
  #
  RA, Dec, check, name = [], [], [], []
  for ii in Deimos_mask.listSlits:
  	RA.append(ii.RA)
  	Dec.append(ii.Dec)
  	check.append(ii.check)
  	name.append(ii.name)
  #
  fig = figure(num=1, figsize=(6,5))
  ax = subplot(111)
  #
  ax.scatter(Deimos_mask.gal_RA, Deimos_mask.gal_Dec, marker='x', color='r')
  radiuses = numpy.array([1,2,3])
  ells = [Ellipse(xy=[Deimos_mask.gal_RA, Deimos_mask.gal_Dec], width=(2.*ii*Deimos_mask.gal_Re/numpy.sqrt(Deimos_mask.gal_ba)), 
                  height=(2.*ii*Deimos_mask.gal_Re*numpy.sqrt(Deimos_mask.gal_ba)), angle=90-Deimos_mask.gal_PA0, 
                  edgecolor = 'k', facecolor = 'none', fill = False, linestyle = 'dashed') for ii in radiuses]
  #
  for ee in ells: 
    ax.add_artist(ee)
  #
  ax.set_xlim([40.12999999,40.0399999])
  #
  ax.scatter(RA, Dec, color='k')
  selGood = numpy.nonzero(numpy.array(check))
  ax.scatter(numpy.array(RA)[selGood], numpy.array(Dec)[selGood], color='g')
  #
  ax.scatter(Deimos_mask.RA_c, Deimos_mask.Dec_c, color='b', marker='x')
  #
  for ii in Deimos_mask.listSlits:
    ax.annotate(ii.name[12:-11], xy=(ii.RA, ii.Dec))
  #
  draw()
  savefig(outputname, bbox_inches='tight')
  return True


# Finding Heliocentric corrections for the mask
def findHeliocentricCorr(Deimos_mask):
  year, month, day = Deimos_mask.header['DATE-OBS'][0:4], Deimos_mask.header['DATE-OBS'][5:7], Deimos_mask.header['DATE-OBS'][8:10]
  RAt, Dect, UTCt = Deimos_mask.header['RA'], Deimos_mask.header['Dec'], Deimos_mask.header['UTC']
  RAh, Decd, UTh = convAngCoord(RAt)[4], convAngCoord(Dect)[4], convAngCoord(UTCt)[4]
  #
  pyraf.iraf.rv()
  pyraf.iraf.rv.rvcorrect(epoch='2000', observa='keck', vsun='20', ra_vsun='18', dec_vsu='30', epoch_vsu='1900', 
                        year=year, month=month, day=day, ut=UTh, ra=RAh, dec=Decd)
  heliocentric_corr = float(raw_input('write heliocentric correction: '))
  #
  print "\n\n\n"
  return heliocentric_corr



def run_pPXF(Deimos_mask, templates='deimos', cutRegion = [8450, 8750]):	
  #
  lib_path = os.path.abspath('/Users/npastorello/Desktop/CaTPros/pPXF_files')
  sys.path.append(lib_path)
  import ppxf, ppxf_util
  #
  heliocentric_corr = findHeliocentricCorr(Deimos_mask)
  #
  if templates == 'deimos':	#Standard stars observed with DEIMOS
    #
    # Saving templates
    listTemplates = glob.glob('/Users/npastorello/Desktop/CaTPros/DEIMOS_StellarTemplate/*.fits')  #DEIMOS TEMPLATES
  elif templates == 'cenarro':
  	listTemplates = glob.glob('/Volumes/G2/exLustre/REDUCTION/JacobReductionSkims/stellartemplates/cat_spec_fits/scan*.fits')  
  #
  tmp = pyfits.open(listTemplates[0])
  lenTemp = len(tmp[0].data)
  wvTemp = numpy.arange(lenTemp)*tmp[0].header['CDELT1']+tmp[0].header['CRVAL1']
  #selCut = numpy.nonzero((wvTemp < 8750.) & (wvTemp > 8450.))
  #wvTemp_cut, fluxTemp_cut = wvTemp[selCut], tmp[0].data[selCut]
  wvTemp_cut, fluxTemp_cut = wvTemp, tmp[0].data
  #
  sspTempl, logLam2, velscale = ppxf_util.log_rebin([wvTemp_cut[0],wvTemp_cut[-1]], fluxTemp_cut)
  templatesArray = numpy.empty((len(sspTempl),len(listTemplates)))
  #
  counter = 0
  for ii in listTemplates:
    tmp = pyfits.open(ii)
    fluxTemp_cut = tmp[0].data#[selCut]
    sspNew, logLam2, velscale = ppxf_util.log_rebin([wvTemp_cut[0],wvTemp_cut[-1]], fluxTemp_cut, velscale=velscale)
    templatesArray[:,counter] = sspNew
    counter += 1
  #
  if not(os.path.exists('FitSpectra')): os.mkdir('FitSpectra')
  #
  if not(os.path.exists('Output')): os.mkdir('Output')
  #
  listFail = []
  for ii in Deimos_mask.listSlits:
    try:
      dicOutput = {}
      #
      raSlit, decSlit, nameSlit = ii.RA, ii.Dec, ii.name
      ii.heliocentric_corr = heliocentric_corr
      #
      # Reading spectrum
      wv, flux, ivar = ii.wv, ii.flux, ii.ivar
      #
      # Applying heliocentric correction with IRAF
      #
      wv = wv+(heliocentric_corr*wv/c)
      #
      #rescaling wv onto same wv scale as templates
      newWV = numpy.arange(numpy.min(wv), numpy.max(wv), tmp[0].header['CDELT1'])
      functFlux = scipy.interpolate.interp1d(wv, flux)
      newFlux = functFlux(newWV)
      functIvar = scipy.interpolate.interp1d(wv, ivar)
      newIvar = functIvar(newWV)
      #
      # Cutting spectrum in CaT region
      #
      selCut = numpy.nonzero((newWV < cutRegion[1]) & (newWV > cutRegion[0]))
      wv_cut, flux_cut, ivar_cut = newWV[selCut], newFlux[selCut], newIvar[selCut]
      #
      #The variance of the original sky file is underestimated (it doesn't take in account the sky-subtraction)  
      #
      lamRange1 = numpy.array([wv_cut[0], wv_cut[-1]])
      galaxy, logLam1, velscale = ppxf_util.log_rebin(lamRange1, flux_cut)
      #
      # Normalizing spectra and noise
      #
      fac = numpy.median(galaxy) 
      galaxy /= fac # Normalize spectrum to avoid numerical issues
      variance, logLam1, velscale = ppxf_util.log_rebin(lamRange1, 1./ivar_cut)
      noise = numpy.sqrt(variance)/fac
      #
      if templates == 'cenarro': 
      	print 'Not yet implemented.'
      	return False
      	#resTempl = 1.5 #22km/s 1.5 AA (FWHM), sigma~22 km/s, R~5700
      elif templates == 'deimos': # NOT CONVOLVING BY THE RESOLUTION, BECAUSE THE INSTRUMENT IS THE SAME 
      #
        dummy = True
      #
      dv = (logLam2[0]-logLam1[0])*c # km/s
      vel = vel0[Deimos_mask.galaxy]
      goodPixels = ppxf_util.determine_goodpixels(logLam1, [wvTemp[0],wvTemp[-1]], vel)
      #
      start = [vel, 180.] # (km/s), starting guess for [V,sigma]
      #
      fig = figure(0, figsize=(6,5))
      clf()
      pp2 = ppxf.ppxf(templatesArray, galaxy, noise, velscale, start,
                     goodpixels=goodPixels, 
                     plot=False, moments=2, 
                     vsyst=dv, degree=4, oversample=False, quiet=True)
      #
      pp4 = ppxf.ppxf(templatesArray, galaxy, noise, velscale, pp2.sol[0:2],
                     goodpixels=goodPixels, 
                     plot=True, moments=4, 
                     vsyst=dv, degree=4, oversample=False)
      #
      ii.vel, ii.sigma, ii.h3, ii.h4 = pp4.sol
      ii.ppxf_obj = pp4
      ii.wv_bestfit = numpy.e**loglam1
      #
      pp4_MC = runMC(pp4, dv, velscale, pp2.sol[0:2], nReal=100, quiet=True)
      ii.errvel, ii.errsigma, ii.errh3, ii.errh4 = extractErrorKin(pp4_MC)
      #
      suptitle(ii.name)
      savefig('FitSpectra/'+ii.name+'.pdf',bbox_inches='tight')
      #
      #
      # inputTab = numpy.loadtxt('DATA2.DAT', dtype={'names': ('filename', 'action', 'RA', 'Dec', 'dummy1', 'dummy2', 'dummy3', 'dummy4'),
      #                    'formats': ('S50', 'S3', 'f30', 'f30', 'f30', 'f30', 'f30', 'f30')})
      # #
      # dicOutput['nameSlit'] = nameSlit
      # dicOutput['RASlit'], dicOutput['DECSlit'] = raSlit, decSlit
      # dicOutput['pPXF_obj'] = pp4
      # #
      # dicOutput['vel'] = v_slit
      # dicOutput['sigma'] = sigma_slit
      # dicOutput['h3'] = h3_slit
      # dicOutput['h4'] = h4_slit
      # #
      # dicOutput['evel'] = v_err
      # dicOutput['esigma'] = sigma_err
      # dicOutput['eh3'] = h3_err
      # dicOutput['eh4'] = h4_err
      # #
      # #    
      # fileOut = open('Output/'+nameSlit+'.dat', 'wb')
      # pickle.dump(dicOutput, fileOut)
      # fileOut.close()
    except:
      listFail.append(ii.name)
  return listFail


def run_saveKinPlot(Deimos_mask, outputFile='kinplot.pdf', 
					xlimits=[300,-300]):	#In arcsec 
  ### Create 2d kinematic plot
  fig = figure(num=1, figsize=(12.6, 10))
  clf()
  plt.subplots_adjust(hspace = .000, wspace = .000)
  ax_vel = subplot(2,2,1, aspect='equal')
  ax_sigma = subplot(2,2,2, aspect='equal')
  ax_h3 = subplot(2,2,3, aspect='equal')
  ax_h4 = subplot(2,2,4, aspect='equal')
  #
  RA0, Dec0 = CentreCoordinates[Deimos_mask.galaxy]
  RA0deg, Dec0deg = Deimos_mask.gal_RA, Deimos_mask.gal_Dec    
  pa0, ba0, reff = Deimos_mask.gal_PA0, Deimos_mask.gal_ba, Deimos_mask.gal_Re
  #
  ax_vel.scatter(0, 0, marker='x', color='r')
  ax_sigma.scatter(0, 0, marker='x', color='r')
  ax_h3.scatter(0, 0, marker='x', color='r')
  ax_h4.scatter(0, 0, marker='x', color='r')
  #
  RA, Dec, vel, sigma, h3, h4, check = [], [], [], [], [], [], []
  for ii in Deimos_mask.listSlits:
    RA.append((ii.RA-RA0deg)*3600.) #in arcsec
    Dec.append((ii.Dec-Dec0deg)*3600.)
    vel.append(ii.vel-vel0[Deimos_mask.galaxy])
    sigma.append(ii.sigma)
    h3.append(ii.h3)
    h4.append(ii.h4)
    check.append(ii.check)
  #
  RA, Dec = numpy.array(RA), numpy.array(Dec)
  vel, sigma = numpy.array(vel), numpy.array(sigma)
  h3, h4 = numpy.array(h3), numpy.array(h4)
  #selgal = numpy.nonzero((RA <250) & (RA > -300) & ~(numpy.isnan(RA)))
  selgal = numpy.nonzero(numpy.array(check) == True)
  vpoints = ax_vel.scatter(RA[selgal], Dec[selgal], c=vel[selgal])
  posAx = ax_vel.get_position()
  posAx.y0, posAx.y1 = posAx.y1, posAx.y0+0.02
  axCB = fig.add_axes(posAx)
  CB = colorbar(vpoints, orientation='horizontal', cax=axCB)
  CB.ax.xaxis.set_ticks_position('top')
  #
  spoints = ax_sigma.scatter(RA[selgal], Dec[selgal], c=numpy.log10(sigma[selgal]))
  posAx = ax_sigma.get_position()
  posAx.y0, posAx.y1 = posAx.y1, posAx.y0+0.02
  axCB = fig.add_axes(posAx)
  CB = colorbar(spoints, orientation='horizontal', cax=axCB)
  CB.ax.xaxis.set_ticks_position('top')
  #
  h3points = ax_h3.scatter(RA[selgal], Dec[selgal], c=h3[selgal])
  posAx = ax_h3.get_position()
  posAx.y1, posAx.y0 = posAx.y0 - 0.04, posAx.y0 - 0.06
  axCB = fig.add_axes(posAx)
  CB = colorbar(h3points, orientation='horizontal', cax=axCB)
  CB.ax.xaxis.set_ticks_position('bottom')
  #
  h4points = ax_h4.scatter(RA[selgal], Dec[selgal], c=h4[selgal])
  posAx = ax_h4.get_position()
  posAx.y1, posAx.y0 = posAx.y0 - 0.04, posAx.y0 - 0.06
  axCB = fig.add_axes(posAx)
  CB = colorbar(h4points, orientation='horizontal', cax=axCB)
  CB.ax.xaxis.set_ticks_position('bottom')
  #
  ax_vel.xaxis.set_ticklabels([])
  ax_vel.set_xlim(xlimits)
  ax_vel.annotate('Velocity', xy=(-200,40))
  ax_sigma.yaxis.set_ticklabels([])
  #tmpax = ax_sigma.twinx()
  #ax_sigma.xaxis.set_ticklabels([])
  #tmpax.set_ylim(ax_vel.set_ylim())
  ax_sigma.set_xlim(xlimits)
  ax_sigma.annotate('Velocity dispersion', xy=(-100,40))
  #
  ax_h3.set_xlim(xlimits)
  ax_h3.annotate('h3', xy=(-200,40))
  #
  ax_h4.yaxis.set_ticklabels([])
  #tmpax = ax_h4.twinx()
  #tmpax.set_ylim(ax_h3.set_ylim())
  ax_h4.set_xlim(xlimits)
  ax_h4.annotate('h4', xy=(-200,40))
  #
  radiuses = numpy.array([1,2,3])
  ells = [Ellipse(xy=[0, 0], width=(2.*ii*Deimos_mask.gal_Re/numpy.sqrt(Deimos_mask.gal_ba)), 
                  height=(2.*ii*Deimos_mask.gal_Re*numpy.sqrt(Deimos_mask.gal_ba)), angle=90-Deimos_mask.gal_PA0, 
                  edgecolor = 'k', facecolor = 'none', fill = False, linestyle = 'dashed') for ii in radiuses]
  for ee in ells: 
    ax_vel.add_artist(ee) 
  #
  savefig(outputFile)
  return True


    



































