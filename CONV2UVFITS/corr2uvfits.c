/* Program to convert the binary output format of the 32T on-site correlators
 * UV FITS format. Written by Randall Wayth. Feb, 2008.
 *
 * August 12, 2011 - precess u, v, w, data to J2000 frame (Alan Levine)
$Rev: 4135 $:     Revision of last commit
$Author: rwayth $:  Author of last commit
$Date: 2011-10-18 22:53:40 +0800 (Tue, 18 Oct 2011) $:    Date of last commit
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include <fitsio.h>
#include <assert.h>
#include "slalib.h"
#include "uvfits.h"
#include "convutils.h"


/* private function prototypes */
void printusage(const char *progname);
int readScan(FILE *fp_ac, FILE *fp_cc,int scan,Header *header, InpConfig *inps,uvdata *data);
void parse_cmdline(const int argc, char * const argv[], const char *optstring);
void initData(uvdata *data);
int applyHeader(Header *header, uvdata *data);
void checkInputs(Header *header,uvdata *data,InpConfig *inputs);

/* private global vars */
int debug=0,do_flag=0,lock_pointing=0;
int pol_index[4];
FILE *fpd=NULL;
char *stationfilename="antenna_locations.txt",*outfilename=NULL;
char *configfilename="instr_config.txt";
char *header_filename="header.txt";
char *crosscor_filename=NULL;
char *autocorr_filename=NULL;
char *flagfilename=NULL;
double arr_lat_rad=MWA_LAT*(M_PI/180.0),arr_lon_rad=MWA_LON*(M_PI/180.0),height=MWA_HGT;

/************************
************************/
int main(const int argc, char * const argv[]) {
  const char optstring[] = "vldS:a:c:o:I:H:A:F:f:";
  FILE *fpin_ac=NULL,*fpin_cc=NULL;
  int i,scan=0,res=0;
  Header header;
  InpConfig inputs;
  uvdata *data;
  array_data *arraydat;
  source_table *source;
  ant_table *antennas;
  uvWriteContext *iter=NULL;

  fpd=stderr;

  if(argc < 2) printusage(argv[0]);
  parse_cmdline(argc,argv,optstring);

  /* initialise some values for the UV data array and antennas*/
  data = calloc(1,sizeof(uvdata));
  source = calloc(1,sizeof(source_table));
  arraydat = calloc(1,sizeof(array_data));
  antennas = calloc(MAX_ANT,sizeof(ant_table));
  assert(antennas!=NULL && source!=NULL && arraydat!=NULL && data !=NULL);
  data->array    = arraydat;
  data->array->antennas = antennas;
  data->source   = source;
  initData(data);

  /* get the mapping of inputs to anntena numbers and polarisations */
  if ((res = readInputConfig(configfilename, &inputs)) != 0) {
      fprintf(stderr,"readInputConfig failed with code %d. exiting\n",res);
  }

  /* get the number of antennas and their locations relative to the centre of the array */
  if ((res = readArray(stationfilename, arr_lat_rad, antennas,&(arraydat->n_ant))) != 0) {
      fprintf(stderr,"readArray failed with code %d. exiting\n",res);
  }

  /* read the header/metadata  */
  res = readHeader(header_filename,&header);
  if (res != 0) {
    fprintf(stderr,"Error reading main header. exiting.\n");
    exit(1);
  }

  checkInputs(&header,data,&inputs);

  /* open input files */
  if (header.corr_type!='A' && (fpin_cc=fopen(crosscor_filename,"r"))==NULL) {
    fprintf(stderr,"cannot open cross correlation input file <%s>\n",crosscor_filename);
    exit(1);
  }
  if (header.corr_type!='C' && (fpin_ac=fopen(autocorr_filename,"r"))==NULL) {
    fprintf(stderr,"cannot open auto correlation input file <%s>\n",autocorr_filename);
    exit(1);
  }

  /* assign vals to output data structure from inputs */
  res = applyHeader(&header, data);

  /* populate antenna info */
  if (debug) fprintf(fpd,"there are %d antennas\n",arraydat->n_ant);
  for (i=0; i<arraydat->n_ant; i++){
    //sprintf(antennas[i].name,"ANT%03d",i+1);
    // FIXME: update this for correct pol type
    sprintf(antennas[i].pol_typeA,"X");
    sprintf(antennas[i].pol_typeB,"Y");
    antennas[i].pol_angleA = 0.0;
    antennas[i].pol_angleB = 90.0;
    antennas[i].pol_calA = 0.0;
    antennas[i].pol_calB = 0.0;
    antennas[i].mount_type = 0;
  }

  /* assign XYZ positions of the array for the site. */
  Geodetic2XYZ(arr_lat_rad,arr_lon_rad,height,&(arraydat->xyz_pos[0]),&(arraydat->xyz_pos[1]),&(arraydat->xyz_pos[2]));
  if (debug) fprintf(fpd,"converted array location to XYZ\n");

  /* open the iterator for writing. Need to have set the date in  */
  res = writeUVFITSiterator(outfilename, data, &iter);
  if (res!=0) {
    fprintf(stderr,"writeUVFITSiterator returned %d.\n",res);
    return res;
  }
  if (debug) fprintf(fpd,"created UVFITS writer iterator\n");

  /* read each scan, populating the data structure. */
  scan=0;
  while ((res = readScan(fpin_ac,fpin_cc,scan, &header, &inputs, data))==0) {
    int status;

    if (flagfilename != NULL) {
      if (debug) fprintf(fpd,"Applying flags...\n");
      res = applyFlagsFile(flagfilename,data);
      if(res!=0) {
        fprintf(stderr,"Problems in applyFlagsFile. exiting\n");
        exit(1);
      }
    }

    /* write this chunk of data to the output file */
    status = writeUVinstant(iter, data, data->date[0]-iter->jd_day_trunc,0);
    if (status) {
        fprintf(stderr,"ERROR: code %d when writing time instant %d\n",status,i);
        exit(EXIT_FAILURE);
    }
    scan++;
    if (scan ==header.n_scans) break;   // don't read more than is specified
  }
  if(res < 0) {
      fprintf(stderr,"Problems in readScan(). exiting\n");
      exit(1);
  }

  if (debug) fprintf(fpd,"Read %d time steps\n",scan);
  if (abs(header.n_scans-scan) > 4) {
    fprintf(stderr,"WARNING: expected to read %d scans, but actually read %d. Are you sure you have the correct number of freqs, inputs and timesteps? Carrying on and hoping for the best...\n",header.n_scans,scan);
  }
  else if (res > 0) {
      fprintf(stderr,"WARNING: Wanted to read %d time steps, but actually read %d. Carrying on...\n",header.n_scans, scan);
  }

  if (do_flag) {
    if (debug) fprintf(fpd,"Auto flagging...\n");
    res = autoFlag(data,5.0,AUTOFLAG_N_NEIGHBOURS,do_flag);
    if(res!=0) {
      fprintf(stderr,"Problems in autoflag. exiting\n");
      exit(1);
    }
  }

  /* finish up and write the antenna table */
  writeUVFITSfinalise(iter, data);

  if(debug) fprintf(fpd,"finished writing UVFITS file\n");
  if(fpin_ac !=NULL) fclose(fpin_ac);
  if(fpin_cc !=NULL) fclose(fpin_cc);

  freeUVFITSdata(data);

  return 0;
}


/***************************
 ***************************/
int readScan(FILE *fp_ac, FILE *fp_cc,int scan_count, Header *header, InpConfig *inps, uvdata *uvdata) {

  static int bl_ind_lookup[MAX_ANT][MAX_ANT], init=0;
  static float vis_weight=1.0;
  static double date_zero=0.0;  // time of zeroth scan

  array_data *array;    /*convenience pointer */
  double ha=0,mjd,lmst,ra_app,dec_app;
  double ant_u[MAX_ANT],ant_v[MAX_ANT],ant_w[MAX_ANT]; //u,v,w for each antenna, in meters
  double u,v,w,cable_delay=0;
  double x,y,z, xprec, yprec, zprec, rmatpr[3][3], rmattr[3][3];
  double lmst2000, ha2000, newarrlat, ant_u_ep, ant_v_ep, ant_w_ep;
  double ra_aber, dec_aber;
  int i,inp1,inp2,ant1,ant2,pol1,pol2,bl_index=0,visindex,pol_ind,chan_ind,n_read;
  int temp,baseline_reverse,scan=0;
  float *visdata=NULL;

  array = uvdata->array;

  /* allocate space to read binary correlation data. Size is complex float * n_channels */
  visdata = calloc(2*uvdata->n_freq,sizeof(float));
  if (visdata==NULL) {
      fprintf(stderr,"ERROR: no malloc for visdata array in readScan\n");
      exit(1);
  }

  if (!init) {
    /* count the total number of antennas actually present in the data */
    if (debug) fprintf(fpd,"Init readscan.\n");

    /* make a lookup table for which baseline corresponds to a correlation product */
    if(header->corr_type=='A') {
        /* autocorrelations only */
        for(ant1=0; ant1 < uvdata->array->n_ant; ant1++) {
          if (checkAntennaPresent(inps,ant1) == 0) continue;
          if(debug) fprintf(fpd,"AUTO: bl %d is for ant %d\n",bl_index,ant1);
          bl_ind_lookup[ant1][ant1] = bl_index++;
      }
    }
    else if(header->corr_type=='C') {
        /* this for cross correlations only */
        for (ant1=0; ant1 < uvdata->array->n_ant-1; ant1++) {
          if (checkAntennaPresent(inps,ant1) == 0) continue;
            for(ant2=ant1+1; ant2 < uvdata->array->n_ant; ant2++) {
              if (checkAntennaPresent(inps,ant2) == 0) continue;
              if(debug) fprintf(fpd,"CROSS: bl %d is for ants %d-%d\n",bl_index,ant1,ant2);
              bl_ind_lookup[ant1][ant2] = bl_index++;
            }
        }
    }
    else {
        /* this for auto and cross correlations */
        for (ant1=0; ant1 < uvdata->array->n_ant; ant1++) {
          if (checkAntennaPresent(inps,ant1) == 0) continue;
            for(ant2=ant1; ant2 < uvdata->array->n_ant; ant2++) {
              if (checkAntennaPresent(inps,ant2) == 0) continue;
              if(debug) fprintf(fpd,"BOTH: bl %d is for ants %d-%d\n",bl_index,ant1,ant2);
              bl_ind_lookup[ant1][ant2] = bl_index++;
            }
        }
    }

    /* increase size of arrays for the new scan */
    // date and n_baselines should already be allocated
    assert(uvdata->date != NULL);
    assert(uvdata->n_baselines[0] > 0);
    uvdata->n_vis=1;
    uvdata->visdata=calloc(1,sizeof(double *));
    uvdata->weightdata=calloc(1,sizeof(double *));
    uvdata->u=calloc(1,sizeof(double *));
    uvdata->v=calloc(1,sizeof(double *));
    uvdata->w=calloc(1,sizeof(double *));
    uvdata->baseline = calloc(1,sizeof(float *));

    /* make space for the actual visibilities and weights */
    if (debug) fprintf(fpd,"readScan: callocing array of %d floats\n",uvdata->n_baselines[scan]*uvdata->n_freq*uvdata->n_pol*2);
    uvdata->visdata[scan]    = calloc(uvdata->n_baselines[scan]*uvdata->n_freq*uvdata->n_pol*2,sizeof(float));
    uvdata->weightdata[scan] = calloc(uvdata->n_baselines[scan]*uvdata->n_freq*uvdata->n_pol  ,sizeof(float));
    uvdata->u[scan] = calloc(uvdata->n_baselines[scan],sizeof(double));
    uvdata->v[scan] = calloc(uvdata->n_baselines[scan],sizeof(double));
    uvdata->w[scan] = calloc(uvdata->n_baselines[scan],sizeof(double));
    uvdata->baseline[scan] = calloc(uvdata->n_baselines[scan],sizeof(float));
    if(uvdata->visdata[scan]==NULL || uvdata->weightdata[scan]==NULL || uvdata->visdata[scan]==NULL
       || uvdata->visdata[scan]==NULL || uvdata->visdata[scan]==NULL || uvdata->baseline[scan]==NULL) {
      fprintf(stderr,"readScan: no malloc for BIG arrays\n");
      exit(1);
    }
    /* set a weight for the visibilities based on integration time */
    if(header->integration_time > 0.0) vis_weight = header->integration_time;

    date_zero = uvdata->date[0];    // this is already initialised in applyHeader

    init=1;
  }
  
  /* set time of scan. Note that 1/2 scan time offset already accounted for in date[0]. */
  if (scan_count > 0) uvdata->date[scan] = date_zero + scan_count*header->integration_time/86400.0;

  /* set default ha/dec from header, if HA was specified. Otherwise, it will be calculated below */
  dec_app = header->dec_degs*(M_PI/180.0);
  mjd = uvdata->date[scan] - 2400000.5;  // get Modified Julian date of scan.
  if (lock_pointing==0) {   // special case for RTS output which wants a phase centre fixed at an az/el, not ra/dec
    ha = (header->ha_hrs_start+(scan_count+0.5)*header->integration_time/3600.0*1.00274)*(M_PI/12.0);
  } else {
    ha = (header->ha_hrs_start)*(M_PI/12.0);
    mjd = date_zero - 2400000.5;  // get Modified Julian date of scan.
  }

  lmst = slaRanorm(slaGmst(mjd) + arr_lon_rad);  // local mean sidereal time, given array location

  /* convert mean RA/DEC of phase center to apparent for current observing time. This applies precession,
     nutation, annual abberation. */
  slaMap(header->ra_hrs*(M_PI/12.0), header->dec_degs*(M_PI/180.0), 0.0, 0.0, 0.0, 0.0, 2000.0, mjd, &ra_app, &dec_app);
  if (debug) fprintf(fpd,"Precessed apparent coords (radian): RA: %g, DEC: %g\n",ra_app,dec_app);
  /* calc apparent HA of phase center, normalise to be between 0 and 2*pi */
  ha = slaRanorm(lmst - ra_app);

  /* I think this is correct - it does the calculations in the frame with current epoch 
  * and equinox, nutation, and aberrated star positions, i.e., the apparent geocentric
  * frame of epoch. (AML)
  */

  if(debug) fprintf(fpd,"scan %d. lmst: %g (radian). HA (calculated): %g (radian)\n",scan_count, lmst,ha);

  /* calc el,az of desired phase centre for debugging */
  if (debug) {
      double az,el;
      slaDe2h(ha,dec_app,arr_lat_rad,&az,&el);
      fprintf(fpd,"Phase cent ha/dec: %g,%g. az,el: %g,%g\n",ha,dec_app,az,el);
  }

  /* Compute the apparent direction of the phase center in the J2000 coordinate system */
  aber_radec_rad(2000.0,mjd,header->ra_hrs*(M_PI/12.0),header->dec_degs*(M_PI/180.0), &ra_aber,&dec_aber);

  /* Below, the routines "slaPrecl" and "slaPreces" do only a precession correction,
   * i.e, they do NOT do corrections for aberration or nutation.
   *
   * We want to go from apparent coordinates at the observation epoch
   * to J2000 coordinates which do not have the nutation or aberration corrections
   * (and since the frame is J2000 no precession correction is needed).
   */

  // slaPrecl(slaEpj(mjd),2000.0,rmatpr);  /* 2000.0 = epoch of J2000 */
  slaPrenut(2000.0,mjd,rmattr);
  mat_transpose(rmattr,rmatpr);
  /* rmatpr undoes precession and nutation corrections */
  ha_dec_j2000(rmatpr,lmst,arr_lat_rad,ra_aber,dec_aber,&ha2000,&newarrlat,&lmst2000);

  if (debug) {
    fprintf(fpd,"Dec, dec_app, newarrlat (radians): %f %f %f\n",
        header->dec_degs*(M_PI/180.0),dec_app,newarrlat);
    fprintf(fpd,"lmst, lmst2000 (radians): %f %f\n",lmst,lmst2000);
    fprintf(fpd,"ha, ha2000 (radians): %f %f\n",ha,ha2000);
  }
  /* calc u,v,w at phase center and reference for all antennas relative to center of array */
  for(i=0; i<array->n_ant; i++) {
    // double x,y,z;   /* moved to front of this function (Aug. 12, 2011) */
      x = array->antennas[i].xyz_pos[0];
      y = array->antennas[i].xyz_pos[1];
      z = array->antennas[i].xyz_pos[2];
      /* value of lmst at current epoch - will be changed to effective value in J2000 system 
       * 
       * To do this, need to precess "ra, dec" (in quotes on purpose) of array center
       * from value at current epoch 
       */
      precXYZ(rmatpr,x,y,z,lmst,&xprec,&yprec,&zprec,lmst2000);
      calcUVW(ha,dec_app,x,y,z,&ant_u_ep,&ant_v_ep,&ant_w_ep);
      calcUVW(ha2000,dec_aber,xprec,yprec,zprec,ant_u+i,ant_v+i,ant_w+i);
      if (debug) {
        /* The w value should be the same in either reference frame. */
          fprintf(fpd,"Ant: %d, u,v,w: %g,%g,%g.\n",i,ant_u[i],ant_v[i],ant_w[i]);
          fprintf(fpd,"Ant at epoch: %d, u,v,w: %g,%g,%g.\n",i,ant_u_ep,ant_v_ep,ant_w_ep);
      }
  }
 
  for(inp1=0; inp1 < header->n_inputs ; inp1++) {
    for(inp2=inp1; inp2 < header->n_inputs ; inp2++) {

        /* decode the inputs into antennas and pols */
        baseline_reverse=0;
        ant1 = inps->ant_index[inp1];
        ant2 = inps->ant_index[inp2];
        pol1 = inps->pol_index[inp1];
        pol2 = inps->pol_index[inp2];
        /* UVFITS by convention expects the index of ant2 to be greater than ant1, so swap if necessary */
        if (ant1>ant2) {
            temp=ant1;
            ant1=ant2;
            ant2=temp;
            temp=pol1;
            pol1=pol2;
            pol2=temp;
            baseline_reverse=1;
        }
        pol_ind = decodePolIndex(pol1, pol2);
        bl_index = bl_ind_lookup[ant1][ant2];
        
        /* cable delay: the goal is to *correct* for differential cable lengths. The inputs include a delta (offset)
           of cable length relative to some ideal length. (positive = longer than ideal)
           Call the dot product of the baseline (ant2-ant1) and look direction 'phi'.
           Then if ant1 has more delay than ant2, then this is like having phi be positive where
           the visibility is V = Iexp(-j*2*pi*phi)
           Hence we want to add the difference ant2-ant1 (in wavelengths) to phi to correct for the length difference.
         */
        cable_delay = (inps->cable_len_delta[inp2] - inps->cable_len_delta[inp1]);
        
        /* only process the appropriate correlations */
        if (header->corr_type=='A' && inp1!=inp2) continue;
        if (header->corr_type=='C' && inp1==inp2) continue;

        /* There is now one block of channels for each correlation product, read it. */
        if (inp1 != inp2) {
            /* read a block of cross-correlations */
            n_read = fread(visdata,sizeof(float)*2,uvdata->n_freq,fp_cc);
        }
        else {
            /* read a block of auto-correlations */
            n_read = fread(visdata,sizeof(float),uvdata->n_freq,fp_ac);
        }
        if (n_read != uvdata->n_freq) {
            fprintf(stderr,"EOF: inps %d,%d. expected to read %d channels, only got %d\n",inp1, inp2,uvdata->n_freq,n_read);
            return 1;
        }
 
         /* throw away cross correlations from different pols on the same antenna if we only want cross products */
         /* we do this here to make sure that the data is read and the file advanced on to data we want */
         if (header->corr_type=='C' && ant1==ant2) continue;
 
        /* calc u,v,w for this baseline in meters */
        u=v=w=0.0;
        if(ant1 != ant2) {
            u = ant_u[ant1] - ant_u[ant2];
            v = ant_v[ant1] - ant_v[ant2];
            w = ant_w[ant1] - ant_w[ant2];
        }

        /* populate the baseline info. Antenna numbers start at 1 in UVFITS.  */
        uvdata->baseline[scan][bl_index] = 0; // default: no baseline. useful to catch bugs.
        EncodeBaseline(ant1+1, ant2+1, uvdata->baseline[scan]+bl_index);
        /* arcane units of UVFITS require u,v,w in nanoseconds */
        uvdata->u[scan][bl_index] = u/VLIGHT;
        uvdata->v[scan][bl_index] = v/VLIGHT;
        uvdata->w[scan][bl_index] = w/VLIGHT;

        if (debug) {
            fprintf(fpd,"doing inps %d,%d. ants: %d,%d pols: %d,%d, polind: %d, bl_ind: %d, w (m): %g, delay (m): %g, blrev: %d\n",
                    inp1,inp2,ant1,ant2,pol1,pol2,pol_ind,bl_index,w,cable_delay,baseline_reverse);
        }
        
        /* if not correcting for geometry, don't apply w */
        if(!header->geom_correct) {
            w = 0.0;
        }

        /* populate the visibility arrays */
        for(chan_ind=0; chan_ind<uvdata->n_freq; chan_ind++) {
            double freq,lambda;
            float complex vis,phase=1.0;

            /* calc wavelen for this channel in meters. header freqs are in MHz*/
            freq = (header->cent_freq + (header->invert_freq? -1.0:1.0)*(chan_ind - uvdata->n_freq/2.0)/uvdata->n_freq*header->bandwidth);
            lambda = (VLIGHT/1e6)/freq;
            phase = cexp(I*(-2.0*M_PI)*(w+cable_delay*(baseline_reverse? -1.0:1.0))/lambda);
            visindex = bl_index*uvdata->n_pol*uvdata->n_freq + chan_ind*uvdata->n_pol + pol_ind;
            vis = visdata[chan_ind*2] + I*(header->conjugate ? -visdata[chan_ind*2+1]: visdata[chan_ind*2+1]);

            if(debug && chan_ind==uvdata->n_freq/2) {
                fprintf(fpd,"Chan %d, w: %g (wavelen), vis: %g,%g. ph: %g,%g. rot vis: %g,%g\n",chan_ind,w/lambda,creal(vis),cimag(vis),creal(phase),cimag(phase),creal(vis*phase),cimag(vis*phase));
            }

            if (baseline_reverse) vis = conj(vis);
            vis *= phase;

            if (inp1 != inp2) {
                /* cross correlation, use imaginary and real */
                uvdata->visdata[scan][visindex*2   ] = crealf(vis);
                uvdata->visdata[scan][visindex*2 +1] = cimagf(vis);
            }
            else {
                /* auto correlation, set imag to zero */
                uvdata->visdata[scan][visindex*2   ] = visdata[chan_ind];
                uvdata->visdata[scan][visindex*2 +1] = 0.0;
            }
            uvdata->weightdata[scan][visindex] = vis_weight;
            // apply input-based flags if necessary
            if ( (inps->inpFlag[inp1] || inps->inpFlag[inp2]) && vis_weight > 0) {
                uvdata->weightdata[scan][visindex] = -vis_weight;
            }
        }
    }

  }
  if (visdata != NULL) free(visdata);
  return 0;
}


/****************************
*****************************/
void parse_cmdline(const int argc,char * const argv[], const char *optstring) {
    int result=0;
    char arrayloc[80],*lon,*lat;

    arrayloc[0]='\0';

    while ( (result = getopt(argc, argv, optstring)) != -1 ) {
        switch (result) {
          case 'S': stationfilename = optarg;
            break;
          case 'o': outfilename = optarg;
            break;
          case 'a': autocorr_filename = optarg;
            break;
          case 'c': crosscor_filename = optarg;
            break;
          case 'd': debug = 1;
            fprintf(fpd,"Debugging on...\n");
            break;
          case 'I': configfilename = optarg;
            break;
          case 'H': header_filename = optarg;
            break;
          case 'f': do_flag=atoi(optarg);
            break;
          case 'F': flagfilename=optarg;
            break;
          case 'l': lock_pointing=1;
            fprintf(fpd,"Locking phase center to initial HA/DEC\n");
            break;
          case 'A':
            strncpy(arrayloc,optarg,80);
            break;
          case 'v':
            fprintf(stdout,"corr2uvfits revision $Rev: 4135 $\n");
            exit(1);
            break;
          default:
              fprintf(stderr,"unknown option: %c\n",result);
              printusage(argv[0]);
        }
    }

    /* convert array lon/lat */
    if(arrayloc[0]!='\0') {
        lon = arrayloc;
        lat = strpbrk(arrayloc,",");
        if (lat ==NULL) {
            fprintf(stderr,"Cannot find comma separator in lon/lat. Typo?\n");
            printusage(argv[0]);
        }
        /* terminate string for lon, then offset for lat */
        *lat = '\0';
        lat++;
        arr_lat_rad = atof(lat);
        arr_lon_rad = atof(lon);
        fprintf(fpd,"User specified array lon,lat: %g, %g (degs)\n",arr_lon_rad,arr_lat_rad);
        arr_lat_rad *= (M_PI/180.0); // convert to radian
        arr_lon_rad *= (M_PI/180.0);
    }

    /* do some sanity checks */
    if(outfilename==NULL) {
        fprintf(stderr,"ERROR: no output file name specified\n");
        exit(1);
    }
    /* auto flagging requires autocorrelations */
    if (autocorr_filename==NULL && do_flag) {
        fprintf(stderr,"ERROR: auto flagging requires the autocorrelations to be used\n");
        exit(1);
    }
}


/***************************
 ***************************/
void printusage(const char *progname) {
  fprintf(stderr,"Usage: %s [options]\n\n",progname);
  fprintf(stderr,"options are:\n");
  fprintf(stderr,"-a filename\tThe name of autocorrelation data file. no default.\n");
  fprintf(stderr,"-c filename\tThe name of cross-correlation data file. no default.\n");
  fprintf(stderr,"-o filename\tThe name of the output file. No default.\n");
  fprintf(stderr,"-S filename\tThe name of the file containing antenna name and local x,y,z. Default: %s\n",stationfilename);
  fprintf(stderr,"-I filename\tThe name of the file containing instrument config. Default: %s\n",configfilename);
  fprintf(stderr,"-H filename\tThe name of the file containing observing metadata. Default: %s\n",header_filename);
  fprintf(stderr,"-A lon,lat \tSpecify East Lon and Lat of array center (degrees). Comma separated, no spaces. Default: MWA\n");
  fprintf(stderr,"-l         \tLock the phase center to the initial HA/DEC\n");
  fprintf(stderr,"-f mode    \tturn on automatic flagging. Requires autocorrelations\n");
  fprintf(stderr,"\t\t0:\tno flagging\n");
  fprintf(stderr,"\t\t1:\tgeneric flagging\n");
  fprintf(stderr,"\t\t2:\tspecial treatment for MWA edge channels, 40kHz case\n");
  fprintf(stderr,"\t\t3:\tspecial treatment for MWA edge channels, 10kHz case\n");
  fprintf(stderr,"-F filename\tOptionally apply global flags as specified in filename.\n");
  fprintf(stderr,"-d         \tturn debugging on.\n");
  fprintf(stderr,"-v         \treturn revision number and exit.\n");
  exit(1);
}


/*******************************
*********************************/
void initData(uvdata *data) {
  data->date = calloc(1,sizeof(double));
  data->n_pol=0;
  data->n_baselines = calloc(1,sizeof(int));
  data->n_freq=0;
  data->n_vis=0;
  data->cent_freq=0.0;
  data->freq_delta = 0.0;
  data->u=NULL;
  data->v=NULL;
  data->w=NULL;
  data->baseline=NULL;
  data->weightdata=NULL;
  data->visdata=NULL;
  data->pol_type=1;    /* default is Stokes pol products */
  data->array->n_ant=0;
  strcpy(data->array->name,"UNKNOWN");
}


/***********************
***********************/
int applyHeader(Header *header, uvdata *data) {

  double jdtime_base=0,mjd,lmst;
  int n_polprod=0,i,res;

  data->n_freq = header->n_chans;
  data->cent_freq = header->cent_freq*1e6;
  data->freq_delta = header->bandwidth/(header->n_chans)*1e6*(header->invert_freq ? -1.0: 1.0);
  strncpy(data->array->name,header->telescope,SIZ_TELNAME);
  strncpy(data->array->instrument,header->instrument,SIZ_TELNAME);

  /* discover how many pol products there are */
  createPolIndex(header->pol_products, pol_index);
  for(i=0; i<4; i++) if (header->pol_products[2*i] !='\0') n_polprod++;
  if(n_polprod<1 || n_polprod > 4) {
    fprintf(stderr,"bad number of stokes: %d\n",n_polprod);
    exit(1);
  }
  data->n_pol = n_polprod;

  /* set the polarisation product type. linear, circular or stokes */
  if (toupper(header->pol_products[0]) == 'X' || toupper(header->pol_products[0]) == 'Y') data->pol_type=-5;
  if (toupper(header->pol_products[0]) == 'R' || toupper(header->pol_products[0]) == 'L') data->pol_type=-1;
  if (debug) fprintf(fpd,"Found %d pol products. pol_type is: %d\n",n_polprod,data->pol_type);

  /* calculate the JD of the beginning of the data */
  slaCaldj(header->year, header->month, header->day, &jdtime_base, &res); // get MJD for calendar day of obs
  jdtime_base += 2400000.5;  // convert MJD to JD
  jdtime_base += header->ref_hour/24.0+header->ref_minute/1440.0+header->ref_second/86400.0; // add intra-day offset
  data->date[0] = jdtime_base+0.5*(header->integration_time/86400.0);
  if (debug) fprintf(fpd,"JD time is %.2lf\n",jdtime_base);

  memset(data->source->name,0,SIZE_SOURCE_NAME+1);
  strncpy(data->source->name,header->field_name,SIZE_SOURCE_NAME);

  mjd = data->date[0] - 2400000.5;  // get Modified Julian date of scan.
  lmst = slaRanorm(slaGmst(mjd) + arr_lon_rad);  // local mean sidereal time, given array location

  /* if no RA was specified in the header,then calculate the RA based on lmst and array location
     and update the ra */
  if (header->ra_hrs < -98.0 ) {
    // set the RA to be for the middle of the scan
//    header->ra_hrs = lmst*(12.0/M_PI) - header->ha_hrs_start + header->n_scans*header->integration_time*1.00274/(3600.0*2);   // include 1/2 scan offset
    header->ra_hrs = lmst*(12.0/M_PI) - header->ha_hrs_start;  // match existing code. RA defined at start of scan
    if (debug) fprintf(fpd,"Calculated RA_hrs: %g of field centre based on HA_hrs: %g and lmst_hrs: %g\n",
                        header->ra_hrs,header->ha_hrs_start,lmst*(12.0/M_PI));
  }

  /* extract RA, DEC from header. Beware negative dec and negative zero bugs. */
  data->source->ra  = header->ra_hrs;
  data->source->dec = header->dec_degs;

  /* calcualte the number of baselines, required to be constant for all data */
  data->n_baselines[0] = (data->array->n_ant)*(data->array->n_ant+1)/2; //default: both auto and cross
  if (header->corr_type=='A') data->n_baselines[0] = data->array->n_ant;
  if (header->corr_type=='C') data->n_baselines[0] = data->array->n_ant*(data->array->n_ant-1)/2;
  if (debug) fprintf(fpd,"Corr type %c, so there are %d baselines\n",header->corr_type,data->n_baselines[0]);

  return 0;
}


/**************************
***************************/
void checkInputs(Header *header,uvdata *data,InpConfig *inputs) {
    int total_ants;

    if(inputs->n_inputs != header->n_inputs) {
        fprintf(stderr,"ERROR: mismatch between the number of inputs in %s (%d) and header (%d)\n",
                configfilename,inputs->n_inputs,header->n_inputs);
        exit(1);
    }

    total_ants = countPresentAntennas(inputs);

    if (total_ants > data->array->n_ant) {
        fprintf(stderr,"ERROR: mismatch between the number of antennas in %s (%d) and %s (%d)\n",
                stationfilename,data->array->n_ant,configfilename,total_ants);
    }
    
    if (do_flag && header->corr_type == 'C') {
        fprintf(stderr,"ERROR: CORRTYPE must be auto or both for autoflagging\n");
        exit(1);
    }
    if (header->ha_hrs_start == -99.0 && lock_pointing) {
        fprintf(stderr,"ERROR: HA must be specified in header if -l flag is used.\n");
        exit(1);
    }

}

