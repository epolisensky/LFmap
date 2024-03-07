/*
/* UtilLFmap -- functions for LFmap.c v1.0
/*
/* Refer to header of LFmap.c file
/*
/* */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "fitsio.h"
#include "UtilLFmap.h"

#define PI 3.14159265359
#define NLINES 4 //number of lines in config file
#define STRLEN 256
#define NAX 2 //number of axes in 04de_ecp.fit & bss_ecp.fit

//user specificed parameters in LFmap.config file
struct params
{
  double FINALFREQ,BENDFREQ,ABSFREQ;
  int OUTFORM;
};

//map generating function
int LFmapGen(void)
{
  int j;
  long totpix,ii;
  double *pixNORD,*pix408,*pixIND,*ra408,*dec408,*area408,*raIND,*decIND,*areaIND,*tempCONV;
  double *pix22,*ra22,*dec22,*area22,*tempNEW;
  double minra,maxra,mindec,maxdec,inde,CASA,CYGA;
  FILE *fin,*outfile,*infile,*configfile;
  char fname[STRLEN],keywrd[STRLEN],line[STRLEN],fits_name[STRLEN];
  struct params par;

  //read config file
  sprintf(fname,"LFmap.config");
  printf("Opening config file [%s]\n",fname);
  configfile=fopen(fname,"r");
  for (j=0;j<NLINES;j++){
    if(j==0)
      fscanf(configfile,"%s %lf",&keywrd,&par.FINALFREQ);
    if(j==1)
      fscanf(configfile,"%s %lf",&keywrd,&par.BENDFREQ);
    if(j==2)
      fscanf(configfile,"%s %lf",&keywrd,&par.ABSFREQ);
    if(j==3)
      fscanf(configfile,"%s %d",&keywrd,&par.OUTFORM);
  } 
  fclose(configfile);
  printf("FINALFREQ = %lf\n",par.FINALFREQ);
  printf("BENDFREQ  = %lf\n",par.BENDFREQ);
  printf("ABSFREQ  = %lf\n",par.ABSFREQ);
  printf("OUTFORM  = %d\n",par.OUTFORM);

  totpix=1024*512;
  //allocate memory for arrays
  pix408   = (double *) malloc(totpix * sizeof(double)); /* memory for T408 values of all pixels */
  tempCONV = (double *) malloc(totpix * sizeof(double)); /* memory for convolved 408 map pixels */
  ra408    = (double *) malloc(totpix * sizeof(double)); /* memory for ra values of all pixels */
  dec408   = (double *) malloc(totpix * sizeof(double)); /* memory for dec values of all pixels */
  area408  = (double *) malloc(totpix * sizeof(double)); /* memory for sky area values of all pixels */
  pixIND   = (double *) malloc(totpix * sizeof(double)); /* memory for index values of all pixels */
  raIND    = (double *) malloc(totpix * sizeof(double)); /* memory for ra values of all pixels */
  decIND   = (double *) malloc(totpix * sizeof(double)); /* memory for dec values of all pixels */
  areaIND  = (double *) malloc(totpix * sizeof(double)); /* memory for sky area values of all pixels */
  pix22    = (double *) malloc(totpix * sizeof(double)); /* memory for T22 values of all pixels */
  ra22     = (double *) malloc(totpix * sizeof(double)); /* memory for ra values of all pixels */
  dec22    = (double *) malloc(totpix * sizeof(double)); /* memory for dec values of all pixels */
  area22   = (double *) malloc(totpix * sizeof(double)); /* memory for sky area values of all pixels */
  tempNEW  = (double *) malloc(totpix * sizeof(double)); /* memory for scaled temp values of all pixels */
  pixNORD  = (double *) malloc(totpix * sizeof(double)); /* memory for T_Nord values of all pixels */
  
  printf("starting code\n");
  //read destriped Haslam map in EQ coords and J2000.0
  sprintf(fits_name,"04de_ecp.fit");
  readPlataniaMaps(pix408,ra408,dec408,area408,totpix,fits_name);
  //read spectral index map (indices of Galactic emission)
  sprintf(fits_name,"bss_ecp.fit");
  readPlataniaMaps(pixIND,raIND,decIND,areaIND,totpix,fits_name);
  //read resampled 22 MHz map
  sprintf(fits_name,"drao22resampled.fits");
  readPlataniaMaps(pix22,ra22,dec22,area22,totpix,fits_name);
  printf("read fits maps\n");
  //subtract CMB and extragal emission from Halsam and 22 MHz maps
  subtractEG(pix408,408.0,totpix);
  subtractEG(pix22,22.0,totpix);
  //read Mike's T74 map (for Cas A and Cyg A corrections later)
  sprintf(fname,"T74_B.dat");
  infile=fopen(fname,"r");
  if(infile==NULL){
    fprintf(stderr,"Error opening file %s\n",fname);
    exit(94);  }
  if (pixNORD == NULL) {
    printf("Memory allocation error\n");
    return(1);  }
  ii=0;
  while (fgets(line,100,infile) != NULL) {
    sscanf(line,"%lf",&pixNORD[ii]);
    ii++;
  }
  fclose(infile);
  ////////////////
  //now begin generating final map
  if (par.BENDFREQ > par.FINALFREQ){//scale to BENDFREQ first
    if (par.BENDFREQ < par. ABSFREQ){
      fprintf(stderr,"error! ABSFREQ=%lf  BENDFREQ=%lf\nbut ABSFREQ should be less than BENDFREQ!\n",par.ABSFREQ,par.BENDFREQ);
      exit(1); }
    for(ii=0;ii<totpix;ii++){
      if (decIND[ii] < -79.6) //for incomplete portion of bss_ecp.fit map
	inde=2.695;
      else
	inde=pixIND[ii];
      tempNEW[ii]=pix408[ii]*pow((408.0/par.BENDFREQ),inde);
    }
    //convolve Platania map to 22 map resolution (only convolve overlap region)
    convolve(tempNEW,tempCONV,ra408,dec408,totpix);
    //now calc new indices with 22 map 
    calcnewIND(tempCONV,pix22,pixIND,par.BENDFREQ,22.0,totpix);
    //deal with incomplete portions of map
    fillincomplete(pixIND,pix22,raIND,decIND,totpix);
    //deal with HII absorption, if necessary
    if (par.FINALFREQ < par.ABSFREQ){//scale to ABSFREQ
      for(ii=0;ii<totpix;ii++){
	tempNEW[ii]=tempNEW[ii]*pow((par.BENDFREQ/par.ABSFREQ),pixIND[ii]);
      }
      //calc new indices for absorption regions
      calcabsIND(tempNEW,pix22,pixIND,par.ABSFREQ,22.0,totpix,raIND,decIND);
      //now scale to FINALFREQ
      for(ii=0;ii<totpix;ii++){
	tempNEW[ii]=tempNEW[ii]*pow((par.ABSFREQ/par.FINALFREQ),pixIND[ii]);
      }
    }
    else{//no need to scale to ABSFREQ, just scale to FINALFREQ
      for(ii=0;ii<totpix;ii++){
	tempNEW[ii]=tempNEW[ii]*pow((par.BENDFREQ/par.FINALFREQ),pixIND[ii]);
      }
    }   
  }
  else{//scale to FINALFREQ with just Platania indicies
    for(ii=0;ii<totpix;ii++){
      if (decIND[ii] < -79.6) //for incomplete portion of map
	inde=2.695;
      else
	inde=pixIND[ii];
      tempNEW[ii]=pix408[ii]*pow((408.0/par.FINALFREQ),inde);
    }
  }
  
  //add extragal emission and CMB to final map
  addEG(tempNEW,par.FINALFREQ,totpix);

  //scale CAS A & CYG A using Nord's T74 values and these indicies:
  CASA=2.77;CYGA=2.66; //from Whitfield, 1957
  for(ii=0;ii<totpix;ii++){
    if (ii<=379566 && ii>=379559)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CYGA);
    if (ii<=380590 && ii>=380584)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CYGA);
    if (ii<=381614 && ii>=381607)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CYGA);
    if (ii<=382637 && ii>=382633)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CYGA);
    if (ii<=430619 && ii>=430616)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CASA);
    if (ii<=431643 && ii>=431638)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CASA);
    if (ii<=432668 && ii>=432659)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CASA);
    if (ii<=433697 && ii>=433684)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CASA);
    if (ii<=434718 && ii>=434708)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CASA);
    if (ii<=435737 && ii>=435735)
      tempNEW[ii]=pixNORD[ii]*pow((74.0/par.FINALFREQ),CASA);
  }

  //write outfile
  if (par.OUTFORM==3 || par.OUTFORM==4){
      printf("OUTFORM = %d, writing FITS file\n",par.OUTFORM);
      writefitsout(tempNEW,totpix,par.FINALFREQ,par.BENDFREQ,par.ABSFREQ,par.OUTFORM);
  }
  else if (par.OUTFORM==1 || par.OUTFORM==2){
      printf("OUTFORM = %d, writing txt file\n",par.OUTFORM);
      writetxtout(ra408,dec408,tempNEW,totpix,par.FINALFREQ,par.OUTFORM);
  }
  else
      printf("error! OUTFORM = %d not recognized!\n",par.OUTFORM);
  
  free(pix22);   free(area22);  free(ra22);   free(dec22);
  free(pix408);  free(area408); free(ra408);  free(dec408);
  free(pixIND);  free(raIND);   free(decIND); free(areaIND);
  free(pixNORD); free(tempNEW); free(tempCONV);

  return 0;
}//end of LFmapGen


/* routine to calculate the surface area of a pixel, given the RA and DEC of its corners */
/* assumes RA & DEC values are in degrees */
double surfareaDEG(double ra1, double ra2, double dec1, double dec2)
{
  double area,R;
  R = 1.0; // unit sphere --> surfarea=solidangle
  area = R*R*fabs( sin(dec1*PI/180.0) - sin(dec2*PI/180.0) )*fabs((ra1-ra2)*PI/180.0);
  return area;
}

// wrtie output map in FITS format
int writefitsout(double *tempNEW, long totpix, double ffreq, double tfreq, double afreq, int OUTFORM)
{
  fitsfile *fptr;
  long naxes[2]={ 1024, 512 }; //image is 1024 pixels wide by 512 rows
  long fpixel[2],ii,jj,lpixel[2],naxis=2;
  int status=0,keytype,bitpix=FLOAT_IMG; //32 bit single precision floating point
  char filename[FLEN_CARD];
  char card[FLEN_CARD],comment[FLEN_COMMENT],keyword[FLEN_KEYWORD];
  float *buffer,value;
  
  buffer = (float *) malloc(naxes[0]*sizeof(float));

  if (totpix != naxes[0]*naxes[1])
    printf("error! totpix doesn't equal size of FITS image!\n");

  sprintf(filename,"LFmap_%.1lf.fits",ffreq);
  printf("%s\n",filename);

  //create and open new FITS file
  remove(filename); //delete old file if already exists
  if (fits_create_file(&fptr,filename,&status))
    printf("open error\n");

  //create image
  if ( fits_create_img(fptr, bitpix, naxis, naxes, &status) )
    printerror( status );

  //write image
  fpixel[0]=1; 
  lpixel[0]=naxes[0];

  for(jj=0;jj<naxes[1];jj++){    
    for(ii=0;ii<naxes[0];ii++){
      buffer[ii]=(float)(tempNEW[(jj*naxes[0])+ii]);
    }
    fpixel[1]=jj+1;
    lpixel[1]=fpixel[1];
    if ( fits_write_subset(fptr, TFLOAT, fpixel, lpixel, buffer, &status) )
      printerror( status );
  }
  free(buffer);

  //write header
  if (OUTFORM==3)
      strcpy(card,"RA---CAR");
  if (OUTFORM==4)
      strcpy(card,"RA");
  ffpkys(fptr,"CTYPE1",card,"X-axis type",&status);
  if (OUTFORM==3)
      strcpy(card,"DEC--CAR");
  if (OUTFORM==4)
      strcpy(card,"DEC");
  ffpkys(fptr,"CTYPE2",card,"Y-axis type",&status);
  value=0.0;
  strcpy(keyword,"CRVAL1");
  strcpy(comment,"X-axis Reference pixel value");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  value=0.0;
  strcpy(keyword,"CRVAL2");
  strcpy(comment,"Y-axis Reference pixel value");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  value=512.5;
  strcpy(keyword,"CRPIX1");
  strcpy(comment,"X-axis Reference pixel");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  value=256.5;
  strcpy(keyword,"CRPIX2");
  strcpy(comment,"Y-axis Reference pixel");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  value=-0.351562;
  strcpy(keyword,"CDELT1");
  strcpy(comment,"X-axis degrees/pixel");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  value=0.351562;
  strcpy(keyword,"CDELT2");
  strcpy(comment,"Y-axis degrees/pixel");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  value=1.0;
  strcpy(keyword,"CD001001");
  strcpy(comment,"no rotation or skew");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  value=1.0;
  strcpy(keyword,"CD002002");
  strcpy(comment,"no rotation or skew");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  strcpy(card,"FK5");
  ffpkys(fptr,"RADESYS",card,"",&status);
  value=2000.0;
  strcpy(keyword,"EQUINOX");
  strcpy(comment,"Equinox of coordinate system (years)");
  fits_write_key(fptr,TFLOAT,keyword,&value,comment,&status);
  strcpy(card,"K");
  ffpkys(fptr,"BUNIT",card,"Brightness temperature units are Kelvins",&status);
  sprintf(card,"%lf MHz",ffreq);
  ffpkys(fptr,"FINALFRQ",card,"Final freq map is scaled to.",&status);
  sprintf(card,"%lf MHz",tfreq);
  ffpkys(fptr,"BENDFREQ",card,"Freq where new spectral indices are calculated",&status);
  sprintf(card,"%lf MHz",afreq);
  ffpkys(fptr,"ABSFREQ",card,"Freq below which HII abs is taken into account",&status);
  //write comments
  sprintf(card,"COMMENT Frequency - %.1lfMHz",ffreq);
  fits_write_record(fptr,card,&status);
  fits_write_record(fptr,"COMMENT Pixel Scale - 0.3515625 degrees/pixel",&status);
  fits_write_record(fptr,"COMMENT Pixel Units - Kelvins",&status);
  fits_write_record(fptr,"COMMENT Coordinate System - Equatorial",&status);
  fits_write_record(fptr,"COMMENT Projection - Plate carree (Rectangular)",&status);
  fits_write_record(fptr,"COMMENT ---------------------------------",&status);
  fits_write_record(fptr,"COMMENT   This map created by scaling the destriped Haslam 408 MHz",&status);
  fits_write_record(fptr,"COMMENT   map with the spectral indices from Platania et al,2003,A&A,410,847,",&status);
  fits_write_record(fptr,"COMMENT   to the BENDFREQ. At the BENDFREQ the resampled 22 MHz Map is used ",&status);
  fits_write_record(fptr,"COMMENT   to calculate new spectral indices and scale to FINALFREQ.",&status);
  fits_write_record(fptr,"COMMENT   If BENDFREQ < FINALFREQ the map is scaled using ",&status);
  fits_write_record(fptr,"COMMENT   just the Platania indicies.",&status);
  fits_write_record(fptr,"COMMENT   If ABSFREQ < FINALFREQ  HII absorption is not taken",&status);
  fits_write_record(fptr,"COMMENT   into account.",&status);
  fits_write_record(fptr,"COMMENT ---------------------------------",&status);
  strcpy(card,"This map created by LFmap.c v1.0");
  fits_write_history(fptr,card,&status);
  fits_write_date(fptr,&status);

  //close map file
  if (fits_close_file(fptr,&status)){
    fits_report_error(stderr, status); /* print any error message */
    printf("close error\n");
  }

  return 0;
}

int writetxtout(double *ra, double *dec, double *tempNEW, long totpix, double ffreq, int OUTFORM)
{
  char filename[STRLEN];
  FILE *fout;
  long ii;
  double k1,*gl,*gb,tmp1,tmp2,tmp3;

  if (OUTFORM==1){    //write txt file w/ EQ coords
      sprintf(filename,"LFmap_%.1lf_EQ.txt",ffreq);
      fout=fopen(filename,"w");
      if(fout==NULL){
	  fprintf(stderr,"Error opening file %s\n",filename);
	  exit(94);  }
      for(ii=0;ii<totpix;ii++){
	  if (ra[ii] > 180.0)
	      ra[ii]=ra[ii]-360.0;
	  fprintf(fout,"%d %lf %lf %lf\n",ii,ra[ii],dec[ii],tempNEW[ii]);
      }
      fclose(fout);
  }
  else if (OUTFORM==2){    //write txt file w/ GAL coords
      gl =(double *) malloc(totpix*sizeof(double));
      gb =(double *) malloc(totpix*sizeof(double));
      k1=62.6*PI/180.0;
      sprintf(filename,"LFmap_%.1lf_GAL.txt",ffreq);
      fout=fopen(filename,"w");
      if(fout==NULL){
	  fprintf(stderr,"Error opening file %s\n",filename);
	  exit(94);  }
      for(ii=0;ii<totpix;ii++){
	  //convert to galactic coords (l,b)
	  tmp1=(sin(dec[ii]*PI/180.0)*cos(k1)) - (cos(dec[ii]*PI/180.0)*sin((ra[ii]-282.25)*PI/180.0)*sin(k1));
	  gb[ii]=(180.0/PI)*(asin(tmp1));
	  tmp2=(sin(dec[ii]*PI/180.0)*sin(k1)) + (cos(dec[ii]*PI/180.0)*sin((ra[ii]-282.25)*PI/180.0)*cos(k1));
	  tmp3=cos(dec[ii]*PI/180.0)*cos((ra[ii]-282.25)*PI/180.0);
	  gl[ii]=33.0+((180.0/PI)*atan2(tmp2,tmp3));
	  if (gl[ii] > 180.0)
	      gl[ii]=gl[ii]-360.0;
	  fprintf(fout,"%d %lf %lf %lf\n",ii,gl[ii],gb[ii],tempNEW[ii]);
      }
      fclose(fout);
      free(gl);free(gb);
  }
  else{
      printf("error in writetxtout! OUTFORM = %d not recognized!\n",OUTFORM);
      return 1;
  }

  printf("%s\n",filename);
  return 0;
}

void printerror(int status)
{
    /* Print out cfitsio error messages and exit program */
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

int readPlataniaMaps(double *pix, double *ra, double *dec, double *area, int totpix, char *fits_name)
{
  fitsfile *fptr;
  int status=0;
  int hdutype,naxis;
  double nulval = -1.0;
  double anynul=0.0;
  char errstr[STRLEN], commstr[STRLEN];
  char *statlist;
  long naxes[NAX], fpixel[NAX], lpixel[NAX], inc[NAX], i;

  long ii,column;
  double refptra,refptdec;
  double *ip,dum1,deltara,deltadec,ra2,ra1,dec2,dec1,tmp3,tarea2=0.0;
  double minra,maxra,mindec,maxdec;
  int row,k;
  FILE *fin,*infile;
  char fname[STRLEN],line[STRLEN];

  ip = malloc(sizeof(double));

  //  printf("reading Platania map...\n");
  //read all pixels, store in a 1D array, calculate RA and DEC for each pixel center
  fits_open_file(&fptr,fits_name,READONLY,&status); 
  if (fits_get_hdu_type(fptr, &hdutype, &status) || hdutype != IMAGE_HDU) { 
    printf("Error: this program only works on images, not tables\n");
    return(1);
  }
  fits_get_img_dim(fptr, &naxis, &status);
  fits_get_img_size(fptr, naxis, naxes, &status);

  if (status || naxis != NAX) { 
    printf("NAXIS = %d.  Only %d-D images supported\n", naxis,NAX);
    return(1);
  }

  for (k=0;k<NAX;k++){
    fpixel[k]=1;
    inc[k]=1;
    lpixel[k]=naxes[k];
  }

  //  printf("reading pixels...\n");
  fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &nulval, pix, 0, &status);
  //printf("pixels read\n");
  fits_close_file(fptr,&status);

  refptra=512.5; //ra=0 here
  refptdec=256.5;//dec=0 here
  deltara = -0.3515625; //step size (both axes) in decimal degrees
  deltadec = 0.3515625;

  if (status)  {
    fits_report_error(stderr, status); /* print any error message */
  }

  minra=1000;maxra=0;maxdec=-90.0;mindec=90.0;
  for (ii=0; ii<totpix; ii++){
    //convert ii to pixel row and column
    dum1=(double)ii/naxes[0];
    tmp3=((double)naxes[0]*(modf(dum1,ip)))+1.0;
    column=(long)tmp3;
    row=(int)(1+(*ip));
    //convert pixel ROW and COL to ra and dec
    ra[ii]=(((double)tmp3)-refptra)*deltara;
    dec[ii]=(((double)row)-refptdec)*deltadec;

    if (ra[ii] < 0.0)
      ra[ii]+=360.0;

    if (ra[ii] > maxra)
      maxra=ra[ii];
    if(ra[ii] < minra)
      minra=ra[ii];
    if (dec[ii] > maxdec)
      maxdec=dec[ii];
    if(dec[ii] < mindec)
      mindec=dec[ii];

    //calculate area
    ra2=ra[ii]+0.5*deltara;
    ra1=ra[ii]-0.5*deltara;
    dec1=dec[ii]+0.5*deltadec;
    dec2=dec[ii]-0.5*deltadec;
    
    area[ii]=surfareaDEG(ra1,ra2,dec1,dec2);
    tarea2+=area[ii];

  } //end of ii loop

  return 0;
}

//method of calculating extragalactic emssion from Lawson et al. 1987
int subtractEG(double *tempNEW, double freq, long totpix)
{
  double Tcmb=2.73,T150=50.0,indEG=2.75,TEG;
  long ii;

  TEG=T150*pow((150.0/freq),indEG);
  for(ii=0;ii<totpix;ii++){
    tempNEW[ii]=tempNEW[ii]-Tcmb-TEG;
  }
  return 0;
}

int addEG(double *tempNEW, double freq, long totpix)
{
  double Tcmb=2.73,T150=50.0,indEG=2.75,TEG;
  long ii;

  TEG=T150*pow((150.0/freq),indEG);
  for(ii=0;ii<totpix;ii++){
    tempNEW[ii]=tempNEW[ii]+Tcmb+TEG;
  }
  return 0;
}

int calcnewIND(double *tempNEW, double *pix22, double *pixIND, double turnfreq, double reffreq, long totpix)
{
  long ii;

  for(ii=0;ii<totpix;ii++){
    if (pix22[ii] > 0.0)
      pixIND[ii]=(log10(pix22[ii]/tempNEW[ii]))/(log10(turnfreq/reffreq));
  }
  return 0;
}

//calc indices for HII absorbed regions only
int calcabsIND(double *tempNEW, double *pix22, double *pixIND, double absfreq, double reffreq, long totpix, double *ra, double *dec)
{
  long ii;
  double *gl,*gb,tmp1,tmp2,tmp3,k1;

  gl =(double *) malloc(totpix*sizeof(double));
  gb =(double *) malloc(totpix*sizeof(double));

  k1=62.6*PI/180.0;
  for(ii=0;ii<totpix;ii++){
    //convert to galactic coords (l,b)
    tmp1=(sin(dec[ii]*PI/180.0)*cos(k1)) - (cos(dec[ii]*PI/180.0)*sin((ra[ii]-282.25)*PI/180.0)*sin(k1));
    gb[ii]=(180.0/PI)*(asin(tmp1));
    tmp2=(sin(dec[ii]*PI/180.0)*sin(k1)) + (cos(dec[ii]*PI/180.0)*sin((ra[ii]-282.25)*PI/180.0)*cos(k1));
    tmp3=cos(dec[ii]*PI/180.0)*cos((ra[ii]-282.25)*PI/180.0);
    gl[ii]=33.0+((180.0/PI)*atan2(tmp2,tmp3));
    if (gl[ii] > 180.0)
      gl[ii]=gl[ii]-360.0;
  }

  for(ii=0;ii<totpix;ii++){
    if (gb[ii] < 6.0 && gb[ii] > -4.0 && gl[ii] > 5.0 && gl[ii] < 55.0 && pix22[ii] > 0.0)//region contaminated by HII abs.
      pixIND[ii]=(log10(pix22[ii]/tempNEW[ii]))/(log10(absfreq/reffreq));
  }

  free(gl);  free(gb);
  return 0;
}

int fillincomplete(double *pixIND, double *pix22, double *ra, double *dec, long totpix)
{
  long ii,c1,c2,c3,cTAU,cCYG,cVIR,c4;
  double tmp1,tmp2,tmp3,*gb,*gl,k1,sum1,sum2,sum3,ind1,ind2,ind3;
  double indCYG,indTAU,indVIR;
  double sumCYG,sumTAU,sumVIR;

  gl =(double *) malloc(totpix*sizeof(double));
  gb =(double *) malloc(totpix*sizeof(double));

  k1=62.6*PI/180.0;

  for(ii=0;ii<totpix;ii++){
    //convert to galactic coords (l,b)
    tmp1=(sin(dec[ii]*PI/180.0)*cos(k1)) - (cos(dec[ii]*PI/180.0)*sin((ra[ii]-282.25)*PI/180.0)*sin(k1));
    gb[ii]=(180.0/PI)*(asin(tmp1));
    tmp2=(sin(dec[ii]*PI/180.0)*sin(k1)) + (cos(dec[ii]*PI/180.0)*sin((ra[ii]-282.25)*PI/180.0)*cos(k1));
    tmp3=cos(dec[ii]*PI/180.0)*cos((ra[ii]-282.25)*PI/180.0);
    gl[ii]=33.0+((180.0/PI)*atan2(tmp2,tmp3));
    if (gl[ii] > 180.0)
      gl[ii]=gl[ii]-360.0;
  }

  c1=c2=c3=c4=0;
  sum1=sum2=sum3=0.0;
  cCYG=cTAU=cVIR=0;
  sumCYG=sumTAU=sumVIR=0.0;
  //break up incomplete regions into subregions and set indicies in subregions according to class
  //
  //first break up known regions of map into classes and calc average index in these regions
  for(ii=0;ii<totpix;ii++){
    //GC & Bulge region
    if (gb[ii] < 6.0 && gb[ii] > -4.0 && gl[ii] > 5.0 && gl[ii] < 55.0){//ignore this region, contaminated by HII abs.
      c4++;
    }
    else if (fabs(gb[ii]) < 10.0 && gl[ii] > -40.0 && gl[ii] < 60.0 && pix22[ii] > 0.0 ){
      sum1+=pixIND[ii]; c1++;
    }
    else if (gb[ii] > -20.0 && gb[ii] < 30.0 && gl[ii] > -40.0 && gl[ii] < -20.0 && pix22[ii] > 0.0 ){
      sum1+=pixIND[ii]; c1++;
    }
    else if (fabs(gb[ii]) < 20.0 && gl[ii] > -20.0 && gl[ii] < 50.0 && pix22[ii] > 0.0 ){
      sum1+=pixIND[ii]; c1++;
    }
    //min/cool regions
    else if (gb[ii] > -70.0 && gb[ii] < -18.0 && gl[ii] > -130.0 && gl[ii] < -120.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    else if (gb[ii] > -70.0 && gb[ii] < -25.0 && gl[ii] > -150.0 && gl[ii] < -130.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    else if (gb[ii] > -70.0 && gb[ii] < -35.0 && gl[ii] > -170.0 && gl[ii] < -150.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    else if (gb[ii] > 40.0 && gb[ii] < 70.0 && gl[ii] > -120.0 && gl[ii] < -100.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    else if (gb[ii] > 15.0 && gb[ii] < 70.0 && gl[ii] > -150.0 && gl[ii] < -120.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    else if (gb[ii] > 30.0 && gb[ii] < 70.0 && gl[ii] < -150.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    else if (gb[ii] > 40.0 && gb[ii] < 70.0 && gl[ii] > 170.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    else if (gb[ii] > 45.0 && gb[ii] < 65.0 && gl[ii] > 160.0 && gl[ii] < 170.0 && pix22[ii] > 0.0){
      sum2+=pixIND[ii]; c2++;
    }
    //loops, spurs, and other Gal emission
    else if (dec[ii] > 17.0 && dec[ii] < 26.0 && ra[ii] < 88.0 && ra[ii] > 77.0 && pix22[ii] > 0.0){//Tau A
      sumTAU+=pixIND[ii]; cTAU++;
    }
    else if (gb[ii] > -20.0 && gb[ii] < 20.0 && gl[ii] > 120.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    else if (gb[ii] > -40.0 && gb[ii] < -10.0 && gl[ii] > 50.0 && gl[ii] < 90.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    else if (gb[ii] > -40.0 && gb[ii] < -20.0 && gl[ii] > 0.0 && gl[ii] < 50.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    else if (gb[ii] > -20.0 && gb[ii] < 20.0 && gl[ii] < -160.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    else if (gb[ii] > -10.0 && gb[ii] < 10.0 && gl[ii] > -160.0 && gl[ii] < -110.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    else if (gb[ii] > 30.0 && gb[ii] < 50.0 && gl[ii] > -60.0 && gl[ii] < 10.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    else if (gb[ii] > 20.0 && gb[ii] < 30.0 && gl[ii] > -20.0 && gl[ii] < 10.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    else if (gb[ii] > -55.0 && gb[ii] < -25.0 && gl[ii] > 145.0 && gl[ii] < 165.0 && pix22[ii] > 0.0){
      sum3+=pixIND[ii]; c3++;
    }
    //discrete sources
    else if (dec[ii] > 34.0 && dec[ii] < 47.0 && ra[ii] < 307.0 && ra[ii] > 292.0 && pix22[ii] > 0.0){//Cyg A
      sumCYG+=pixIND[ii]; cCYG++;
    }
    else if (dec[ii] > 8.0 && dec[ii] < 17.0 && ra[ii] < 192.5 && ra[ii] > 182.5 && pix22[ii] > 0.0){//Vir A
      sumVIR+=pixIND[ii]; cVIR++;
    }

  }//end for
  ind1=sum1/(double)c1;  ind2=sum2/(double)c2;  ind3=sum3/(double)c3;
  //printf("%d %d %d \n%.2lf %.2lf %.2lf   c4=%d\n",c1,c2,c3,ind1,ind2,ind3,c4);

  if (cCYG==0 || cTAU==0 || cVIR==0){
    printf("error! cCYG or cTAU or cVIR = 0!!\n");
    printf("%d %d %d \n",cCYG,cTAU,cVIR);
  }
  else{
    indCYG=sumCYG/(double)cCYG;    indTAU=sumTAU/(double)cTAU;    indVIR=sumVIR/(double)cVIR;
    //printf("cCYG=%d indCYG=%lf   cVIR=%d indVIR=%lf   cTAU=%d indTAU=%lf\n",cCYG,indCYG,cVIR,indVIR,cTAU,indTAU);
  }

  //set indicies in incomplete regions
  for(ii=0;ii<totpix;ii++){
      //GC and Bulge
    if (gb[ii] < 6.0 && gb[ii] > -4.0 && gl[ii] > 5.0 && gl[ii] < 55.0)//region contaminated by HII abs, set to ind1
      pixIND[ii]=ind1;
    else if (fabs(gb[ii]) < 10.0 && gl[ii] > -100.0 && gl[ii] < -60.0 && pix22[ii] < 0.0)
      pixIND[ii]=ind1;
    else if (gb[ii] > -20.0 && gb[ii] < 30.0 && gl[ii] > -60.0 && gl[ii] < 5.0 && pix22[ii] < 0.0)
      pixIND[ii]=ind1;
    else if (gb[ii] > -20.0 && gb[ii] < -4.0 && gl[ii] > 5.0 && gl[ii] < 20.0 && pix22[ii] < 0.0)
      pixIND[ii]=ind1;
    //min, SCP, and NCP regions
    else if (gb[ii] < -60.0 && gl[ii] < -40.0 && gl[ii] > -90.0 && pix22[ii] < 0.0)
      pixIND[ii]=ind2;
    else if (gb[ii] < -50.0 && gl[ii] < 35.0 && gl[ii] > -40.0 && pix22[ii] < 0.0)
      pixIND[ii]=ind2;
    else if (gb[ii] < -20.0 && gl[ii] > -150.0 && gl[ii] < -90.0 && pix22[ii] < 0.0)
      pixIND[ii]=ind2;
    else if (dec[ii] > 75.0 && pix22[ii] < 0.0) //NCP
      pixIND[ii]=ind2;
    //loops, spurs, and other Gal regions
    else if (dec[ii] < -27.0 && pix22[ii] < 0.0)
      pixIND[ii]=ind3;
    //discrete sources
    else if (dec[ii] > 5.0 && dec[ii] < 20.0 && ra[ii] < 195.0 && ra[ii] > 180.0 && pix22[ii] < 0.0)//Virgo A
      pixIND[ii]=indVIR;
    else if (dec[ii] > 30.0 && dec[ii] < 50.0 && ra[ii] < 310.0 && ra[ii] > 290.0 && pix22[ii] < 0.0)//Cyg A
      pixIND[ii]=indCYG;
    else if (dec[ii] > 15.0 && dec[ii] < 30.0 && ra[ii] < 95.0 && ra[ii] > 75.0 && pix22[ii] < 0.0)//Tau A
      pixIND[ii]=indTAU;
  }
  free(gl);  free(gb);
  return 0;
}

int convolve(double *temp, double *tempCONV, double *ra, double *dec, long totpix)
{
  int i,j,k,p,q,case1,case2,case3,wx=4,wy=9,npts=(1+2*wx)*(1+2*wy);
  int pixlist[npts];
  double WIDTH,HEIGHT,ZEN,ZA,dRA,dDEC,DECMIN,DECMAX,area,sum,totalarea;
  double declow,dechigh,ralow,rahigh;

  dRA=dDEC=0.3515625;
  DECMIN=-28.4;
  DECMAX=80.4;
  ZEN=48.8;
  WIDTH=1.1;

  printf("convolving map...\n");
  //convolve map to 22 MHz map resolution
  for(i=0;i<totpix;i++){
    if ((dec[i] < DECMIN) || (dec[i] > DECMAX))
      tempCONV[i]=temp[i];//default value for pixels outside precessed 22 MHz map
    else{
      //determine which j pixels are contained partly or wholly in i pixel.
      totalarea=sum=0.0;
      ZA=abs(ZEN-dec[i]);
      HEIGHT=1.7/cos(ZA*PI/180.0);

      k=0;
      //generate list of neighboring pixels
      for(q=(-1*wy);q<=wy;q++){
	for(p=(-1*wx);p<=wx;p++){
	  pixlist[k]=(i+p)+(q*1024);
	  k++;
	}
      }
      for(j=0;j<npts;j++){
	case1=case2=case3=0;

	if (fabs(ra[i]-ra[pixlist[j]]) < 0.5*(dRA+WIDTH))
	  case1=1;
	if ( (ra[i]-ra[pixlist[j]]) < (-1.0*(360.0-0.5*(dRA+WIDTH))) )
	  case2=1; //small ra[i], wrap around pixlist[j]
	if ( (ra[i]-ra[pixlist[j]]) > 360.0-0.5*(dRA+WIDTH) )
	  case3=1; //large ra[i], wrap around pixlist[j]
	if (case1 && (case2 || case3))
	  printf("error! cannot be both case 1 and case 2 or 3! %d %d %d\n",case1,case2,case3);

	if (case1 || case2 || case3){
	    if ( fabs(dec[i]-dec[pixlist[j]]) < 0.5*(dDEC+HEIGHT)){
	      //then at least some of pixlist[j] is in i
	      //is pixlist[j] wholly in i or only partly? find area of pix pixlist[j] contained in pix i (find coordinates of corners)
	      if (case1){
		if ((ra[pixlist[j]]-(0.5*dRA)) >= (ra[i]-(0.5*WIDTH)))
		  ralow=ra[pixlist[j]]-(0.5*dRA);
		else
		  ralow=ra[i]-(0.5*WIDTH);
		if ((ra[pixlist[j]]+(0.5*dRA)) <= (ra[i]+(0.5*WIDTH)))
		  rahigh=ra[pixlist[j]]+(0.5*dRA);
		else
		  rahigh=ra[i]+(0.5*WIDTH);
	      }
	      else if (case2){//case2 = small ra[i]  
		rahigh=ra[pixlist[j]]+(0.5*dRA);
		if ( (ra[pixlist[j]]-(0.5*dRA)) >= (360.0+(ra[i]-(0.5*WIDTH))) )
		  ralow=ra[pixlist[j]]-(0.5*dRA);
		else
		  ralow=360.0+(ra[i]-(0.5*WIDTH));
	      }
	      else if (case3){//case3 = large ra[i]
		ralow=ra[pixlist[j]]-(0.5*dRA);
		if ((ra[pixlist[j]]+(0.5*dRA)) <= (360.0-(ra[i]+(0.5*WIDTH))) )
		  rahigh=ra[pixlist[j]]+(0.5*dRA);
		else
		  rahigh=360.0-(ra[i]+(0.5*WIDTH));
	      }
	      else
		printf("error! neither case 1 or case2 or case3!\n");

	      if ((dec[pixlist[j]]-(0.5*dDEC)) >= (dec[i]-(0.5*HEIGHT)))
		declow=dec[pixlist[j]]-(0.5*dDEC);
	      else
		declow=dec[i]-(0.5*HEIGHT);
	      if ((dec[pixlist[j]]+(0.5*dDEC)) <= (dec[i]+(0.5*HEIGHT)))
		dechigh=dec[pixlist[j]]+(0.5*dDEC);
	      else
		dechigh=dec[i]+(0.5*HEIGHT);
	      
	      area=surfareaDEG(ralow,rahigh,declow,dechigh);
	      sum+=temp[pixlist[j]]*area;
	      totalarea+=area;
	    }
	  }
      }//end j for loop

      tempCONV[i]=sum/totalarea;

    }//end else

  }//end ii loop

  printf("convolved map\n");

  return 0;
}
