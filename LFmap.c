/*
/* LFmap v1.0
/*
/* Program to generate low frequency (20-400 MHz) sky maps
/* 
/* Reads input parameters from LFmap.config file
/*
/* Reads the EQ coord destriped 408 MHz Haslam map (04de_ecp.fit) and the spectral index map of Platania 
/* (bss_ecp.fit) and, if BENDFREQ > FINALFREQ, makes a new map scaled to BENDFREQ. New spectral indices
/* are calculated with this scaled map and the precessed & resampled 22 MHz map (drao22resampleded.fits) of Roger et al. 1999.
/* HII regions seen in absorption in the 22 MHz map are masked and scaled with an average index for the surrounding, unabsorbed
/* region. The map is scaled to ABSFREQ. At ABSFREQ the HII region is unmasked and spectral indices are calculated against the 
/* 22 MHz map. The map is then scaled to FINALFREQ.
/*
/* The output map is written in FITS and txt format. ("LFmap_FINALFREQ.*")
/*
/* If BENDFREQ < FINALFREQ Haslam's map is scaled directly to FINALFREQ with Platania's spectral indicies only, 
/* no corrections for spectral bending, no corrections for HII absorption.
/*
/* If ABSFREQ < FINALFREQ no HII absorption will be accounted for.
/*
/* created by Emil Polisensky, July 2007 
/* */

#include <stdio.h>
#include <stdlib.h>
#include "UtilLFmap.h"

main ()
{
  //call map generating function
  LFmapGen();
  printf("\nLFmap_1.0 finished!\n");
}
