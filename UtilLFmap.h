/* UtilLFmap.h */
//
int readPlataniaMaps(double *pix, double *ra, double *dec, double *area, int totpix, char *fits_name);
double surfareaDEG(double ra1, double ra2, double dec1, double dec2);
int writefitsout(double *tempNEW, long totpix, double ffreq, double tfreq, double afreq, int OUTFORM);
int writetxtout(double *ra, double *dec, double *tempNEW, long totpix, double ffreq, int OUTFORM);
void printerror(int status);
int subtractEG(double *tempNEW, double freq, long totpix);
int addEG(double *tempNEW, double freq, long totpix);
int calcnewIND(double *tempNEW, double *pix22, double *pixIND, double turnfreq, double reffreq, long totpix);
int calcabsIND(double *tempNEW, double *pix22, double *pixIND, double absfreq, double reffreq, long totpix, double *ra, double *dec);
int fillincomplete(double *pixIND, double *pix22, double *ra, double *dec, long totpix);
int LFmapGen(void);
int convolve(double *temp, double *tempCONV, double *ra, double *dec, long totpix);
