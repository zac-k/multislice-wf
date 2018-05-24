/*      *** image.cpp ***

------------------------------------------------------------------------
Copyright 1998-2012 Earl J. Kirkland

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 

------------------------------------------------------------------------

  ANSI C and TIFF version
  this version uses FFTW 3

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -o image image.cpp slicelib.o
                       floatTIFF.cpp -lfftw3f


  rewritten in ANSI-C 29-sep-1995 E. Kirkland
  convert to TIFF file format 20-may-1996 ejk
  remove commas in input format 11-july-1997 ejk
  fixed sign error 23-jan-1998 ejk
  add 2-fold and 3-fold astigmatism to coherent mode 24-jan-1998 ejk
  add "\n" to one printf() format 28-jan-1998 ejk
  changed scaling in diffraction pattern mode 2-feb-1998 ejk
  fixed sqrt() in diff patt scaling 19-feb-1998 ejk
  changed variable PI in tcross() to pi because Linux gcc
      predefines PI  9-feb-1999 ejk
  update memory allocation routines 13-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  change data type of nxl,nyl to long32 for compatibility with new
      tiffsubs.c  17-jul-2007 ejk
  convert to GPL 4-jul-2008 ejk
  get return value of scanf() to remove warnings from gcc 4.4
      and convert to 4 char TAB size formatting 10-apr-2010 ejk
  convert to faster FFTW 10,24-apr-2010, 8-may-2010 ejk
  parameterize extra param[] offset so they can be moved easily
      and not conflict with slicelib.h offsets 21-dec-2010 ejk
  change Cs to Cs3 and add Cs5 on 21-dec-2010 ejk
  fix sign error that caused PC images to be flipped in tcross
       17-jan-2011 ejk
  convert to C++ and floatTIFF 21-mar-2012 ejk

  This file is formatted for a tab size of 4 characters

*/
#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

//#include "fftw3.h"         /* FFT routines from FFTW 3*/
//#include "slicelib.cpp"    /* misc. routines for multislice */
//#include "floatTIFF.hpp"   /* file I/O routines in TIFF format */

//const int NCMAX=   132; /* max number of characters per line */

const int MCPIX=   0;   /* coherent image mode */
const int MPCPIX=  1;   /* partial coherent mode */
const int MDIFFP=  2;   /* diffraction mode */

/* define more parameter keys not in slicelib.h */
#define pB      24     /* Debye-Waller factor * 0.25 in Ang */
#define pCS5    35     /* 5th order Cs5 */
#define p51 36
#define p52 37

#define p31 51  /* make symbolic offsets to merge with slicelib.h def's */
#define p32 52
#define p33 53
#define p34 54
#define p35 55
#define p36 56
#define p37 57
#define p38 58
#define p39 59


/*  define subroutines at end of file */
void invert2DW( fftwf_complex* pix, int nx, int ny );
void tcrossW( fftwf_complex* pixin, fftwf_complex* pixo,
        int nx, int ny, float *p );


int imageFromExitWave( char baseFilename[], const double &defoc, const double &a, const double &accel_pot,
										 const bool isMicro,
										const vector<vector<complex<double> > > &exitWave, 
										 vector<vector<double > > &TS, imagingArgs imaging )//  defocus is in metres
{

	My1TransferProjectedWF::Form1::TheInstance->progressBar1->Value = 0;
	if( defoc < 0 )
		My1TransferProjectedWF::Form1::TheInstance->progressBarDescription->Text = "Creating Underfocus Micrograph:";
	else if( defoc == 0 )
		My1TransferProjectedWF::Form1::TheInstance->progressBarDescription->Text = "Creating In-focus Micrograph:";
	else
		My1TransferProjectedWF::Form1::TheInstance->progressBarDescription->Text = "Creating Overfocus Micrograph:";
	
    char fileout[NCMAX];

	char filein[NCMAX];

	strcpy( filein, baseFilename );
	strcat( filein, ".tif" );

	strcpy( fileout, baseFilename );
	strcat( fileout, "(image)" );
	strcat( fileout, ".tif" );
	
    const char version[] = "16-apr-2012 (ejk)";

    int ix, iy, nx, ny, ixmid, iymid, npix, ns,
        itens, mode, nsum, NPARAM;
    int lcenter=0, lapert=0;

    float kx,ky, kx2,ky2, k2, k2max, v0, wavlen, scale, pixc,
        ax, by, cz, rx, ry, pi, rmin, rmax, aimin, aimax,
        Cs3, Cs5, df,  alpha0, ddf, objlx, objly,
        tr, ti, wr, wi, dfa2, dfa2phi, dfa3, dfa3phi;
    float *param;
    float  **pixr, **pixi;
    double sum, time, chi, chi1, chi2, chi3, phi, clog;

    fftwf_complex *cpix, *cpix2;  /* complex arrays for fft */
    fftwf_plan  planTf, planTi;   /* FFTW plans */

	bool isInputMicrographs = false;

    floatTIFF myFile;


	float aA = float( a*1e10 );

/*  echo version date */
/*
		printf( "image version dated %s\n", version );
		printf("Copyright (C) 1998-2012 Earl J. Kirkland\n" );
    printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
        " under the GNU general public license\n\n" );

    printf( "calculate TEM images with defocus, using FFTW\n\n");
	*/
	setcolor(0,15); cout << endl << "Calculate images..." << endl;setcolor(15,0);cout<< endl;
    pi = (float) (4.0 * atan( 1.0 ));

/*  get input file name */

    //printf("Name of file with input multislice result:\n");
    //ns = scanf("%s", filein );

/*  open input file and check sizes etc. */

	if( isInputMicrographs )
		if( myFile.tFloatTest( filein ) != 1 ) 
		{
        printf( "Cannot open input file %s\n", filein );
        exit( 0 );
		}

/*  ask mode */
/*
    printf("Type %d for coherent real space image,\n"
        "  or %d for partially coherent real space image,\n"
        "  or %d for diffraction pattern output:\n",
        MCPIX, MPCPIX, MDIFFP );
    ns = scanf("%d", &mode );
	*/


/*  get imaging parameters if required */

	mode = imaging.mode; // 0 for real space image, 1 for partially coherent, 2 for diffraction pattern
	Cs3 = imaging.Cs3;//1.2F;//Spherical abberation in mm
    Cs5 = imaging.Cs5;//Spherical abberation in mm
	df = float( defoc*1e10 ); //convert defocus to angstroms and float
    k2max = imaging.k2max; // objective aperture in mrad
	dfa2 = 0.0F; // magnitude of two-fold astigmatism (angstroms)
	dfa2phi = 0.0F; // angle of two-fold astigmatism (degrees)
	dfa3 = 0.0F; // magnitude of three-fold astigmatism (angstroms)
	dfa3phi = 0.0F; // angle of three-fold astigmatism (degrees)

    objlx = imaging.objlx;// "Objective lens and aperture center x,y in mrad (i.e. non-zero for dark field) [x-value]
    objly = imaging.objly;// "Objective lens and aperture center x,y in mrad (i.e. non-zero for dark field) [y-value]















    if( (mode == MCPIX) || (mode == MPCPIX) ) {
        //printf("Name of file to get defocused output:\n");
        //ns = scanf("%s", fileout );

        //printf("Spherical aberration Cs3, Cs5 (in mm.):\n");
        //ns = scanf("%f %f", &Cs3, &Cs5);
        Cs3 *= 1.0e7F;   /*  convert to Angstroms */
        Cs5 *= 1.0e7F;

        //printf("Defocus in Angstroms:\n");
        //ns = scanf("%f", &df );

        //printf("Objective aperture size in mrad:\n");
        //ns = scanf("%f", &k2max );

        if( mode == MPCPIX ) {
            printf("Illumination semi-angle in mrad:\n");
            ns = scanf("%f", &alpha0 );
            alpha0 = alpha0 * 0.001F;
            printf("Defocus spread in Angstroms:\n");
            ns = scanf("%f", &ddf );
        } else if( mode == MCPIX ) {
            /*printf( "Magnitude and angle of two-fold astigmatism"
                " (in Angst. and degrees):\n");
           ns = scanf( "%f %f", &dfa2, &dfa2phi); */
            dfa2phi = dfa2phi * pi /180.0F;
            /*printf( "Magnitude and angle of three-fold astigmatism"
                " (in Angst. and degrees):\n");
            ns = scanf( "%f %f", &dfa3, &dfa3phi);*/
            dfa3phi = dfa3phi * pi /180.0F;
            /*printf("Objective lens and aperture"
                " center x,y in mrad\n"
                " (i.e. non-zero for dark field):\n");
            ns = scanf("%f %f", &objlx, &objly);*/
        }

/*  get diffraction pattern parameters if required  */

    } else {
        printf("Name of file to get diffraction pattern:\n");
        ns = scanf("%s", fileout );

        lcenter = askYN("Do you want to include central beam");
        lapert = askYN("Do you want to impose the aperture");
        if ( lapert == 1 ) {
            printf("Aperture size in mrad:\n");
            ns = scanf("%f", &k2max );
            printf("Objective lens and aperture"
                " center x,y in mrad\n"
                "  (i.e. non-zero for dark field):\n");
            ns = scanf("%f %f", &objlx, &objly );
        }

        printf( "Type 0 for linear scale,\n"
            "  or 1 to do logarithmic intensity scale:\n"
            "  or 2 to do log(1+c*pixel) scale:\n");
        ns = scanf("%d", &itens );
        if( itens == 2 ) {
            printf( "Type scaling constant c:\n");
            ns = scanf( "%lf", &clog );
        }
    }  /* end else... */

    /* ---- read in specimen parameters and multislice result ---- */

    time = cputim();    /* get CPU time for comparison */

    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
    for( ix=0; ix<NPARAM; ix++) param[ix] = 0.0F;

   
	/*
	if( myFile.read( filein ) != 1 ) {
        printf( "Cannot open input file %s\n", filein );
        exit( 0 );
    }
	*/
    //nx = myFile.nx();
	nx = exitWave[0].size();

    //nx = nx/2;      /* should be complex */
    //ny = myFile.ny();
    ny = exitWave.size();
	//npix = myFile.getnpix();
    ixmid = nx/2;
    iymid = ny/2;

	/*
    if( npix != 2 ) {
        printf("Input file %s must be complex, can't continue.\n",
            filein );
        exit( 0 );
    }
	*/

    pixr = (float**) malloc2D( 2*nx, ny, sizeof(float), "pixr" );
    pixi = pixr + nx;

    //  crude port of old code using param[]
    for( ix=0; ix<NPARAM; ix++) param[ix] = myFile.getParam(ix);

	
    //ax = param[pDX] * ((float) nx);
    //by = param[pDY] * ((float) ny);
    ax = aA;
	by = aA;
	
	//if( param[pDY] <= 0.0F ) by = param[pDX] * ((float) ny);
    //cz = param[pC];
    cz = aA;

	rx  = 1.0F / ax;
    ry  = 1.0F / by;

	
    //v0 = param[pENERGY];
	v0 = float( accel_pot );

    wavlen = (float) wavelength( v0 );
    //printf("Starting pix energy = %8.2f keV\n", v0);
	
	/*
    rmax  = param[pRMAX];
    aimax = param[pIMAX];
    rmin  = param[pRMIN];
    aimin = param[pIMIN];
	
    printf("Starting pix range %g %g real\n"
           "                   %g %g imag\n", rmin,rmax,aimin,aimax);
*/
    param[pDEFOCUS] = df;
    param[pASTIG] = 0.0F;
    param[pTHETA] = 0.0F;
    param[pOAPERT] = k2max * 0.001F;
    param[pCS] = Cs3;
    param[pCS5] = Cs5;
    param[pWAVEL] = wavlen;

    k2max = k2max*0.001F/wavlen;
    k2max = k2max * k2max;

    objlx = objlx * 0.001F / wavlen;   /* convert to spatial freq. */
    objly = objly * 0.001F / wavlen;

    /*------  make FFTW arrays and plans ------- */
    cpix = (fftwf_complex*) fftwf_malloc( nx*ny * sizeof(fftwf_complex) );
    if( NULL == cpix ) {
        printf("Cannot allocate wave array cpix\n");
        exit( EXIT_FAILURE );
    }
    /* remember FFTW has inverse sign convention */
    planTi = fftwf_plan_dft_2d( nx, ny, cpix, cpix, 
        FFTW_FORWARD, FFTW_ESTIMATE );  /* inverse in place */
    planTf = fftwf_plan_dft_2d( nx, ny, cpix, cpix, 
        FFTW_BACKWARD, FFTW_ESTIMATE );  /* forward in place */

    /*------  copy to FFTW style array -------
        not very elegant but... */
    for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            //cpix[iy + ix*ny][0] = myFile(ix,iy);    // real
            //cpix[iy + ix*ny][1] = myFile(ix+nx,iy); // imag
			cpix[iy + ix*ny][0] = exitWave[iy][ix].real();    // real
            cpix[iy + ix*ny][1] = exitWave[iy][ix].imag(); // imag
    }

    /*  every option uses FT so do it here */
    fftwf_execute_dft( planTf, cpix, cpix );

/*------  coherent real space image mode ---------------------

    Convolute multislice result with objective lens
    aberration function (coherent case)
    with shifted objective lens for dark field

   NOTE zero freq. is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for FFT
*/

    if( mode == MCPIX ) {
        printf("calculate coherent image....\n");
        chi1  = pi * wavlen;
        chi2  = 0.5 * Cs3 * wavlen * wavlen;
        chi3  = Cs5 * wavlen*wavlen*wavlen*wavlen /3.0;

        ns = 0;
        for( ix=0; ix<nx; ix++) {
			My1TransferProjectedWF::Form1::TheInstance->progressBar1->Value = ix*100/(nx-1);
            kx = (float) ix;
            if( ix > ixmid ) kx = (float) (ix-nx);
            kx = kx*rx - objlx;
            kx2 = kx*kx;
            for( iy=0; iy<ny; iy++) {
                ky = (float) iy;
                if( iy > iymid ) ky = (float)(iy-ny);
                ky = ky*ry - objly;
                k2 = ky*ky + kx2;
                if ( k2 < k2max ) {
                    ns++;
                    phi = atan2( ky, kx );
                    chi = chi1*k2* ( (chi2 + chi3*k2)*k2 - df 
                        + dfa2*sin( 2.0*(phi-dfa2phi) ) 
                        + 2.0F*dfa3*wavlen*sqrt(k2)*
                        sin( 3.0*(phi-dfa3phi) )/3.0 );
                    tr = (float)  cos( chi );
                    ti = (float) -sin( chi );
                    wr = cpix[iy + ix*ny][0];
                    wi = cpix[iy + ix*ny][1];
                    cpix[iy + ix*ny][0] = tr*wr - ti*wi;  /* real */
                    cpix[iy + ix*ny][1] = tr*wi + ti*wr;  /* imag */
                 } else {
                    cpix[iy + ix*ny][0] = 0.0F;  /* real */
                    cpix[iy + ix*ny][1] = 0.0F;  /* imag */
                }
            }  /* end for( ix.... */
        }  /* end for(iy... */

		//Form1::TheInstance->imageConfirm->ForeColor = System::Drawing::Color::Lime;

        fftwf_execute_dft( planTi, cpix, cpix );
        scaleW( cpix, nx, ny );
        //printf("There were %d pixels inside the obj. apert.\n", ns );

        /*   copy back to old style array for output */
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            pixr[ix][iy] = cpix[iy + ix*ny][0];
            pixi[ix][iy] = cpix[iy + ix*ny][1];
        }

/* -------- calculate partially coherent image here  ----------------
    NOTE this assumes that the parameter offsets in txcoef()
    match those in slicelib.h !!!! (this is bad practice but....)
*/
    } else if( mode == MPCPIX ) {
        printf("calculate partially coherent image....\n");

        param[pCAPERT] = alpha0;
        param[pDDF] = ddf;
        param[ pB ] = 0.0F;

        /*------  make another FFTW array ------- */
        cpix2 = (fftwf_complex*) fftwf_malloc( nx*ny * sizeof(fftwf_complex) );
        if( NULL == cpix2 ) {
            printf("Cannot allocate wave array cpix2\n");
            exit( EXIT_FAILURE );
        }

        invert2DW( cpix, nx, ny );
        tcrossW( cpix, cpix2, nx, ny, param );
        invert2DW( cpix2, nx, ny );

        fftwf_execute_dft( planTi, cpix2, cpix2 );
        scaleW( cpix2, nx, ny );

        /*  copy back to old style arrays for output */
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            pixr[ix][iy] = cpix2[iy + ix*ny][0];
            pixi[ix][iy] = cpix2[iy + ix*ny][1];
        }

/*  do diffraction pattern here  */

    } else {

         printf("calculate diffraction pattern....\n");

        invert2DW( cpix, nx, ny );

       /*  copy back to old style arrays to use old subroutine */
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            pixr[ix][iy] = cpix[iy + ix*ny][0];  /* real */
            pixi[ix][iy] = cpix[iy + ix*ny][1];  /* imag */
        }

        if( lcenter == 0 ) 
            pixr[ixmid][iymid] = pixi[ixmid][iymid] = 0.0F;

        if ( lapert ) {
            for( iy=0; iy<ny; iy++) {
                ky = (iy-iymid)*ry - objly;
                ky2 = ky*ky;
                for( ix=0; ix<nx; ix++) {
                    kx = (ix-ixmid)*rx - objlx;
                    k2 = kx*kx + ky2;
                    if ( k2 > k2max )
                        pixr[ix][iy] = pixi[ix][iy] = 0.0F;
                } /* end for ix */
            } /* end for iy... */
        }  /* end if( lapert )... */

    }  /* end else... = last mode */


	
	
	/*Create vector to hold high resolution image*/
	vector<vector<double> > hires = CreateWF (  nx,ny, 0 ) ;
	
	/*  Output results and find min and max to echo */

    printf("output image\n");
    scale = 1.0F / ( ((float)nx) * ((float)ny) );
    
    myFile.resize( nx, ny );
    myFile.setnpix( 1 );
    sum = 0.0;
    nsum = 0;
    for( iy=0; iy<ny; iy++) {
        for( ix=0; ix<nx; ix++) {
            tr = pixr[ix][iy];
            ti = pixi[ix][iy];
            if( mode == MPCPIX ) pixc = tr;
            else  pixc = tr*tr + ti*ti;
            if ( (mode == MDIFFP) &&  (itens == 1)) {
                if( pixc > 1.e-30F)  pixc = (float) log( (double) pixc );
                else pixc = -30.0F;
            } else if ( (mode == MDIFFP) && (itens == 2)) {
                pixc = (float) log( 1.0 + clog*sqrt((double)pixc) );
            } 

			/* Normalise beam intensity to 1 */
			pixc /= 2;

            if( (ix == 0) && (iy == 0) ) {
                rmin = pixc;
                rmax = rmin;
            } else if( (ix != ixmid) && (iy != iymid) ) {
                if( pixc < rmin ) rmin = pixc;
                if( pixc > rmax ) rmax = pixc;
            }
            if( (ix>(3*nx)/8) && (ix<(5*nx)/8) &&
                (iy>(3*ny)/8) && (iy<(5*ny)/8) ) {
                sum = sum + pixc;
                nsum += 1;
            }
			


            myFile(ix,iy) = pixc;   
			/* Fill hires vector */
			hires[ix][iy] = double(pixc);
        }  /* end for ix... */

		
    } /* end for iy... */
	

    param[pRMAX] = rmax;
    param[pIMAX] = 0.0F;
    param[pRMIN] = rmin;
    param[pIMIN] = 0.0F;
    if ( mode == MDIFFP ) {
        param[pRMIN] = (float) (0.05*rmin + 0.95*sum/nsum);
        param[pDX] = 1.0F / (nx*param[pDX]);
        param[pDY] = 1.0F / (ny*param[pDY]);
    }

    //  crude port of old code using param[]
    for(ix=0; ix<NPARAM; ix++) myFile.setParam(ix, param[ix] );
	
	
	if( isMicro )
		myFile.write( fileout, rmin, rmax, 0.0F, 0.0F, 1, 1 );//myFile.write( fileout, rmin, rmax, 0.0F, 0.0F, param[pDX], param[pDY] );
		
		downsample( hires, TS );
		

		



   // printf("Pix range %f to %f\n",  param[pRMIN], rmax );

    time = cputim() - time;
    printf("Elapsed time = %f sec\n", time );


	/* Free dynamically allocated memory */
	free( pixr );
	fftwf_free( cpix );

	if( mode == MPCPIX )
		fftwf_free( cpix2 );

    return EXIT_SUCCESS;

} /* end main() */






void downsample( vector<vector<double> > hires, vector<vector<double> > &TS )

{
	size_t M1 = TS.size();
	size_t M2 = hires.size();
	size_t r = M2/M1;
	double sum;

	for( size_t i = 0 ; i < M1 ; i++ )
	{
		for( size_t j = 0 ; j < M1 ; j++ )
		{
			sum = 0;
			for( size_t t1 = 0 ; t1 < r ; t1++ )
			{
				for( size_t t2 = 0 ; t2 < r ; t2++ )
				{

					sum += hires[r*i+t1][r*j+t2];
				}
			}
			TS[i][j] = sum/(r*r);
		}
	}






}






/*------------------- tcrossW() ------------------------------*/
/*  this version uses FFTW style arrays

 Subroutines to do exact nonlinear partial coherence transfer as in
      [1]  M.A. O'Keefe, 37th EMSA (1979) p.556
      [2]  K. Ishizuka, Ultramicroscopy 5 (1980) p.55-65

  both input and output are in Fourier space with zero frequency
  in the center (NX/2 , NY/2)

 perform 2D weighted convolution with transmission cross coefficient
  TXCOEF to calculate partially coherent CTEM images
    NOTE: this version uses Friedels law

  rpixin, ipixin  : real,imag. input pix array dimensioned as [ix][iy]
             range used = P(17) = aperture
  rpixo, ipixo    : real,imag. output pix array dimensioned as [ix][iy]
             range output= 2.0*P(17) = 2 x aperture
  nx, ny : actual image size to work with
  p[43]  : parameter array ( dx,dy, defocus etc.) explained in TXCOEF

  started 14-NOV-85 on CIPRES2 Earl J. Kirkland
      started by modifying XFERFN.FTN (form RSX CIPRES)
  changed TXCOEF to be more general (i.e. include alignment and leading
      factors so that it is complete - little or no speed lose)
     28-NOV-85 EJK
  fixed typo (DT3 to D3T) and sign (-D4) in TXCOEF 17-jun-1986 ejk
  changed sign convention to be consistent with forward propagation
    of electrons (see Self et.al. UM 11 (1983) p.35-52 
      - only TXCOEF effected   EJK 2-JUL-86
  converted to C 8-nov-1995 ejk
  fixed sign error 23-jan-1998 ejk
  convert to FFTW style arrays 8-may-2010 ejk
  fix order of psi/psi* indexes to fix flipped image error
      17-jan-2011 ejk

*/

void tcrossW( fftwf_complex* pixin, fftwf_complex* pixo,
        int nx, int ny, float *p )
{
    int j, ix, iy, ixmid, iymid, ixmin, ixmax, iymin, iymax,
        ix2min, ix2max, iy0, iy2min, iy2max, ix1, ix2, iy1, iy2,
        *ixi, *ixf, ixi0, ixf0;
    float scale, k2p, k2pp, rx, ry, txcoefr, txcoefi,
        xr, xi, yr, yi, zr, zi, PI, *kxp, *kyp, *kxp2, *kyp2;

    void txcoef( float kxp, float kyp, float kxpp, float kypp,
         float p[], float *txcoefr, float *txcoefi );


/*  get scratch arrays and init params */

    ixi = (int*) malloc1D( nx, sizeof(int), "ixi" );
    ixf = (int*) malloc1D( nx, sizeof(int), "ixf" );
    kxp = (float*) malloc1D( nx, sizeof(float), "kxp" );
    kyp = (float*) malloc1D( ny, sizeof(float), "kyp" );
    kxp2 = (float*) malloc1D( nx, sizeof(float), "kxp" );
    kyp2 = (float*) malloc1D( ny, sizeof(float), "kyp2" );

    PI = (float) (4.0 * atan( 1.0 ));
    xr = p[pOAPERT]/p[pWAVEL];
    p[p31] = xr*xr;
    xi = p[pWAVEL];
    p[p32] = 0.5F*PI* p[pCS]* xi*xi*xi;
    p[p33] = PI* p[pWAVEL] * p[pDEFOCUS];

    p[p51] = PI * p[pCS5]* xi*xi*xi*xi*xi/3.0F;
    p[p52] = PI * p[pCS5]* xi*xi*xi*xi;

    xr = PI*p[pCAPERT]* p[pDDF];
    p[p34] = xr*xr;
    p[p35] = p[pCS]* p[pWAVEL]* p[pWAVEL];
    xr = PI*p[pCAPERT];
    p[p36] = xr*xr;
    xr =  PI* p[pWAVEL]* p[pDDF];
    p[p37] = xr*xr/4.0F;
    xr = PI*PI*p[pCAPERT]*p[pCAPERT]*p[pDDF] ;
    p[p38] = xr*xr;
    xr =  PI*p[pCAPERT]*p[pDDF] ;
    p[p39] = PI*( xr*xr )*p[pWAVEL];

    /* initialize  */

    rx = p[pDX]* nx;
    ry = p[pDY]* ny;
    if( ry <= 0.0F) ry = rx;
    rx = 1.0F/rx;
    ry = 1.0F/ry;

    scale=1.0F/( ((float)nx) * ((float)ny) );

    ixmid = nx/2;
    iymid = ny/2;

    /* find range of convolution */

    j = (int) ( p[pOAPERT] / (rx*p[pWAVEL]) +1.0F );
    ixmin  = ixmid - j;
    ixmax  = ixmid + j;
    ix2min = ixmid - 2*j;
    ix2max = ixmid + 2*j;
    j = (int) ( p[pOAPERT] / (ry*p[pWAVEL]) +1.0F );
    iymin  = iymid - j;
    iymax  = iymid + j;
    iy2min = iymid - 2*j;
    iy2max = iymid + 2*j;
    
    if( ixmin < 0 ) ixmin = 0;
    if( ixmax > (nx-1) ) ixmax = (nx-1);
    if( ix2min < 0 ) ix2min = 0;
    if( ix2max > (nx-1) ) ix2max = (nx-1);
    if( iymin < 0 ) iymin = 0;
    if( iymax > (ny-1) ) iymax = (ny-1);
    if( iy2min < 0 ) iy2min = 0;
    if( iy2max > (ny-1) ) iy2max = (ny-1);

    for ( ix=ix2min; ix<=ix2max; ix++) {
       ix1= ixmin - ix + ixmid;
       if( ixmin > ix1 ) ixi[ix] = ixmin;  else  ixi[ix] = ix1;
       ix1= ixmax - ix + ixmid;
       if( ixmax < ix1 ) ixf[ix] = ixmax;   else  ixf[ix] = ix1;
    }

    for( ix=ix2min; ix<=ix2max; ix++) {
        kxp[ix] = (ix-ixmid) * rx;
        kxp2[ix]= kxp[ix] * kxp[ix];
    }
    for( iy=iy2min; iy<=iy2max; iy++) {
        kyp[iy] = (iy-iymid) * ry;
        kyp2[iy]= kyp[iy] * kyp[iy];
    }

    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
        pixo[iy + ix*ny][0] = pixo[iy + ix*ny][1] = 0.0F;
/*???? 13-jan-2011        pixo[iy + ix*ny][0] = pixo[iy + ix*ny][0] = 0.0F; */

/*  do convolution in bottom half plane 
    the origin is at (ixmid,iymid)
    (ix,iy)   = output point = k
    (ix1,iy1) = input point from psi* = kp
    (ix2,iy2) = input point from psi  = kp + k
*/

    for( iy=iy2min; iy<=iymid; iy++ )
        for( iy1=iymin; iy1<=iymax; iy1++ ) {
            iy2 = iy1 + (iy - iymid);
            if( (iy2>=iymin) && (iy2<=iymax)) {
                for( ix=ix2min; ix<=ix2max; ix++ ) {
                    ixi0 = ixi[ix];
                    ixf0 = ixf[ix];
                    for( ix1=ixi0; ix1<=ixf0; ix1++ ) {
                        ix2= ix1 + (ix - ixmid);
                        k2p = kxp2[ix1] + kyp2[iy1];
                        k2pp= kxp2[ix2] + kyp2[iy2];
                        if ( (k2p<=p[p31]) && (k2pp<=p[p31]) ) {
                            txcoef( kxp[ix1],kyp[iy1],kxp[ix2],kyp[iy2],
                                p, &txcoefr, &txcoefi);
                            xr = pixin[iy1 + ix1*ny][0];  /* rpixin[ix1][iy1]; */
                            xi = -pixin[iy1 + ix1*ny][1]; /* -ipixin[ix1][iy1]; */
                            yr = pixin[iy2 + ix2*ny][0];  /* rpixin[ix2][iy2]; */
                            yi = pixin[iy2 + ix2*ny][1];  /* ipixin[ix2][iy2]; */
                            zr = xr*yr - xi*yi;
                            zi = xr*yi + xi*yr;
                            pixo[iy + ix*ny][0] += zr*txcoefr - zi*txcoefi;
                            pixo[iy + ix*ny][1] += zr*txcoefi + zi*txcoefr;
                        }  /* end if( ( k2p... */
                    }  /* end for( ix1=... */
                }  /* end for( ix=... */
            }  /* end if( (iy2... */
        }  /* end for iy1... */

/*  scale result */

    for( ix= ix2min; ix<=ix2max; ix++) 
    for( iy= iy2min; iy<=iymid; iy++) {
        pixo[iy + ix*ny][0] *= scale;
        pixo[iy + ix*ny][1] *= scale;
    }

/*  Invoke Friedel's law to get top half plane  */

    iy0 = iymid+1;
    for( iy=iy0; iy<=iy2max; iy++) {
        iy1= iymid - (iy-iymid);
        for( ix=ix2min; ix<=ix2max; ix++) {
            ix1= ixmid - (ix-ixmid);
            pixo[iy + ix*ny][0]=  pixo[iy1+ix1*ny][0]; /* rpixo[ix1][iy1]; */
            pixo[iy + ix*ny][1]= -pixo[iy1+ix1*ny][1]; /* -ipixo[ix1][iy1]; */
        }
       pixo[iy][0] = pixo[iy][1] = 0.0F;  /* ??? */
    }

    free( ixi );
    free( ixf );
    free( kxp );
    free( kyp );
    free( kxp2 );
    free( kyp2 );

    return;
}     /*   end tcrossW() */

/*------------------------ txcoef -------------------------*/
/*
   The cross spectral transfer function as in:
    [1] M.A. O'Keefe, 37th EMSA (1979) p.556 EJK
    [2] K. Ishizuka, Ultramic. 5 (1980) p.55-65.

   switched to my derivation 5-DEC-83 EJK
   changed sign convention to be consistent with forward propagation
     of electrons (see Self et.al. UM 11 (1983) p.35-52 
       EJK 2-JUL-86
   converted to C  6-nov-1995 ejk
   fixed sign error 23-jan-1998 ejk
   fix sign to fix flipped image error
      17-jan-2011 ejk

  the following are the array index definitions
  (the order is historical - don't ask why! )

   P(11) = defocus
   P(12) = defocus astigmatism
   P(13) = angle of astigmatism
   P(14) = dx in A pixel dimension
   P(15) = dy in A pixel dimension
   P(17) = aperture in radians
   P(18) = Cs3 in A
   P(19) = wavelength in A
   P(21) = illumination semiangle in rad
   P(22) = defocus spread in A
   P(24) = Debye Waller temp factor/4 in A
   P(31) = P(17)/P(19) **2 maximum k^2 value
 
   P(32) = PI*P(18)* (P(19)**3) /2    ; coherent part
   P(33) = PI* P(19) * P(11)

   p(pCS5) = Cs5 in A         ; add 5th order 21-dec-2010 ejk
   p(p51) = PI*P(pCS5)* (P(19)**5) /3.0 
 
   P(34) = (PI* P(21)* P(22))**2
   P(35) = P(18)* P(19)* P(19)
   P(36) = (PI* P(21)) **2
   P(37) = (PI* P(19)* P(22))**2 /4.
   P(38) = PI^4 * P(21)^4 *P(22)^2
   P(39) = PI^3 * P(22)^2 * P(21)^2 *P(19)
      
   kp  = kprime      = argument of psi*
   kpp = kprime + k  = argument of psi
 
   txcoefr, txcoefi = returned real and imag. parts of function
        (i.e. C can't return a complex value )
*/

void txcoef( float kxp, float kyp, float kxpp, float kypp,
         float p[], float *txcoefr, float *txcoefi )
{
    double kx, ky, k2 ,kp2, kpp2, chip, chipp,
        d1, d2, d3, d3t, d4, d4t, vx, vy, w, u, xr;

    kp2  = kxp*kxp + kyp*kyp;
    kpp2 = kxpp*kxpp + kypp*kypp;

    if( ( kp2 > p[p31] ) || ( kpp2 > p[p31] ) ) {
      *txcoefr = *txcoefi = 0.0F;
      return;
    }

    kx = kxpp - kxp;   /* output freq. k */
    ky = kypp - kyp;
    k2 = kx*kx + ky*ky;

    /*  v = Wc1/(2*pi*lambda)
        w = Wc2/(-pi*lambda)
        u = 1 + pi^2 * beta^2 * delta0^2 k^2
    */
/* -- add Cs5 below 21-dec-2010 ejk */
    vx = p[p52]*( kp2*kp2*kxp - kpp2*kpp2*kxpp )
        + p[p35]*( kp2*kxp - kpp2*kxpp ) + p[pDEFOCUS]*kx;
    vy = p[p52]*( kp2*kp2*kyp - kpp2*kpp2*kypp )
        + p[p35]*( kp2*kyp - kpp2*kypp ) + p[pDEFOCUS]*ky;
    w  = kp2 - kpp2;
    u  = 1.0 + p[p34]*k2;

    d1  = p[p36] * (vx*vx+vy*vy);
    d2  = p[p37] *w*w /u;
    d3t = vx*kx + vy*ky;
    d3  = p[p38]*d3t*d3t/u;
    d4t = p[p39]*w/u;
    d4  = d4t*d3t;

/* ---- coherent chi  --------------------
     add 5th order Cs5 21-dec-2010 ejk */
    chip  = ( (p[p51]*kp2  + p[p32] )*kp2  - p[p33] ) *kp2;
    chipp = ( (p[p51]*kpp2 + p[p32] )*kpp2 - p[p33] ) *kpp2;

    chip = chip - chipp - d4;
    xr = exp( -d1 -d2 +d3 -p[pB]*(kp2+kpp2) ) / sqrt(u);
    *txcoefr = (float) ( cos(chip) * xr );
    *txcoefi = (float) ( sin(chip) * xr );

    return;

}  /* end txcoef() */


/*------------------------- invert2DW() ----------------------*/
/*  from fft2dc.c and converted to FFTW style arrays
     6-may-2010 ejk

        rearrange pix with corners moved to center (a la FFT's)

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void invert2DW( fftwf_complex* pix, int nx, int ny )
{
#define SWAP(a,b)       {t=a; a=b; b=t;}

    int ix, iy, ixmid, iymid;
    float t;

    ixmid = nx/2;
    iymid = ny/2;

    for( ix=0; ix<nx; ix++) 
    for( iy=0; iy<iymid; iy++) {
                SWAP( pix[iy + ix*ny][0], pix[iy+iymid + ix*ny][0] );
                SWAP( pix[iy + ix*ny][1], pix[iy+iymid + ix*ny][1] );
    }

    for( ix=0; ix<ixmid; ix++) 
    for( iy=0; iy<ny; iy++) {
                SWAP( pix[iy+ix*ny][0], pix[iy+(ix+ixmid)*ny][0] );
                SWAP( pix[iy+ix*ny][1], pix[iy+(ix+ixmid)*ny][1] );
    }

#undef SWAP
}  /* end invert2DW()  */

