#include "mex.h"
#include <stdio.h>

/* Generates a report of features and segmentation boundaries for ingestion into a database. */
/* Format is slide \t globalcentroid-x \t globalcentroid-y \t {features} \t globalx1,globaly1; globalx2,globaly2... */

/* input: */
/* Output - (string) Path and filename of output file */
/* Slide - (string) Base slide name (without extensions) for generating output. */

/* X - (N-length float) Horizontal coordinates of segmented object centroid in global frame. */
/* Y - (N-length float) Vertical coordinates of segmented object centroid in global frame. */
/* Features - (N x D float) Features of the segmented objects. */
/* Names - (D-length cell array of strings) Feature names. */
/* bX - (N-length cell array of float) Global horizontal coordinates of segmented object boundaries. */
/*          Each element is a K_i-length vector. */
/*          Bounds should be filtered to first to remove colinear points and compress size. */
/* bY - (N-length cell array of float) Global vertical coordinates of segmented object boundaries. */
/*          Each element is a K_i-length vector. */
/*          Bounds should be filtered to first to remove colinear points and compress size. */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* variables */
    FILE *file;
    int i, j; /* loop iterators */
    char *Output, *Slide, *Name; /*pointer to input string denoting slide*/
    double *X, *Y, *Features, *bX, *bY; /* pointers to input matrices */
    mxArray *Names, *BoundsX, *BoundsY; /* cell array of feature names */
    size_t N, D, Length; /* number of objects, feature dimension, # boundary points */

    /* capture pointers to input arrays */
    X = (double*)mxGetPr(prhs[2]);
    Y = (double*)mxGetPr(prhs[3]);
    Features = (double*)mxGetPr(prhs[4]);
    Names = prhs[5];
    BoundsX = prhs[6];
    BoundsY = prhs[7];

    /* get number of objects */
    if(mxGetM(prhs[2]) > mxGetN(prhs[2])) {
        N = mxGetM(prhs[2]);
    }
    else {
        N = mxGetN(prhs[2]);
    }
    
    /* get feature dimension */
    if(mxGetM(prhs[4]) != N) {
        D = mxGetM(prhs[4]);
    }
    else {
        D = mxGetN(prhs[4]);
    }

    /* copy strings for slide name and output filename */
    Output = mxCalloc(mxGetM(prhs[0])*mxGetN(prhs[0])+1, sizeof(char));
    Slide = mxCalloc(mxGetM(prhs[1])*mxGetN(prhs[1])+1, sizeof(char));
    mxGetString(prhs[0], Output, mxGetM(prhs[0])*mxGetN(prhs[0])+1);
    mxGetString(prhs[1], Slide, mxGetM(prhs[1])*mxGetN(prhs[1])+1);    

    /* open file, write contents */
    file = fopen(Output, "w");
    if(file != NULL) {

        /* print header */
        fprintf(file, "Slide\tX\tY\t");
        for(i = 0; i < D; i++) {
            mxArray *Cell = mxGetCell(Names,i);
            Name = mxCalloc(mxGetM(Cell)*mxGetN(Cell)+1, sizeof(char));
            mxGetString(Cell, Name, mxGetM(Cell)*mxGetN(Cell)+1);
            fprintf(file, "%s\t", Name);
            mxFree(Name);
        }
        fprintf(file, "Boundaries: x1,y1;x2,y2;...\n");

        /* print data */
        for(i = 0; i < N; i++){
            
            /* slide name, centroids */
            fprintf(file, "%s\t%.1f\t%.1f\t", Slide, X[i], Y[i]);
            
            /* features */
            for(j = 0; j < D; j++){
                fprintf(file, "%f\t", Features[i+ j*N]);
            }
            
            /* boundaries */
            bX = mxGetPr(mxGetCell(BoundsX,i));
            bY = mxGetPr(mxGetCell(BoundsY,i));
            Length = mxGetM(mxGetCell(BoundsX,i)) * mxGetN(mxGetCell(BoundsX,i));
            for(j = 0; j < Length; j++) {
                fprintf(file, "%.1f,%.1f;", bX[j], bY[j]);
            }
            fprintf(file, "\n");
        }
    }

    /* close file */
    fclose(file);
    
    /* free mem */
    mxFree(Output);
    mxFree(Slide);

}
