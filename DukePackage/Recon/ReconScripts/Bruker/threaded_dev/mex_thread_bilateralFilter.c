#include "mex.h"
#include "math.h"
#include <pthread.h>
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))


    unsigned int nVox, maxNeighbors;
    double inv_intensity_sigma_sqr;
    mwSize nDims;
    double *unfiltVol, *filtVol, *d_proximity;
    double *filt_radius, *proximity_sigma, *intensity_sigma;
    double *inv_proximity_sigma_sqr;
    const mwSize *dims;
    unsigned int *idx_convert;
    signed int *cur_nhood_idx, *neigh_shift;
   
struct ThreadData {
    int startVox, endVox;
};
    
void* filterVoxels(void *td){
    /* Variable declarations */
    unsigned int iVox, iNeighbor;
    signed int neighVox; //must be signed
    double d, kernel_sum;
    struct ThreadData *threadData = td;
    
    // Loop through each voxel in input volume
    for(iVox=threadData->startVox;iVox<=threadData->endVox;iVox++){
        //Loop through neighborhood
        kernel_sum = 0;
        filtVol[iVox] = 0;
        for(iNeighbor=0;iNeighbor<maxNeighbors;iNeighbor++){
            neighVox = iVox + neigh_shift[iNeighbor];
            if(neighVox>=0 & neighVox<nVox){
                //Calculate intensity weight
                d = (unfiltVol[iVox]-unfiltVol[neighVox]);
                d = d_proximity[iNeighbor]*exp(d*d*inv_intensity_sigma_sqr);
                kernel_sum = kernel_sum + d;
                
                // Accumulate voxel intensity
                filtVol[iVox] = filtVol[iVox] + unfiltVol[neighVox]*d;
            }
        }
        
        //Normalize voxel by kernel sum
        filtVol[iVox] = filtVol[iVox]/kernel_sum;
    }
    return NULL;
}
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    /* Variable declarations */
    unsigned int iNeighbor;
    signed int dim;
    unsigned int *nthreads;
    struct ThreadData *threadData;
    unsigned int iThread, numPtsPerThread;
    pthread_t *thread;
    
    printf("Startind mex...\n");
    
    /* INPUT 0 - UNFILTERED VOLUME */
    dims =  mxGetDimensions(prhs[0]);
    nDims = mxGetNumberOfDimensions(prhs[0]);
    unfiltVol = mxGetPr(prhs[0]);
    
    /* INPUT 1 - FILTERED VOLUME */
    filtVol = mxGetPr(prhs[1]);
    
    /* INPUT 2 - FILTER RADIUS */
    filt_radius = mxGetPr(prhs[2]);
    
    /* INPUT 3 - PROXIMITY SIGMA */
    proximity_sigma = mxGetPr(prhs[3]);
    
    /* INPUT 4 - INTENSITY SIGMA */
    intensity_sigma = mxGetPr(prhs[4]);
    
    /* INPUT 5 - NTHREADS */
    nthreads = (unsigned int*)mxGetPr(prhs[5]);
    threadData = calloc(*nthreads,sizeof(struct ThreadData));
    thread = calloc(*nthreads,sizeof(pthread_t));
    
    // Allocate memory for temp variables
    cur_nhood_idx = calloc(nDims, sizeof(signed int));
    idx_convert = calloc(nDims, sizeof(unsigned int));
    inv_proximity_sigma_sqr = calloc(nDims, sizeof(double));
    
    //Calculate number of voxels and set starting voxel to first voxel
    nVox = dims[0];
    idx_convert[0] = 1;
    maxNeighbors = floor(2*filt_radius[0]+1);
    cur_nhood_idx[0] = -filt_radius[0];
    for(dim=1; dim<nDims; dim++){
        nVox = nVox * dims[dim];
        idx_convert[dim] = idx_convert[dim-1]*dims[dim-1];
        maxNeighbors = maxNeighbors*floor(2*filt_radius[dim]+1);
        cur_nhood_idx[dim] = -filt_radius[dim];
        inv_proximity_sigma_sqr[dim] = -0.5/(proximity_sigma[dim]*proximity_sigma[dim]);
    }
    inv_intensity_sigma_sqr = -0.5/(intensity_sigma[0]*intensity_sigma[0]);
    
    /* allocate */
    d_proximity = calloc(maxNeighbors, sizeof(double));
    neigh_shift = calloc(maxNeighbors, sizeof(signed int));
    for(iNeighbor=0; iNeighbor<maxNeighbors; iNeighbor++){
        d_proximity[iNeighbor] = 0;
        neigh_shift[iNeighbor] = 0;
        for(dim=0; dim<nDims; dim++){
            d_proximity[iNeighbor] = d_proximity[iNeighbor] + 
                    (cur_nhood_idx[dim]*cur_nhood_idx[dim]*inv_proximity_sigma_sqr[dim]);
            neigh_shift[iNeighbor] = neigh_shift[iNeighbor] + cur_nhood_idx[dim]*idx_convert[dim];
        }
        d_proximity[iNeighbor] = exp(d_proximity[iNeighbor]);
        
        //Update neighbor inxex for next neighbor
        for(dim=0; dim<nDims; dim++){
            if(cur_nhood_idx[dim]==(filt_radius[dim])){
                //Reset that dimmension, then the next dimmension will increment
                cur_nhood_idx[dim] = -filt_radius[dim];
            }else{
                cur_nhood_idx[dim]++;
                break;
            }
        }
    }
    
    numPtsPerThread = floor(((float)nVox)/(*nthreads));
    for(iThread=0; iThread<*nthreads; iThread++){
        threadData[iThread].startVox = iThread*numPtsPerThread;
        threadData[iThread].endVox = (iThread+1)*numPtsPerThread-1;
    }
    threadData[*nthreads-1].endVox = nVox-1;
    
    printf("Startind Threads...\n");
    for(iThread=0; iThread<*nthreads; iThread++){
        pthread_create(&thread[iThread], NULL, filterVoxels, &threadData[iThread]);
    }
    
    /* Wait for threads to finish */
    for(iThread=0; iThread<*nthreads; iThread++){
        pthread_join(thread[iThread],NULL);
    }
    printf("Finished mex...\n");
        
    // Free up temp variables
    free(cur_nhood_idx);
    free(idx_convert);
    free(d_proximity);
    free(neigh_shift);
    free(inv_proximity_sigma_sqr);
    free(threadData);
    free(thread);
}


