#ifndef ASSIGN_H
#define ASSIGN_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include "sctm.h"
#include "utils.h"

void assignment(char* odir, sctm_data* cdata, sctm_params* params,
		sctm_latent* latent, sctm_counts* counts, int iter);

// double compute_perplexity(char* odir, sctm_data* cdata, sctm_params* params,
// 		sctm_latent* latent, sctm_counts* counts);
// double compute_perplexity(char* odir, sctm_data* cdata, sctm_params* params, sctm_latent* latent, sctm_counts* counts, double *docLogLikelihood);

// double compute_likelihood(sctm_data* cdata, sctm_params* params,
// 		sctm_latent* latent, sctm_counts* counts, double* result, int* tokens);

// double compute_likelihood_cmnt(sctm_data* cdata, sctm_params* params,
// 		sctm_latent* latent, sctm_counts* counts, double* result, int* tokens);

double compute_perplexity(char* odir, sctm_data* data, double *docLogLikelihood, int *count);

double compute_likelihood(sctm_data* data, sctm_params * params, sctm_latent* latent, sctm_counts* counts, double *docLogLikelihood, int *token);

void retrieveChild4Stn(sctm_data* data, sctm_params*params, sctm_latent*latent, sctm_counts* counts, double ***stnLikelihood);

void printChild4Stn(char* odir, sctm_data* data, sctm_params*params, sctm_latent*latent, sctm_counts* counts);

#endif
