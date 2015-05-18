#ifndef INFER_B_H
#define INFER_B_H

#include "sctm.h"
#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <float.h>

void infer_b(sctm_data* data, sctm_params* params, sctm_latent* latent,
		sctm_counts* counts);

#endif
