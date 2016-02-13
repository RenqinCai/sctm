#include "infer.h"

void infer (char* odir, sctm_data* data, sctm_params* params,
		sctm_latent* latent, sctm_counts* counts) {
	int iter, d;
	double* docLogLikelihood = (double*) malloc(sizeof(double)*(data->D));
	int token = 0;
	int likelihoodCount; 

	likelihoodCount = 0;

	printf("\ninfer---likelihoodcount%d", likelihoodCount);

	for(d=0; d<data->D; d++){
		*(docLogLikelihood+d) = 0;
	}

	for (iter = 1; iter < params->ITER+1; iter++) {

		if (params->word_sparsity > 0 && params->trte == 0)
			infer_phi(data, params, latent, counts);

		if (params->model == 2) {
			infer_b(data, params, latent, counts);
		}

		infer_z(data, params, latent, counts, iter, odir);
		//assignment(odir, data, params, latent, counts, iter);

		if (params->CMNTS) {
			if (params->sents_sparsity) infer_xi(data, params, latent, counts);
			//assignment(odir, data, params, latent, counts, iter);

			infer_y(data, params, latent, counts, iter);
			//assignment(odir, data, params, latent, counts, iter);

			infer_t(data, params, latent, counts);
		}

		// printf("end infer....");
		if (iter % params->save_step == 0 || iter == 1 || (iter >= params->burn_in && (iter-params->burn_in)%params->save_state == 0)) {
			// printf("\nsave....");
			// fflush(stdout);

			// printf("\n%4diterations", iter);
			// printf("\n%dtotal iterations", params->ITER);
			printf("\n%4d of %d iterations done", iter, params->ITER);
			
			// if(params->trte==1){

			// 	// printf("compute_likelihood");
			// 	// fflush(stdout);

				likelihoodCount += 1;
				compute_likelihood(data, params, latent, counts, docLogLikelihood, &token);
				printf("\nlikelihoodcount%d", likelihoodCount);
				fflush(stdout);

			// }

			assignment(odir, data, params, latent, counts, iter);
			fflush(stdout);
//			if (params->trte == 1)
//				compute_perplexity(odir, data, params, latent, counts);
		}

	}

	assignment(odir, data, params, latent, counts, iter);
	if (params->trte == 1){


		// likelihoodCount += 1;
		// compute_likelihood(data, params, latent, counts, docLogLikelihood, &token);
		// printf("\nlikelihoodcount%d", likelihoodCount);
		// fflush(stdout);

			// }
		printf("\nlikelihoodCount%d", likelihoodCount);
		compute_perplexity(odir, data, docLogLikelihood, &likelihoodCount);
	}


	// if (params->trte == 1)
	// 	compute_perplexity(odir, data, params, latent, counts, docLogLikelihood);


}

void infer_phi (sctm_data* data, sctm_params* params, sctm_latent* latent,
		sctm_counts* counts) {
	int k, v;
	double p1, p2, lognorm, r;

	for (k = 0; k < params->K; k++) {
		for (v = 0; v < data->V; v++) {
			counts->phiEta_v[k] -= latent->phi[k][v] * params->eta;
			counts->phi_k[v] -= latent->phi[k][v];

			if (counts->n_dij[k][v] > 0)
				latent->phi[k][v] = 1;
			else {
				p1 = lgamma(counts->phiEta_v[k] + params->eta)
						- lgamma(counts->phiEta_v[k] + params->eta
										+ counts->n_dijv[k]);

				p1 += log(params->nu_lambda / (double) data->V + counts->phi_k[v])
						- log(params->nu_lambda / (double) data->V
										+ params->lambda + params->K);

				if (counts->phiEta_v[k] > 0) {
					p2 = lgamma(counts->phiEta_v[k])
						- lgamma(counts->phiEta_v[k] + counts->n_dijv[k]);

					p2 += log(params->lambda + params->K - counts->phi_k[v] - 1)
						- log(params->nu_lambda / (double) data->V
										+ params->lambda + params->K);
					lognorm = log_sum(p1, p2);
				} else lognorm = p1;			

				r = myrand();

				if (log(r) + lognorm <= p1)
					latent->phi[k][v] = 1;
				else
					latent->phi[k][v] = 0;
			}

			counts->phiEta_v[k] += latent->phi[k][v] * params->eta;
			counts->phi_k[v] += latent->phi[k][v];
		}
	}
}
