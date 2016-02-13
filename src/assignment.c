#include "assignment.h"
#include <math.h>

void assignment(char* odir, sctm_data* data, sctm_params* params,
		sctm_latent* latent, sctm_counts* counts, int iter) {

	int d, i, n, k, a, c, j, v, sum, state;

	char fname[500];
	char dir[1000];

	FILE *fbeta=NULL, *fz_dist;
//	FILE *fz, *fb , *fphi;
//	FILE *fm, *fm_k, *fm_1, *fm_1k;

	sprintf(dir, "%s", odir);

//	sprintf(fname, "%s/b.txt", dir);
//	fb = fopen(fname, "w");
	if (params->trte == 0) {
//		sprintf(fname, "%s/phi.txt", dir);
//		fphi = fopen(fname, "w");
		sprintf(fname, "%s/beta", dir);
		fbeta = fopen(fname, "w");
	}


	if (iter >= params->burn_in && (iter-params->burn_in)%params->save_state == 0) {
		state = (iter-params->burn_in)/params->save_state;
	}
	else state = 0;

//	sprintf(fname, "%s/z.txt", dir);
//	fz = fopen(fname, "w");
	if (params->trte == 0) sprintf(fname, "%s/z_dist.txt", dir);
	else sprintf(fname, "%s/z_dist_test.txt", dir);
	fz_dist = fopen(fname, "w");
	
	double* z_dist = (double*) malloc(sizeof(double)*params->K);
	for (k=0; k < params->K; k++) z_dist[k] = 0.;

	for (d = 0; d < data->D; d++) {
		documents* doc = &(data->docs[d]);

//		fprintf(fz, "%d %d\n", d + 1, doc->S);
//		fprintf(fb, "%d %d\n", d + 1, doc->S);
		fprintf(fz_dist, "%d %d\n", d + 1, doc->S);

		for (i = 0; i < doc->S; i++) {
			sentence* sent = &(doc->sents[i]);
			for (n = 0; n < sent->N; n++) {
				k = latent->z[d][i][n];
				j = latent->b[d][i][n];
//				fprintf(fz, "%d ", k);
//				fprintf(fb, "%d ", j);
				z_dist[k] += 1;
			}
//			fprintf(fz, "\n");
//			fprintf(fb, "\n");
			for (k=0; k < params->K; k++) {
				if (state > 0 || iter == params->burn_in) {
					latent->z_dist[d][i][k] = (state*latent->z_dist[d][i][k] + z_dist[k]/sent->N)/(state+1);
					fprintf(fz_dist, "%.4f ", latent->z_dist[d][i][k]);
				}
				else if (iter > params->ITER) fprintf(fz_dist, "%.4f ", latent->z_dist[d][i][k]);
				else fprintf(fz_dist, "%.4f ", z_dist[k]/sent->N);
				z_dist[k] = 0.; //for next sentence
			}
			fprintf(fz_dist, "\n");
		}
//		fprintf(fz, "\n");
//		fprintf(fb, "\n");
		fprintf(fz_dist, "\n");
	}
//	fclose(fz);
//	fclose(fb);
	fclose(fz_dist);
	free(z_dist);

//	FILE *ft = NULL, *fxi = NULL, *fy = NULL;
	FILE *fxi_prob = NULL, *fy_dist = NULL;
	if (params->CMNTS) {
//		sprintf(fname, "%s/t.txt", dir);
//		ft = fopen(fname, "w");
//		sprintf(fname, "%s/xi.txt", dir);
//		fxi = fopen(fname, "w");
		if (params->trte == 0) sprintf(fname, "%s/xi_prob.txt", dir);
		else sprintf(fname, "%s/xi_prob_test.txt", dir);
		fxi_prob = fopen(fname, "w");
//		sprintf(fname, "%s/y.txt", dir);
//		fy = fopen(fname, "w");
		if (params->trte == 0) sprintf(fname, "%s/y_dist.txt", dir);
		else sprintf(fname, "%s/y_dist_test.txt", dir);
		fy_dist = fopen(fname, "w");

		double* y_dist = (double*) malloc(sizeof(double)*(params->K+1));
		for (k=0; k < params->K+1; k++) y_dist[k] = 0.;

		for (d = 0; d < data->D; d++) {
			documents* doc = &(data->docs[d]);
//			fprintf(fxi, "%d %d\n", d + 1, doc->C);
//			fprintf(ft, "%d %d\n", d + 1, doc->C);
//			fprintf(fy, "%d %d\n", d + 1, doc->C);
			fprintf(fxi_prob, "%d %d\n", d + 1, doc->C);
			fprintf(fy_dist, "%d %d\n", d + 1, doc->C);
			for (c = 0; c < doc->C; c++) {
				comment* cmnt = &(doc->cmnts[c]);
				for (a = 0; a < doc->S; a++) {
//					if (latent->xi[d][c][a])
//						fprintf(fxi, "%d ", a);
					fprintf(fxi_prob, "%.3f ", latent->xi_prob[d][c][a]);
				}
//				fprintf(fxi, "\n");
				fprintf(fxi_prob, "\n");
				for (n = 0; n < cmnt->N; n++) {
//					fprintf(ft, "%d ", latent->t[d][c][n]);
//					fprintf(fy, "%d ", latent->y[d][c][n]);
					k = latent->y[d][c][n];
					if (latent->t[d][c][n] == 0) k = params->K;
					y_dist[k] += 1;
				}
//				fprintf(ft, "\n");
//				fprintf(fy, "\n");
				for (k=0; k < params->K+(int)(params->model!=0); k++) {
					if (state > 0 || iter == params->burn_in) {
						latent->y_dist[d][c][k] = (state*latent->y_dist[d][c][k] + y_dist[k]/cmnt->N)/(state+1);
						fprintf(fy_dist, "%.4f ", latent->y_dist[d][c][k]);
					}
					else if (iter > params->ITER) fprintf(fy_dist, "%.4f ", latent->y_dist[d][c][k]);
					else fprintf(fy_dist, "%.4f ", y_dist[k]/cmnt->N);
					y_dist[k] = 0.; //for next cmnt
				}
				fprintf(fy_dist, "\n");
			}
//			fprintf(fxi, "\n");
			fprintf(fxi_prob, "\n");
//			fprintf(ft, "\n");
//			fprintf(fy, "\n");
			fprintf(fy_dist, "\n");
		}
//		fclose(fxi);
		fclose(fxi_prob);
//		fclose(ft);
//		fclose(fy);
		fclose(fy_dist);

		free(y_dist);

	}

	if (params->trte == 0) {
		sum = 0;
		double *beta;
		beta = (double*) malloc(sizeof(double)*data->V);
		if (params->model != 0) fprintf(fbeta, "%d %d\n", params->K+1, data->V);
		else fprintf(fbeta, "%d %d\n", params->K, data->V);
		for (k = 0; k < params->K; k++) {
			sum = 0;
			for (v = 0; v < data->V; v++) {
				if (latent->phi[k][v])
					sum += 1;
				beta[v] = (latent->phi[k][v] * params->eta + counts->n_dij[k][v])*1.0/ (counts->phiEta_v[k] + counts->n_dijv[k]);
				if (state > 0 || iter == params->burn_in) {
					latent->beta[k][v] = (state*latent->beta[k][v] + beta[v])/(state+1);
					fprintf(fbeta, "%lf ", latent->beta[k][v]);
				}
				else if (iter > params->ITER) fprintf(fbeta, "%lf ", latent->beta[k][v]);
				else fprintf(fbeta, "%lf ", beta[v]);
			}
//			fprintf(fphi, "%d\n", sum);
			fprintf(fbeta,"\n");
		}
		for (v = 0; v < data->V; v++) {
			beta[v] = (params->eta + counts->n_dij[k][v])*1.0/ (data->V*params->eta + counts->n_dijv[k]);
			if (state > 0 || iter == params->burn_in) {
				latent->beta[k][v] = (state*latent->beta[k][v] + beta[v])/(state+1);
				fprintf(fbeta, "%lf ", latent->beta[k][v]);
			}
			else if (iter > params->ITER) fprintf(fbeta, "%lf ", latent->beta[k][v]);
			else fprintf(fbeta, "%lf ", beta[v]);
		}
		//fprintf(fbeta,"\n");
//		fclose(fphi);
		fclose(fbeta);
		free(beta);
	}

}

double compute_perplexity(char* odir, sctm_data* data, double *docLogLikelihood, int *count){
	
	printf("compute perplexity");
	int d = 0, i=0;
	double totalDocLogLikelihood = 0;
	double perplexity = 0;
	double totalWords = 0;

	for(d=0; d<data->D; d++){
		docLogLikelihood[d] = docLogLikelihood[d]-log(*count);
		printf("\ndocLogLikelihood %.3f", docLogLikelihood[d]);
		fflush(stdout);
		totalDocLogLikelihood += docLogLikelihood[d];
		documents *doc = &(data->docs[d]);
		totalWords += doc->N;
		for(i=0; i<doc->C; i++){
			comment *cmnt = &(doc->cmnts[i]);
			totalWords += cmnt->N;
		}
	}

	perplexity = totalDocLogLikelihood/totalWords;
	printf("\nperplexity---- %.3f %.3f", perplexity, totalWords);
	fflush(stdout);
	perplexity = exp(-perplexity);
	printf("\nperplexity %.3f", perplexity);
	fflush(stdout);
	return perplexity;
}

double logSum(double log_a, double log_b){
	if(log_a < 1e-20)
		return log_b;
	else if (log_b < 1e-20)
		return log_a;
	else if(log_a < log_b)
		return log_b+log(1 + exp(log_a-log_b));
	else 
		return log_b+log(1 + exp(log_a-log_b));

}

double compute_likelihood(sctm_data* data, sctm_params * params, sctm_latent* latent, sctm_counts* counts, double *docLogLikelihood, int *token){
	int d=0, i=0, c=0, n, k, j, v, t;
	double likelihood;
	double tempLike; 
	double totalLikelihood;
	int wordNum = 0;
	double beta, like, theta, eps;

	eps = 1e-10;
	
	// for(d=0; d<data->D; d++){
	// 	printf("docLogLikelihood   %.3f\n", *(docLogLikelihood+d));
	// }

	for(d=0; d<data->D; d++){
		documents *doc = &(data->docs[d]);
		likelihood = 0;
		for(i=0; i<doc->S; i++){
			sentence *sent = &(doc->sents[i]);
			for(n=0; n<sent->N; n++){
				wordNum += 1;
				like = 0;
				for (k=0; k < params->K; k++) {
					v = sent->words[n];
					j = latent->b[d][i][n];

					if (params->trte == 1) beta = latent->beta[k][v];
					else beta = counts->n_dij[k][v] *1.0 / counts->n_dijv[k];

					if(counts->n_iv[d][j][k] < 0 || counts->n_ikv[d][j] <= 0) debug("in lik djk 0");

					if (counts->n_ikv[d][j] == 0) continue;
					theta = counts->n_iv[d][j][k]*1.0 / counts->n_ikv[d][j];
					like += beta*theta;
				}
				if (like < eps) {
					printf("like:%lf\n",like);
					//debug("like too small");
				}
				likelihood += log(like + eps);
			}
			
		}


		printf("parent like:%lf\n",likelihood);
		for(i=0; i<doc->C; i++){
			
			comment* cmnt = &(doc->cmnts[i]);
			for(n=0; n<cmnt->N; n++){
				wordNum += 1;
				like = 0;
				
				for(k=0; k<params->K+1; k++){
					// k = latent->y[d][i][n];
					tempLike = 0;
					v = cmnt->words[n];
					t = latent->t[d][i][n];
					
					// if (t == 0) k = params->K;
					if (params->trte == 1) beta = latent->beta[k][v];
					else beta = counts->n_dij[k][v] *1.0 / counts->n_dijv[k];

					if(t==1 && (counts->m[d][i][k] < 0 || counts->m_k[d][i] <= 0)) debug("in lik djk 0");

					if(k==params->K){
						theta = (cmnt->N-counts->m_1k[d][i])*1.0/cmnt->N;
					}	// theta = (cmnt->N-counts->m_1k[d][i]+eps)*1.0/(cmnt->N+eps*(params->K+1));
					else{
						theta = (counts->m_1[d][i][k])*1.0/cmnt->N;
						// theta = (counts->m_1[d][i][k]+eps)*1.0/(cmnt->N+eps*(params->K+1));
					}

					// tempLike = log(theta)+log(beta);
					// printf("tempLike:%lf\n",tempLike);
					// if(like=0)
					// 	like = tempLike;
					// else 
					// 	like = logSum(like, tempLike);
					like += beta*theta;
					// printf("tempLike:%lf \t theta:%lf \t beta:%lf\n", like, theta, beta);

				}
				if (like < eps) {
					printf("like:%lf\n",like);
					debug("like too small");
				}
				// likelihood += like;
				likelihood += log(like + eps);
			}
		}

		// totalLikelihood += likelihood;
		// printf("child like:%lf\n",likelihood);

		if(!((docLogLikelihood[d]>0)||(docLogLikelihood[d]<0))){
			// printf("doc LogLikelihood == 0");
			*(docLogLikelihood+d) = likelihood;
		// *(docLogLikelihood+d) += likelihood;
		}
		else{
			// printf("doc LogLikelihood > 0");
			*(docLogLikelihood+d) = logSum(*(docLogLikelihood+d), likelihood);
			// printf("end value doc LogLikelihood%.lf\n", exp(*(docLogLikelihood+d)));
		}
	}

	// *docLogLikelihood = totalLikelihood;
	*token = wordNum;
	// for(d=0; d<data->D; d++){
	// 	printf("end value doc LogLikelihood%.3f\n", *(docLogLikelihood+d));
	// }

	return 0;
}




// double compute_perplexity(char* odir, sctm_data* cdata, sctm_params* params, sctm_latent* latent, sctm_counts* counts) {
// 	double likelihood_art, likelihood_cmnt, perp;
// 	int DxN_art, DxN_cmnt, totalWords;
// 	double totalLikelihood = 0;
// 	// compute_likelihood(cdata, params, latent, counts, &likelihood_art, &DxN_art);
// 	// compute_likelihood_cmnt(cdata, params, latent, counts, &likelihood_cmnt, &DxN_cmnt);
// // totalWords = DxN_art+DxN_cmnt;
// 	// totalLikelihood = likelihood_art+likelihood_cmnt;
// 	//double perp = exp(-(likelihood));
// 	compute_likelihood(cdata, params, latent, counts, &totalLikelihood, &totalWords);
//  // compute_likelihood(sctm_data* data, sctm_params * params, sctm_latent* latent, sctm_counts* counts, double *docLogLikelihood, int *token){

	
// 	char fname[500];
// 	FILE *fl;
// 	sprintf(fname, "%s/perplexity.txt", odir);
// 	fl = fopen(fname, "a");
	
// 	perp = exp(-(totalLikelihood)/(totalWords));

// 	printf("perplexity %.3f %d\n", perp, totalWords);
// 	fprintf(fl, "perplexity: %.2f\n", perp);
// 	fclose(fl);

// 	return(perp);
// }

// double compute_perplexity(char* odir, sctm_data* cdata, sctm_params* params, sctm_latent* latent, sctm_counts* counts, double *docLogLikelihood) {
// 	double likelihood_art, likelihood_cmnt, perp;
// 	int DxN_art, DxN_cmnt, totalWords, d;
// 	double totalLikelihood = 0;
// 	// compute_likelihood(cdata, params, latent, counts, &likelihood_art, &DxN_art);
// 	// compute_likelihood_cmnt(cdata, params, latent, counts, &likelihood_cmnt, &DxN_cmnt);
// // totalWords = DxN_art+DxN_cmnt;
// 	// totalLikelihood = likelihood_art+likelihood_cmnt;
// 	//double perp = exp(-(likelihood));
// 	compute_likelihood(cdata, params, latent, counts, docLogLikelihood, &totalWords);
//  // compute_likelihood(sctm_data* data, sctm_params * params, sctm_latent* latent, sctm_counts* counts, double *docLogLikelihood, int *token){

// 	for(d=0; d<cdata->D; d++)
// 		totalLikelihood += docLogLikelihood[d];
	
// 	char fname[500];
// 	FILE *fl;
// 	sprintf(fname, "%s/perplexity.txt", odir);
// 	fl = fopen(fname, "a");
	
// 	perp = exp(-(totalLikelihood)/(totalWords));

// 	printf("perplexity %.3f %d\n", perp, totalWords);
// 	fprintf(fl, "perplexity: %.2f\n", perp);
// 	fclose(fl);

// 	return(perp);
// }

// double compute_likelihood(sctm_data* cdata, sctm_params* params, sctm_latent* latent, sctm_counts* counts, double* result, int* tokens) {
// 	int d, i, j, k, v, n;//, t, l;
// 	double likelihood = 0;
// 	double beta,like, theta, eps;

// 	int DxN = 0;

// 	eps = 1e-10;
// 	for (d=0; d<cdata->D; d++) {
// 		documents* doc = &(cdata->docs[d]);

// 		for (i=0; i < doc->S; i++) {
// 			sentence* sent = &(doc->sents[i]);
// 			for (n=0; n < sent->N; n++){
// 				like = 0;
// 				for (k=0; k < params->K; k++) {
// 					v = sent->words[n];
// 					j = latent->b[d][i][n];

// 					if (params->trte == 1) beta = latent->beta[k][v];
// 					else beta = counts->n_dij[k][v] *1.0 / counts->n_dijv[k];

// 					if(counts->n_iv[d][j][k] < 0 || counts->n_ikv[d][j] <= 0) debug("in lik djk 0");

// 					if (counts->n_ikv[d][j] == 0) continue;
// 					theta = counts->n_iv[d][j][k]*1.0 / counts->n_ikv[d][j];
// 					like += beta*theta;
// 				}
// 				if (like < eps) {
// 					printf("like:%lf\n",like);
// 					//debug("like too small");
// 				}
// 				likelihood += log(like + eps);
// 				DxN += 1;
// 			}
// 		}
// 	}
	
// 	*result = likelihood;
// 	*tokens = DxN;

// 	return(likelihood/DxN);
// }

// double compute_likelihood_cmnt(sctm_data* cdata, sctm_params* params, sctm_latent* latent, sctm_counts* counts, double* result, int* tokens) {
// 	int d, i, t, k, v, n;//, t, l;
// 	double likelihood = 0;
// 	double beta,like, eps;

// 	int DxN = 0;

// 	eps = 1e-10;
// 	for (d=0; d<cdata->D; d++) {
// 		documents* doc = &(cdata->docs[d]);
// 		for (i=0; i < doc->C; i++) {
// 			comment* cmnt = &(doc->cmnts[i]);
// 			for (n=0; n < cmnt->N; n++) {
// 				like = 0;
// 				k = latent->y[d][i][n];
// 				v = cmnt->words[n];
// 				t = latent->t[d][i][n];
// 				if (t == 0) k = params->K;
// 				if (params->trte == 1) beta = latent->beta[k][v];
// 				else beta = counts->n_dij[k][v] *1.0 / counts->n_dijv[k];

// 				if(t==1 && (counts->m[d][i][k] < 0 || counts->m_k[d][i] <= 0)) debug("in lik djk 0");

// 				like += beta;

// 				//}
// 				if (like < eps) {
// 					printf("like:%lf\n",like);
// 					debug("like too small");
// 				}
// 				likelihood += log(like + eps);
// 				DxN += 1;
// 			}
// 		}
// 	}
// 	*result = likelihood;
// 	*tokens = DxN;

// 	return(likelihood/DxN);
// }

