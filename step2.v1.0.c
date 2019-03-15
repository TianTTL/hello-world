#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

const long snpEffectNum = 1e5, snpTotalNum = 2e5;
const long phenoNum = 100;
const long simulNum = 50;

double ocorm[1000][1000], ocorm_2[1000][1000], ocorm_3[1000][1000], ocorm_4[1000][1000]; //phenoNum * phenoNum
int qvLen;
long qv[10000]; //resolution, 10000
double E[1000]; //phenoNum
long randomSeed[10000]; //simulNum
double distMt[100000][1000]; //snpEffectNum * phenoNum
double prMt[100000][1000]; //snpEffectNum * phenoNum
double tdff[1000]; //phenoNum
double resMT[100000][1000]; //snpEffectNum * phenoNum
double fdffVect[1000]; //phenoNum

void preCorrection();
void getRandomSeed();
void smltDistMt(long randomSeedCurrent);
void corterm3all(long snpEffectNumCurrent);

void main(){
    long i, j, k;
    int row, col;
    FILE * fout = NULL;
    FILE * fin = NULL;

    time_t startTime,stopTime, smltStartTime, smltEndTime;
    startTime = time(NULL);

    fin = fopen("/liufanGroup/gaoxj/HPC419/2_C/inputFile/100StatsCor", "r");
    for (i = 0; i < phenoNum; i++){
        for (j = 0; j < phenoNum; j++){
            fscanf(fin, "%lf", &ocorm[i][j]);
        }
        fscanf(fin,"\n");
    }
    fclose(fin);
    fin = NULL;

    preCorrection();

    for (i = 0; i < phenoNum; i++){
        E[i] = (double)i + i + 2;
    }
        
    getRandomSeed(simulNum, randomSeed);
    for (i = 0; i < simulNum; i++){
        smltStartTime = time(NULL);
        smltDistMt(randomSeed[i]);

        for (j = 0; j < snpEffectNum; j++)
            for (k = 0; k < phenoNum; k++){
                prMt[j][k] = gsl_cdf_gaussian_P(-fabs(distMt[j][k]), 1.0) * 2; // pnorm
            }

        for (j = 0; j < phenoNum; j++){
            tdff[j] = 0;
        }

        for (j = 0; j < phenoNum; j++)
            for (k = 0; k < phenoNum; k++){
                ocorm_2[j][k] = pow(ocorm[j][k], 2);
                ocorm_3[j][k] = pow(ocorm[j][k], 3);
                ocorm_4[j][k] = pow(ocorm[j][k], 4);
            }

        for (j = 0; j < snpEffectNum; j++){
            corterm3all(j);
        }

        // calc & write fminm
        for (j = 0; j < phenoNum; j++){
            // prepare output file for each pheno
            char fminmOutputPath[255];
            sprintf(fminmOutputPath, "/liufanGroup/gaoxj/HPC419/2_C/outputFile/qutm_%d", j + 1);
            fout = fopen(fminmOutputPath, "a");
            
            double resTransVect[snpEffectNum];
            for (k = 0; k < snpEffectNum; k++){
                resTransVect[k] = log10(resMT[k][j]);  //extract 1 col, transform
            }
            gsl_sort(resTransVect, 1, snpEffectNum);
            for (k = 0; k < qvLen; k++){
                double quantileRlt = gsl_stats_quantile_from_sorted_data(resTransVect, 1, snpEffectNum, 1.0 * (qv[k] - 1) / (snpTotalNum - 1)); //quantile
                double outputRlt = pow(10.0, quantileRlt);
                fprintf(fout, "%15g ", outputRlt);
            }
            fprintf(fout, "\n");
            fclose(fout);
            fout = NULL;
        }

        // calc & write fdffVect
        fout = fopen("/liufanGroup/gaoxj/HPC419/2_C/outputFile/dfm", "a");
        for (j = 0; j < phenoNum; j++){
            fdffVect[j] = tdff[j] / snpEffectNum;
            fprintf(fout, "%15g ", fdffVect[j]);
        }
        fprintf(fout, "\n");
        fclose(fout);
        fout = NULL;

        // calc & write fgpcm
        char fgpcmOutputPath[255];
        sprintf(fgpcmOutputPath, "/liufanGroup/gaoxj/HPC419/2_C/outputFile/qutm_%d", phenoNum + 1);    
        fout = fopen(fgpcmOutputPath, "a");
        
        double unifp[snpEffectNum];
        gsl_rng * r;
        r = gsl_rng_alloc(gsl_rng_default); 
        gsl_rng_set(r, randomSeed[i]);
        for (j = 0; j < snpEffectNum; j++){
            double tmp = gsl_ran_chisq(r, fdffVect[phenoNum - 1]);
            double tmp1 = 1 - gsl_cdf_chisq_P(tmp, fdffVect[phenoNum - 1]);
            unifp[j] = log10(tmp1);
        }
        gsl_rng_free(r);
        gsl_sort(unifp, 1, snpEffectNum);

        for (j = 0; j < qvLen; j++){
            double quantileRlt = gsl_stats_quantile_from_sorted_data(unifp, 1, snpEffectNum, 1.0 * (qv[j] - 1) / (snpTotalNum - 1)); //quantile
            double outputRlt = pow(10.0, quantileRlt);
            fprintf(fout, "%15g ", outputRlt);
        }

        fprintf(fout, "\n");
        fclose(fout);
        fout = NULL;

        smltEndTime = time(NULL);
        printf("%ldth simulation use time: %ld secs\n",i, (smltEndTime - smltStartTime));
    }

    stopTime = time(NULL);
    printf("Use Time: %ld secs\n",(stopTime - startTime));

    return;
}

void preCorrection(){
    int i, j;
    
    qvLen = 0;
    i = 1;
    while (i <= (int)(ceil(snpTotalNum * 0.001))){
        qv[qvLen++] = i++;
    }

    i = (int)(ceil(snpTotalNum * 0.001)) + 1;
    j = (int)(round(snpTotalNum * 0.00001));
    while (i <= (int)(ceil(snpTotalNum * 0.01))){
        qv[qvLen++] = i;
        i += j;   
    }

    i = (int)(ceil(snpTotalNum * 0.01)) + 1;
    j = (int)(round(snpTotalNum * 0.0002));
    while (i <= (int)(ceil(snpTotalNum * 0.999)) - 1){
        qv[qvLen++] = i;
        i += j;   
    }

    i = (int)(ceil(snpTotalNum * 0.999));
    while (i <= snpTotalNum){
        qv[qvLen++] = i++;
    }

    return;
}

void getRandomSeed(){
    long i;
    time_t timeStamp;
    timeStamp = time(NULL);
    gsl_rng * r;
    r = gsl_rng_alloc(gsl_rng_default); 
    gsl_rng_set(r, time(&timeStamp));
    for (i = 0; i < simulNum; i++){
        randomSeed[i] = gsl_rng_uniform_int(r, 2147483648);
    }
    gsl_rng_free(r);

    return;
}

void smltDistMt(long randomSeedCurrent){
    long i, j;

    //cholesky
    gsl_matrix * A = gsl_matrix_alloc(phenoNum, phenoNum);
    for (i = 0; i < phenoNum; i++)
        for (j = 0; j < phenoNum; j++){
            gsl_matrix_set(A, i, j, ocorm[i][j]);
        }
    gsl_linalg_cholesky_decomp1(A);

    //random number generator
    gsl_rng * r;
    r = gsl_rng_alloc(gsl_rng_default); 
    gsl_rng_set(r, randomSeedCurrent);
    for (i = 0; i < snpEffectNum; i++){
        gsl_vector * mu = gsl_vector_calloc(phenoNum);
        gsl_vector * generatorRlt = gsl_vector_alloc(phenoNum);
        gsl_ran_multivariate_gaussian(r, mu, A, generatorRlt);
        for (j = 0; j < phenoNum; j++){
            distMt[i][j] = gsl_vector_get(generatorRlt, j);
        }
    }
    gsl_rng_free(r);

    return;
}

void corterm3all(long snpEffectNumCurrent){
    long i, j;

    long stepw[phenoNum]; //phenoNum    
    double distVect_order[phenoNum]; //phenoNum
    double prVect_order[phenoNum]; //phenoNum
    double ocorm_order[phenoNum][phenoNum], ocorm_order_2[phenoNum][phenoNum], ocorm_order_3[phenoNum][phenoNum], ocorm_order_4[phenoNum][phenoNum]; //phenoNum * phenoNum
    double corm[phenoNum][phenoNum]; //phenoNum * phenoNum
    double statVec[phenoNum]; //phenoNum
    double curcov[phenoNum]; //phenoNum
    double V[phenoNum], scaf[phenoNum], dff[phenoNum]; //phenoNum

    // heap sort index
    gsl_vector * B = gsl_vector_alloc(phenoNum);
    for (i = 0; i < phenoNum; i++)
        gsl_vector_set(B, i, prMt[snpEffectNumCurrent][i]);
    gsl_permutation * perm = gsl_permutation_alloc(phenoNum);
    gsl_sort_vector_index(perm, B);
    for (i = 0; i < phenoNum; i++){
        stepw[i] = gsl_permutation_get(perm, i);
    }

    for (i = 0; i < phenoNum; i++){
        distVect_order[i] = distMt[snpEffectNumCurrent][stepw[i]];
        prVect_order[i] = prMt[snpEffectNumCurrent][stepw[i]];
        for (j = 0; j < phenoNum; j++){
            ocorm_order[i][j] = ocorm[stepw[i]][stepw[j]];
            ocorm_order_2[i][j] = ocorm_2[stepw[i]][stepw[j]];
            ocorm_order_3[i][j] = ocorm_3[stepw[i]][stepw[j]];
            ocorm_order_4[i][j] = ocorm_4[stepw[i]][stepw[j]];
        }
    }

    for (i = 0; i < phenoNum; i++)
        if (distVect_order[i] < 0)
            for (j = 0; j < phenoNum; j++)
            {
                ocorm_order[i][j] = -ocorm_order[i][j];
                ocorm_order[j][i] = -ocorm_order[j][i];
                ocorm_order_3[i][j] = -ocorm_order_3[i][j];
                ocorm_order_3[j][i] = -ocorm_order_3[j][i];
            }

    for (i = 0; i < phenoNum; i++){
        for (j = 0; j < phenoNum; j++){
            corm[i][j] = 3.2630398097 * ocorm_order[i][j] + 0.7095678755 * ocorm_order_2[i][j] + 0.0268257772 * ocorm_order_3[i][j] + 0.0005732151 * ocorm_order_4[i][j];
        }
    }

    // calc statVec
    for (i = 0; i < phenoNum; i++){
        statVec[i] = -2 * log(prVect_order[i]);
    }
    for (i = 1; i < phenoNum; i++){
        statVec[i] = statVec[i] + statVec[i - 1];
    }

    // calc curcov
    curcov[0] = 0;
    for (i = 1; i < phenoNum; i++){
        curcov[i] = curcov[i - 1];
        for (j = 0; j < i; j++){
            curcov[i] += corm[j][i];
        }
    }

    // tdff
    for (i = 0; i < phenoNum; i++){
        V[i] = curcov[i] + curcov[i] + E[i] + E[i];
        scaf[i] = V[i] / (E[i] + E[i]);
        dff[i] = 2 * E[i] * E[i] / V[i];
        tdff[i] += dff[i];
    }

    // calc resMT
    for (i = 0; i < phenoNum; i++){
        resMT[snpEffectNumCurrent][i] = 1 - gsl_cdf_chisq_P(statVec[i]/ scaf[i], dff[i]);
    }
    return;
}