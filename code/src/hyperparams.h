#ifndef _HYPERPARAMS_H
#define _HYPERPARAMS_H

extern int WRITE_STATUS;

#define TIMESTAMP 0

#define LOCAL_SOURCE "local"

#define MISSING_VALUE 999999
#define TIME_p_STS_1 0.8902353
#define TIME_k_STS_1 737.7606011
#define TIME_shift_STS_1 44
#define MAX_time_diff_STS_1 90
#define TIME_p_STA_1 0.1279508
#define TIME_k_STA_1 18.3327414
#define TIME_shift_STA_1 44
#define MAX_time_diff_STA_1 329
#define TIME_p_ATS_1 0.4264933
#define TIME_k_ATS_1 285.5214802
#define TIME_shift_ATS_1 365
#define MAX_time_diff_ATS_1 80
#define TIME_p_AUS_1 0.1798645
#define TIME_k_AUS_1 90.2466005
#define TIME_shift_AUS_1 365
#define MAX_time_diff_AUS_1 262
#define TIME_p_ATA_1 0.1831738
#define TIME_k_ATA_1 93.1413123
#define TIME_shift_ATA_1 365
#define MAX_time_diff_ATA_1 313
#define TIME_p_AUA_1 0.1185086
#define TIME_k_AUA_1 60.6120787
#define TIME_shift_AUA_1 365
#define MAX_time_diff_AUA_1 399
// less than loglikelihood of SI = 0
#define MIN_loglik_time -86.0

#define TIME_p_STS_2 0.1213433
#define TIME_k_STS_2 24.0314304
#define TIME_shift_STS_2 44
#define MAX_time_diff_STS_2 720
#define TIME_p_STS_3 0.06264466
#define TIME_k_STS_3 18.32845876
#define TIME_shift_STS_3 44
#define MAX_time_diff_STS_3 900
#define TIME_p_STA_2 0.06408474
#define TIME_k_STA_2 15.28356907
#define TIME_shift_STA_2 44
#define MAX_time_diff_STA_2 782
#define TIME_p_STA_3 0.04935896
#define TIME_k_STA_3 16.76983464
#define TIME_shift_STA_3 44
#define MAX_time_diff_STA_3 970
#define TIME_p_ATS_2 0.1570167
#define TIME_k_ATS_2 86.4214609
#define TIME_shift_ATS_2 365
#define MAX_time_diff_ATS_2 690
#define TIME_p_ATS_3 0.09451821
#define TIME_k_ATS_3 58.20095270
#define TIME_shift_ATS_3 365
#define MAX_time_diff_ATS_3 870
#define TIME_p_AUS_2 0.08840967
#define TIME_k_AUS_2 49.39713807
#define TIME_shift_AUS_2 365
#define MAX_time_diff_AUS_2 821
#define TIME_p_AUS_3 0.06670912
#define TIME_k_AUS_3 43.42367742
#define TIME_shift_AUS_3 365
#define MAX_time_diff_AUS_3 1004
#define TIME_p_ATA_2 0.1039955
#define TIME_k_ATA_2 58.8211428
#define TIME_shift_ATA_2 365
#define MAX_time_diff_ATA_2 753
#define TIME_p_ATA_3 0.0760096
#define TIME_k_ATA_3 49.6186452
#define TIME_shift_ATA_3 365
#define MAX_time_diff_ATA_3 940
#define TIME_p_AUA_2 0.07097978
#define TIME_k_AUA_2 42.40810585
#define TIME_shift_AUA_2 365
#define MAX_time_diff_AUA_2 891
#define TIME_p_AUA_3 0.0585646
#define TIME_k_AUA_3 40.6748579
#define TIME_shift_AUA_3 365
#define MAX_time_diff_AUA_3 1070

extern double LIK_founder_time_S;
extern double LIK_founder_time_A;
extern double LIK_founder_space;

extern vector<double> sigma_AUA_k_1;
extern vector<double> sigma_AUA_k_2;
extern vector<double> sigma_AUA_k_3;
extern vector<double> sigma_AUS_k_1;
extern vector<double> sigma_AUS_k_2;
extern vector<double> sigma_AUS_k_3;
extern vector<double> sigma_STA_k_1;
extern vector<double> sigma_STA_k_2;
extern vector<double> sigma_STA_k_3;
extern vector<double> sigma_ATA_k_1;
extern vector<double> sigma_ATA_k_2;
extern vector<double> sigma_ATA_k_3;
extern vector<double> sigma_STS_k_1;
extern vector<double> sigma_STS_k_2;
extern vector<double> sigma_STS_k_3;
extern vector<double> sigma_ATS_k_1;
extern vector<double> sigma_ATS_k_2;
extern vector<double> sigma_ATS_k_3;

#define MAX_time_GU 864
#define MAX_time_GT 864

#define MAX_time_IDP_A 365
#define MAX_time_IDP_S 44

#define TIME_p_IDP_A 0.023
#define TIME_k_IDP_A 1.72

#define TIME_p_GU 0.064
#define TIME_k_GU 5.78

#define TIME_p_GT 0.58
#define TIME_k_GT 66.10

#define TIME_lambda_IDP_S 16.63
#define TIME_p_IDP_S 0.053
#define TIME_k_IDP_S 3.04

//#define MIN_loglik_space -20.0
//#define MIN_loglik_space -1836.0
// calculated as 2 * dnorm(74.33613, sd = sqrt(1e-3 * 101.6), log = T)  where 74.33613 is the radius of swaziland
#define MIN_loglik_space -54388.0

//#define INTERCEPT_UNKNOWNLOCAL -2.547
//#define SLOPE_UNKNOWNLOCAL -1.000
//#define INTERCEPT_UNKNOWNLOCAL -6.425
#define INTERCEPT_UNKNOWNLOCAL -7.151
#define SLOPE_UNKNOWNLOCAL -1.0

#define INTERCEPT_UNKNOWNCOORDS -2.531
#define SLOPE_UNKNOWNCOORDS -2.0

//#define INTERCEPT_UNKNOWNPARENT 7.492471
//#define SLOPE_UNKNOWNPARENT -3.959452

#endif // _HYPERPARAMS_H
