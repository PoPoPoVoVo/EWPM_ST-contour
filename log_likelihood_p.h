#include <TMath.h>
struct exp_data{
    std::string name;
    double O_exp;
    double O_sm;
    double sigma;
    double a;
    double b;
};

double log_likelihood_p(double S, double T) {
    // Observable 정보: O_exp, O_sm, uncertainty, a, b 
    exp_data obs[] = {
        //1st try
        {"W/Z", 0.8791, 0.8787, 0.0034, -3.15e-3, 4.86e-3 },
        // {"W", 80.13, 80.09, 0.31, (-3.15e-3)*(91.15), (4.86e-3)*(91.15) },
        // {"W/Z", 80.435, 80.126, 0.195, (-3.15e-3)*(91.174), (4.86e-3)*91.174 },
        {"Gamma_Z",  2.487,  2.484,  0.009, -9.58e-3 , 2.615e-2 },
        // {"Gamma_leptonic", 0.0835,   0.0835,   0.0004,    -1.91e-4,  7.83e-4},
        // {"Gamma_had", 1.7440, 1.7348, 0.01, -9.0e-3, 1.993e-2},
        {"R_Z", 20.94, 20.78, 0.12, -5.99e-2, 4.24e-2},
        {"s_*_at_mz", 0.2317, 0.2337, 0.0030,  3.59e-3 ,  -2.54e-3},
        {"FB_asym_b", 0.135, 0.0848, 0.031,  -1.97e-2 , 1.40e-2  },
        // {"FB_asym_l", 0.016, 0.0126, 0.005, -6.72e-3, 4.76e-3},
        {"tau_polarization", -0.152, -0.1297, 0.045, 2.82e-2 , -2.00e-2 },
        {"g_L", 0.2977, 0.3001, 0.0042, -2.67e-3 , 6.53e-3},
        {"g_R", 0.0317, 0.0302, 0.0034, 9.17e-4, -1.94e-4 },
        // {"R_nu", 0.3098, 0.3126, 0.0433, -2.34e-2, 6.46e-2},
        // {"R_nubar", 0.3831, 0.3824, 0.1336, -2.77e-3, 6.03e-3},
        {"Cs_weak_charge", -71.04, -73.31, 1.808, -0.790, -0.011}
        // {"Cs_weak_charge", -69.4, -73.31, 1.5, -0.790, -0.011}
    };

    int N_data = sizeof(obs)/sizeof(exp_data); // N_data: number of data

    double logL = 0;
    for (int i = 0; i < N_data; ++i) {
        double theo = obs[i].O_sm + obs[i].a *S + obs[i].b *T;
        double res = (obs[i].O_exp - theo) / obs[i].sigma;
        logL += -0.5 * TMath::Power(res,2);
    }

    return logL;
}