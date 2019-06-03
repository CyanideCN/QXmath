struct jwdf_type {
    float jd;
    float wd;
    float fd;
    float ff;
};

typedef struct jwdf_type jwdf_type;

double Tc(double p, double t, double td);
double E_WATER(double td);
double E_ICE(double td);
double Qse(double P, double t, double td);
double Qp(double P, double t);
double Ttdm(double p, double t, double td, double z);
double Ttgk(double p,  double t, double td, double z, double v);
double Etotd(double E);
double qgk(double p, double td);
double rgk(double p, double td);
double Fc(double p, double td);
double Rm(double p, double t);
double sqtl(double p, double td, double v);
void fsfj(float fd, float ff, float *u, float *v);
void fsfjb(float fd, float ff, float *u, float *v);
void fshc(float U, float V, float *fd, float *ff);
double scfsd_zfwg(double u_E, double u_W, double v_N, double v_S, double d);
double scfwd_zfwg(double v_E, double v_W, double u_N, double u_S, double d);
double scfsd_jwwg(double u_E, double u_W, double v_N, double v_S, double wgj,  double wd, double v);
double scfwd_jwwg(double v_E, double v_W, double u_N, double u_S, double wgj, double wd, double u);
double scfsd_yxy(jwdf_type dat[3]);
double scfwd_yxy(jwdf_type dat[3]);
double dzfug_zfwg(double N_H, double S_H, double d, double wd);
double dzfvg_zfwg(double E_H, double W_H, double d, double wd);
double dzfug_jwwg(double N_H, double S_H, double wgj, double wd);
double dzfvg_jwwg(double E_H, double W_H, double wgj, double wd);
double dzfwd_zfwg(double H1, double H2, double H3, double H4, double H0, double d, double wd);
double czsdzl(double Dk_1, double Dk, double P_cha);
double showalter(double t8, double td8, double t5);
double richardson(double pdn, double tdn, double fddn, double ffdn, double pup, double tup, double fdup, double ffup);