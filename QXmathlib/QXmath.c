#include <stdio.h>
#include <math.h>
#include "QXmath.h"

const double Cpd = 0.2403;
const double Cpv = 0.445;
const double Rd = 6.85578 * 0.01;
const double Rw = 11.017874 * 0.01;
const double E0 = 6.1078;
const double T0 = 273.16;
const double C1 = 0.57;
const double L0 = 597.4;
const double Cf = 0.06;
const double Lf = 79.72;
const double A = 2.38844 * 1e-8;
const double g = 980.665;
const double rd = 9.76;
const double PI = 3.14159265;
const double R = 6.371229 * 1e3;
const double C = 1.002;

double Tc(double p, double t, double td){
    double Etd, Z, Z0, step = 10.0, w, m1;
    double out, T, Td;
    T = T0 + t;
    Td = T0 + td;
    Etd = E_WATER(td);
    w = (Rd / Rw) * Etd / (p - Etd);
    m1 = (Cpd * (1 + Cpv * w / Cpd)) / (Rd * (1 + w / (Rd / Rw)));
    Z0 = pow(T, m1) / Etd;
    out = Td;
    Z = pow(out, m1) / E_WATER(out - T0);
    while(fabs(Z - Z0) > 10){
        if(Z < Z0){
            out = out - step;
        }
        else{
            out = out +step;
            step = step / 5;
            out = out -step;
        }
        Z = pow(out, m1) / E_WATER(out - T0);
        if (step < 1e-7) break;
    }
    return out - T0;
}

double E_WATER(double td){
    double T, E;
    T = T0 + td;
    E = ((L0 + C1 * T0) * (T - T0)) / (Rw * T0 * T);
    E = exp(E);
    E = E * E0 * pow(T0 / T, C1 / Rw);
    return E;
}

double E_ICE(double td){
    double Ls, T, E;
    if(td > 0.0){
        printf("td>0, can not calculate E on ice\n");
        return NAN;
    }
    T = T0 + td;
    Ls = L0 + Lf - Cf * (T - T0);
    E = ((Ls + Cf * T0) * (T - T0)) / (Rw * T0 * T);
    E = exp(E);
    E = E * E0 * pow(T0 / T, Cf / Rw);
    return E;
}

double Qse(double P, double t, double td){
    double tc, Lc, T, E, w, qse, Kd;
    Kd = Rd / Cpd;
    T = t + T0;
    tc = Tc(P, t, td) + T0;
    Lc = L0 - C1 * (tc - T0);
    E = E_WATER(td);
    w = 0.622 * E / (P - E);
    qse = T * pow((1000 / (P - E)), Kd) * exp((Lc * w) / (Cpd * tc));
    return qse - T0;
}

double Qp(double p, double t){
    double P, T, Kd, out;
    Kd = Rd / Cpd;
    P = p;
    T = t + T0;
    out = T * pow((1000 / P), Kd);
    return out;
}

double Ttdm(double p, double t, double td, double z){
    double E, Tz, Tl, Tt, q, L;
    E = E_WATER(td);
    Tz = 100 * A * g * z / Cpd;
    L = L0 - C1 * td;
    q = 0.622 * E / (p - 0.378 * E);
    Tl = L * q / Cpd;
    Tt = t + Tz + Tl;
    return Tt;
}

double Ttgk(double p, double t, double td, double z, double v){
    double L, q, out;
    L = L0 - C1 * td;
    q = qgk(p, td) / 1000;
    out = t + 1000 * A * g * z / Cpd + L * q / Cpd + A * v * v / (2 * Cpd);
    return out;
}

double Etotd(double E){
    double td = -60;
    double e, step = 20;
    e = E_WATER(td);
    do{
        if(e > E){
            step = step / 2;
            td = td - step;
        }
        else{
            td = td + step;
        }
        e = E_WATER(td);
    }while(fabs(e-E) > 0.0001);
    return td;
}

double qgk(double p, double td){
    double E;
    double out;
    E = E_WATER(td);
    out = 622 * E / (p - 0.378 * E);
    return out;
}

double rgk(double t, double td){
    double out;
    double Et, Etd;
    Et = E_WATER(t);
    Etd = E_WATER(td);
    out = (Etd / Et) * 100;
    return out;
}

double Fc(double p, double td){
    double L, qs, Etd, Td, Fcd, out;
    Td = T0 + Td;
    Etd = E_WATER(td);
    qs = 1000 * Rd * Etd / (p - 0.378 * Etd) / Rw;
    L = L0 - C1 * (Td - T0);
    Fcd = (qs / (p - 0.378 * Etd)) * (Rd * L/ (Cpd * Rw * Td) - 1);
    out = 100 * Fcd;
    out = out / (1 + (L * L * qs * 0.001 / (Cpd * Rw * Td * Td)) * (p / (p - 0.378 * Etd)));
    return out;
}

double sqtl(double p, double td, double v){
    double q, out;
    q = qgk(p, td);
    out = v * q / g;
    return out;
}

double Rm(double P, double t){
    double L, T, E, R, w, q, out;
    T = t + T0;
    L = L0 - C1 * t;
    E = E_WATER(t);
    w = 0.622 * E / (P - E);
    q = 0.622 * E / (P - 0.378 * E);
    R = Rd * (1 + 0.608 * q);
    out = rd * (1 + (L * w / (R * T)));
    out = out / (1 + (0.622 * L * L * w / (Cpd * Rd * T * T)));
    return out;
}

void fsfj(float fd, float ff, float *u, float *v){
    *u = ff * (cos((270 - fd) * PI / 180));
    *v = ff * (sin((270 - fd) * PI / 180));
}

void fsfjb(float fd, float ff, float *u, float *v){
    *u = ff * (cos((225 - fd) * PI / 180));
    *v = ff * (sin((225 - fd) * PI / 180));    
}

void fshc(float U, float V, float *fd, float *ff){
    *ff = hypot(U, V);
    *fd = atan2(U, V);
    *fd = (*fd) * 180 / PI;
    *fd = (*fd) + 180;
}

double scfsd_zfwg(double u_E, double u_W, double v_N, double v_S, double d){
    return (((u_E - u_W) + (v_N - v_S)) / (2 * d * 1000)) * 1e5;
}

double scfwd_zfwg(double v_E, double v_W, double u_N, double u_S, double d){
    return (((v_E - v_W) - (u_N - u_S)) / (2 * d * 1000)) * 1e5;
}

double scfsd_jwwg(double u_E, double u_W, double v_N, double v_S, double wgj, double wd, double v){
    double d, out;
    d = wgj * 111.2;
    out = ((v / (R * 1000)) * (tan(wd * PI / 180))) * 1e5;
    out = (((u_E - u_W) + (v_N - v_S)) / (2 * d * 1000)) * 1e5 - out;
    return out;
}

double scfwd_jwwg(double v_E, double v_W, double u_N, double u_S, double wgj, double wd, double u){
    double d, out;
    d = wgj * 111.2;
    out = ((u / (R * 1000)) * (tan(wd * PI / 180))) * 1e5;
    out = (((v_E - v_W) + (u_N - u_S)) / (2 * d * 1000)) * 1e5 - out;
    return out;
}

double scfsd_yxy(jwdf_type dat[3]){
    float x1, x2, x3, y1, y2, y3;
    float u1, u2, u3, v1, v2, v3;
    double out;
    x1 = 0; y1 = 0;
    x2 = (dat[1].jd - dat[0].jd) * PI * R * cos((dat[1].wd + dat[0].wd) / 2 * PI / 180) / 180;
    y2 = (dat[1].wd - dat[0].wd) * PI * R / 180;
    x3 = (dat[2].jd - dat[0].jd) * PI * R * cos(((dat[2].wd + dat[0].wd) / 2) * PI /180) / 180;
    y3 = (dat[2].wd - dat[0].wd) * PI * R / 180;
    fsfj(dat[0].fd, dat[0].ff, &u1, &v1);
    fsfj(dat[1].fd, dat[1].ff, &u2, &v2);
    fsfj(dat[2].fd, dat[2].ff, &u3, &v3);
    out = ((u1 - u2) * (y2 - y3) - (u2 - u3) * (y1 - y2) + (x1 - x2) * (v2 - v3) - (x2 - x3) * (v1 - v2)) / ((x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2));
    return out * 1e5;
}

double scfwd_yxy(jwdf_type dat[3]){
    return;
}

double showalter(double t8, double td8, double t5){
    double w, m1, m2, Qc, Q5, P5 = 500;
    double Etd8, Eout, ETa, T8, P8 = 850;
    double Ta, Pa;
    double out, step = 10;
    T8 = T0 + t8;
    Etd8 = E_WATER(t8);
    w = (Rd / Rw) * Etd8 / (P8 - Etd8);
    m1 = (Cpd * (1 + Cpv * w / Cpd)) / (Rd * (1 + w / (Rd / Rw)));
    Ta = Tc(P8, t8, td8) + T0;
    Pa = P8 * (pow(Ta / T8, m1));
    m2 = (Cpd / Rd) * (1 + C * w / Cpd);
    ETa = E_WATER(Ta - T0);
    Qc = log((Pa - ETa) / pow(Ta, m2)) - (0.622 / Rd) * (((L0 + C1 * (T0 - Ta)) / Ta) * (ETa / (Pa - ETa)));
    out = T8;
    Eout = E_WATER(out - T0);
    Q5 = log((P5 - Eout) / pow(out, m2)) - (0.622 / Rd) * (((L0 + C1 * (T0 - out)) / out) * (Eout / (P5 - Eout)));
    while(fabs(Qc - Q5) > 0.0001){
        if(Qc > Q5){
            out = out - step;
        }
        else{
            out = out + step;
            step = step / 5;
            out = out - step;
        }
        Eout = E_WATER(out - T0);
        Q5 = log((P5 - Eout) / pow(out, m2)) - (0.622 / Rd) * (((L0 + C1 * (T0 - out)) / out) * (Eout / (P5 - Eout)));
    }
    return t5 - (out - T0);
}

double richardson(double pdn, double tdn, double fddn, double ffdn, double pup, double tup, double fdup, double ffup){
    float udn, vdn, uup, vup;
    double out;
    double pjp, pjT, cha_p, cha_t, cha_u, cha_v;
    fsfj(fdup, ffup, &uup, &vup);
    fsfj(fddn, ffdn, &udn, &vdn);
    pjp = (pup + pdn) / 2;
    pjT = T0 + (tup + tdn) / 2;
    cha_p = pup - pdn;
    cha_t = tup - tdn;
    cha_u = uup - udn;
    cha_v = vup - vdn;
    out = -1e-4 * (Rd * cha_p / pjp) * (cha_t - (A * Rd * pjT / Cpd) * (cha_p / pjp));
    out = out / (cha_u * cha_u + cha_v * cha_v);
    return out;
}