#ifndef EOS_H_
#define EOS_H_

#include "JetScapeLogger.h"
using namespace Jetscape;

#include <string>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include <cmath>

#ifndef PI
#define PI (3.14159265358979324)
#endif

//#include "PPMutil.h"

static const double N_c = 3.0; // # of color
static const double N_s = 2.0; // # of spin
static const double N_f = 3.0; // # of active flavors
static const double d_q = 2.0*N_s*N_c*N_f; // degree of freedom of quark
static const double d_g = (N_c * N_c -1.0)*N_s; // degree of freedom of gluon
static const double d_qgp = (d_q*7.0/8.0) + d_g;

static const double de1 = 0.00000001;
static const double de2 = 0.0001;
static const int dim1 = 10001;
static const int dim2 = 1000000;
static const double ebound = de1*(dim1-1.0);

class EOS{
    public:
    virtual ~EOS(){};
    virtual double P(double e) { return 0.0; }
    virtual double T(double e) { return 0.0; }
    virtual double SV(double p, double e) { return 0.0; }
    virtual double getTFromS(double s) { return 0.0; }
    virtual double E(double T) { return 0.0; }
    virtual double S(double T) { return 0.0; }
    virtual double getEFromS(double s) { return 0.0; }
    virtual std::string getName(){return "null";}
};

class MasslessIdeal: public EOS{
    public:
    MasslessIdeal(){JSINFO << "<-[PPM] Equation of State: Massless Ideal QGP Gas ->";};
    ~MasslessIdeal(){JSINFO << "<-[PPM] Deleting EoS (Massless Ideal QGP Gas) ->";};
    double P(double e) { return e/3.0; }
    double T(double e) {
        double y = e/( (d_qgp/30.0)*PI*PI );
        return std::pow(y, 1.0 / 4.0);
    }
    double SV(double p, double e) {
        return e <= 0.0 ? 0. : std::sqrt( 1.0/3.0 );
    }
    double getTFromS(double s){
        double a = (2.0 * d_qgp * PI*PI /45.0 );
        return std::pow(s/a, 1.0 / 3.0);
    }
    double E(double T) {
        return d_qgp*PI*PI*T*T*T*T/30.0;
    }
    double S(double T) {
        return d_qgp*2.0*PI*PI*T*T*T/45.0;
    }
    double getEFromS(double s){
        return E(getTFromS(s));
    }
    std::string getName(){return "Massless Ideal Gas";}
};

class BMW: public EOS{
    private:

    double e1[dim1],e2[dim2];
    double t1[dim1],t2[dim2];
    double p1[dim1],p2[dim2];
    double c1[dim1],c2[dim2];
    double s1[dim1],s2[dim2];

    void Load();
    void Load1();
    void Load2();

    double getValueFromE(double e, double* a1, double* a2);
    int getIndexOf(double val, double bound, double* a1, double* a2);
    double getValueFrom(double val, double bound,
                        double* f1, double* f2, double* a1, double* a2);

    double getSBound() { return getValueFromE(ebound, s1,s2); }
    double getTBound() { return getValueFromE(ebound, t1,t2); }

    double getValueFromS(double s, double* a1, double* a2)
    {return getValueFrom(s,getSBound(),s1,s2,a1,a2); }

    double getValueFromT(double t, double* a1, double* a2)
    { return getValueFrom(t,getTBound(),t1,t2,a1,a2); }

    public:
    BMW();
    ~BMW(){JSINFO << "<-[PPM] Deleting EoS (BMW) ->";};
    double P(double e) { return getValueFromE(e,p1,p2); }
    double T(double e) { return getValueFromE(e,t1,t2); }
    double SV(double p, double e) { return getValueFromE(e,c1,c2); }
    double getTFromS(double s){ return getValueFromS(s, t1, t2); }
    double E(double T) { return getValueFromT(T, e1, e2); }
    double S(double T) { return getValueFromT(T,s1,s2); }
    double getEFromS(double s){ return E(getTFromS(s)); }
    std::string getName(){return "BMW";}
};

#endif

