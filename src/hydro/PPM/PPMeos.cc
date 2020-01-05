#include "PPMeos.h"

#include "JetScapeLogger.h"

using namespace Jetscape;

BMW::BMW(){
    JSINFO << "<-[PPM] Equation of State: BMW ->";
    Load();
}

void BMW::Load(){
    Load1();
    Load2();
}

void BMW::Load1(){
    std::string filename = "../src/hydro/PPM/EOS/W_B_EoS1.dat";
    std::ifstream ifs(filename.c_str());
    if(ifs.fail()) {
        JSWARN << "<-[PPM] EoS (BMW,1) File - Not Found ->"; exit(0);
    }

    std::string str;
    while(getline(ifs, str)) {
        int iload;
        double e, t, p, c, s;
        std::sscanf(str.data(), "%d %lf %lf %lf %lf %lf",  &iload, &e, &t, &p, &c, &s);
        e1[iload] = e;
        t1[iload] = t;
        p1[iload] = p;
        c1[iload] = c;
        s1[iload] = s;
    }
}

void BMW::Load2(){
    std::string filename = "../src/hydro/PPM/EOS/W_B_EoS2.dat";
    std::ifstream ifs(filename.c_str());
    if(ifs.fail()) {
        JSWARN << "<-[PPM] EoS (BMW,2) File - Not Found ->"; exit(0);
    }
    std::string str;
    while(getline(ifs, str)) {
        int iload;
        double e, t, p, c, s;
        std::sscanf(str.data(), "%d %lf %lf %lf %lf %lf", &iload, &e, &t, &p, &c, &s);
        e2[iload] = e;
        t2[iload] = t;
        p2[iload] = p;
        c2[iload] = c;
        s2[iload] = s;
    }
}

double BMW::getValueFromE(double e, double* a1, double* a2) {

    double ret;
    const double xtab = e<ebound ? e/de1 : e/de2;
    const int itab = (int)xtab;
    const double res = xtab - itab;
    double* _a = e<ebound ? a1 : a2;
    ret = (_a[itab] + res * (_a[itab+1] - _a[itab]));
    return ret;
}

int BMW::getIndexOf(double val, double bound, double* a1, double* a2){
    int itab=0;
    double* _a = val< bound ? a1 : a2;
    while( val>_a[itab] ){ itab++; }//TODO: Use better search algorithm
    return itab;
}

double BMW::getValueFrom(double val, double bound,
                         double* f1, double* f2, double* a1, double* a2){
    double ret;
    const int itab = getIndexOf(val, bound, f1, f2);
    double* _f = val<bound ? f1 : f2;
    const double df = _f[itab+1] - _f[itab];
    const double res = ( val - _f[itab] )/df;
    double* _a = val< bound ? a1 : a2;
    ret = _a[itab] + res * (_a[itab+1] - _a[itab]);
    return ret;
}

