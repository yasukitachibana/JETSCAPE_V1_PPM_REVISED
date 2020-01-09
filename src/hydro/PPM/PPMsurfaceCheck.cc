#include "PPMsurfaceCheck.h"
#include "JetScapeLogger.h"

#include <string>


using namespace Jetscape;
using namespace std;

SurfaceCheck::SurfaceCheck( std::string filename_in ){
    
    JSINFO << "<-[PPM] Creating SurfaceCheck ->";
    
    filename = filename_in;
    
    total_energy = 0.0;
    total_momx   = 0.0;
    total_momy   = 0.0;
    total_momz   = 0.0;
    
}

SurfaceCheck::~SurfaceCheck(){
    JSINFO << "<-[PPM] Deleting SurfaceCheck ->";
}// destructor

void SurfaceCheck::DoCheck(){
    
    JSINFO << "<-[PPM] Start SurfaceCheck ->";
    
    std::ifstream surface_file( filename.c_str() );
    
    if(!surface_file)
    {
        cout << "<-[PPM] Surface file not present. ->" << endl;
        exit(1);
    }
    
    string text_string;
    double tau, x, y, eta, SU0, SU1, SU2, SU3, u0, u1, u2, u3;
    double eps, TFO, muB, epsPlusPOverT, W00, W01, W02, W03;
    double W11, W12, W13, W22, W23, W33, PiB;
    string sepsPlusPOverT;
    double crap;
    
    
    int counter = 0;
    getline(surface_file, text_string);
    
    while (!surface_file.eof()) {
        std::stringstream text_stream(text_string);
        text_stream >> tau >> x   >> y   >> eta
        >> SU0 >> SU1 >> SU2 >> SU3
        >> u0  >> u1  >> u2  >> u3
        >> eps >> TFO >> muB >> sepsPlusPOverT
        >> W00 >> W01 >> W02 >> W03
        >> W11 >> W12 >> W13 >> W22
        >> W23 >> W33
        //            >> PiB;
        >> crap >> crap >> crap >> crap >> crap >> crap;

        //cout << "sepsPlusPOverT= " << sepsPlusPOverT;
        if (sepsPlusPOverT=="-nan" || sepsPlusPOverT=="nan" ) epsPlusPOverT=0.;
        else epsPlusPOverT=atof(sepsPlusPOverT.c_str());
        //cout << epsPlusPOverT << endl;
        

        double pressure = epsPlusPOverT*TFO - eps;//+PiB
        double T_tau_tau = (eps+pressure)*u0*u0 - pressure + W00;
        double T_tau_x   = (eps+pressure)*u0*u1            + W01;
        double T_tau_y   = (eps+pressure)*u0*u2            + W02;
        double T_tau_eta = (eps+pressure)*u0*u3            + W03;
        double T_x_x     = (eps+pressure)*u1*u1 + pressure + W11;
        double T_x_y     = (eps+pressure)*u1*u2            + W12;
        double T_x_eta   = (eps+pressure)*u3*u1            + W13;
        double T_y_y     = (eps+pressure)*u2*u2 + pressure + W22;
        double T_y_eta   = (eps+pressure)*u3*u2            + W23;
        double T_eta_eta = (eps+pressure)*u3*u3 + pressure + W33;

        double cheta = cosh(eta);
        double sheta = sinh(eta);

        //Here we are writing, T_t_tau, T_t_x, T_t_y and T_t_eta
        double T00 = T_tau_tau*cheta + T_tau_eta*sheta;
        double T01 = T_tau_x  *cheta + T_x_eta  *sheta;
        double T02 = T_tau_y  *cheta + T_y_eta  *sheta;
        double T03 = T_tau_eta*cheta + T_eta_eta*sheta;

        double T10 = T_tau_x;
        double T11 = T_x_x  ;
        double T12 = T_x_y  ;
        double T13 = T_x_eta;

        double T20 = T_tau_y;
        double T21 = T_x_y  ;
        double T22 = T_y_y  ;
        double T23 = T_y_eta;

        //Here we are writing, T_z_tau, T_z_x, T_z_y and T_z_eta
        double T30 = T_tau_eta*cheta + T_tau_tau*sheta;
        double T31 = T_x_eta  *cheta + T_tau_x  *sheta;
        double T32 = T_y_eta  *cheta + T_tau_y  *sheta;
        double T33 = T_eta_eta*cheta + T_tau_eta*sheta;

        double energy = 0.0;
        double momx   = 0.0;
        double momy   = 0.0;
        double momz   = 0.0;


        energy = T00*SU0 + T01*SU1 + T02*SU2 + T03*SU3;
        momx   = T10*SU0 + T11*SU1 + T12*SU2 + T13*SU3;
        momy   = T20*SU0 + T21*SU1 + T22*SU2 + T23*SU3;
        momz   = T30*SU0 + T31*SU1 + T32*SU2 + T33*SU3;

        total_energy += energy;
        total_momx   += momx  ;
        total_momy   += momy  ;
        total_momz   += momz  ;

        getline(surface_file, text_string);
        counter++;

    }
    surface_file.close();
    EndMessage();

}




void SurfaceCheck::EndMessage(){
    
    cout << "total energy = " << total_energy << endl;
    cout << "total momx   = " << total_momx << endl;
    cout << "total momy   = " << total_momy << endl;
    cout << "total momz   = " << total_momz << endl;
    
}

