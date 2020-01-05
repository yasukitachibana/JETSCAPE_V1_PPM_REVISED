#include "PPMcoord.h"

#include "JetScapeLogger.h"
using namespace Jetscape;

Coordinates::Coordinates(InitData &DATA_in):DATA(DATA_in){
    JSINFO << "<-[PPM] Creating Coordinates ->";
    grid_nx = DATA.nx;
    grid_ny = DATA.ny;
    grid_neta = DATA.neta;
}


Coordinates::~Coordinates() {
    JSINFO << "<-[PPM] Deleting Coordinates ->";
}// destructor

void Coordinates::SetCartesian(){

    DATA.tau0 = 0.0;
    tau0 = DATA.tau0/hbarc;
    tau = tau0;
    JSINFO << "<-[PPM] NOTICE! Initial lab time in lab. frame is set to "<< tau0*hbarc <<" fm ->";
    
    dtau = DATA.delta_tau/hbarc;
    dx[0] = DATA.delta_x/hbarc; // space step x [1/GeV]
    dx[1] = DATA.delta_y/hbarc; // space step y [1/GeV]
    dx[2] = DATA.delta_eta/hbarc; // space step z [1/GeV]

    SetInitEnd();
    
}

void Coordinates::SetTauEta(){
    
    tau0 = DATA.tau0/hbarc;
    tau = tau0;
    JSINFO << "<-[PPM] Initial proper time is set to "<< tau0*hbarc <<" fm ->";
    dtau = DATA.delta_tau/hbarc;
    dx[0] = DATA.delta_x/hbarc; // space step x [1/GeV]
    dx[1] = DATA.delta_y/hbarc; // space step y [1/GeV]
    dx[2] = DATA.delta_eta;

    SetInitEnd();
    
}

void Coordinates::SetInitEnd(){

    dV = dx[0]*dx[1]*dx[2];
    dxtilde[0] = dx[0];
    dxtilde[1] = dx[1];
    dxtilde[2] = dx[2];
    
}

void Coordinates::CountTau( const double i ){
    tau += i*dtau;
}

void Coordinates::CountTau( const int i ){
    tau += double(i)*dtau;
}

double Coordinates::GetX( const int ix){
    return ( ix - 0.5*(grid_nx-1) ) * dx[0];
}

double Coordinates::GetY( const int iy){
    return ( iy - 0.5*(grid_ny-1) ) * dx[1];
}

double Coordinates::GetEta( const int ieta){
    return ( ieta - 0.5*(grid_neta-1) ) * dx[2];
}

int Coordinates::GetIx( const double x){
    return int(x/dx[0] + grid_nx/2.0);
}

int Coordinates::GetIy( const double y){
    return int(y/dx[1] + grid_ny/2.0);
}

int Coordinates::GetIeta( const double eta){
    return int(eta/dx[2] + grid_neta/2.0);
}





