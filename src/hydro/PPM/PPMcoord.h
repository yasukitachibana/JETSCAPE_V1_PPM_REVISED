#ifndef COORD_H
#define COORD_H

#include "PPMdata.h"
#include "PPMutil.h"

#include <array>

class Coordinates
{
private:

    InitData &DATA;
    
    int grid_nx, grid_ny, grid_neta;
    
    void SetInitEnd();






public:
    Coordinates(InitData &DATA);
    ~Coordinates();  //destructor
    
    void SetCartesian();
    void SetTauEta();

    double tau0;// in [GeV^-1]
    double tau;// in [GeV^-1]
    double dtau;// in [GeV^-1]
    std::array<double, 3> dx;// in [GeV^-1]
    std::array<double, 3> dxtilde;// in [GeV^-1]
    double dV;
    
    double GetX( const int ix);
    double GetY( const int iy);
    double GetEta(const int ieta);
    
    double GetZeroComp( double U0, double U3, double eta);
    double GetThirdComp( double U0, double U3, double eta);

    void CountTau( const double i );
    void CountTau( const int i );
    
    int GetIx( const double x);
    int GetIy( const double y);
    int GetIeta( const double eta);
    
};

#endif
