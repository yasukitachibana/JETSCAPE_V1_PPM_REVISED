#ifndef FLUID_VAL_H_
#define FLUID_VAL_H_

#include "PPMutil.h"
#include "PPMdata.h"
#include "PPMgrid.h"
#include "PPMeos.h"
#include "PPMcoord.h"

#include <cfloat>

class FluidValuables{
    
private:
    
    int grid_nx, grid_ny, grid_neta;
    
    std::shared_ptr<EOS> eos;
    std::shared_ptr<Coordinates> coord;
    InitData &DATA;
    SCGrid &arena;
    
    int max_times;
    double err;
    double miss_energy;
    double U0standard;
    
    double initTotalE, initTotalPx, initTotalPy, initTotalPz;
    
    double (FluidValuables::*U_0)( double e, double p, double u0);
    double (FluidValuables::*U_i)( double e, double p, double u0, double ux);
    double (FluidValuables::*U_c)( double rhob, double u0 );
    double (FluidValuables::*ZeroComp)( double U0, double U3, double eta);
    double (FluidValuables::*ThirdComp)( double U0, double U3, double eta);
    double (FluidValuables::*Difference)(double U, double vx, double p);
    double (FluidValuables::*GetEfromU)( double Uabs, double U0, double vtilde);
    std::string (FluidValuables::*AskDirection3)( int direction );
    std::string (FluidValuables::*AskDirection5)( int d5 );

    
    double U_tau( double e, double p, double u0);
    double U_x( double e, double p, double u0, double ux);
    double U_current( double rhob, double u0);
    double TauComp( double U0, double U3, double eta);
    double EtaComp( double U0, double U3, double eta);
    double Diff(double Uabs, double U0, double x);
    double EfromU( double Uabs, double U0, double vtilde);
    std::string AskDirection3TauEta( int direction );
    std::string AskDirection5TauEta( int d5 );

    
    double U_tCart( double e, double p, double u0);
    double U_xCart( double e, double p, double u0, double ux);
    double U_currentCart( double rhob, double u0);
    double TComp( double U0, double U3, double eta);
    double ZComp( double U0, double U3, double eta);
    double DiffCart(double Uabs, double U0, double x);
    double EfromUCart( double Uabs, double U0, double vtilde);
    std::string AskDirection3Cart( int direction );
    std::string AskDirection5Cart( int d5 );
    
    void CopyArena( std::array<int, 3> i_copy, std::array<int, 3> i_og);

    void CalcThermalVal( const std::array<double, 5> &U, const double Uabs,
                         std::array<double, 4> &u,
                         double &e, double &rhob, double &p, double &temp );
    
public:
    FluidValuables(std::shared_ptr<EOS> eos_in,
                   std::shared_ptr<Coordinates> coord_in,
                   InitData &DATA_in,
                   SCGrid &arena_in);
    ~FluidValuables();// destructor
    
    void SetCartesian();
    void SetTauEta();
    
    void SetInitialProfile();
    
    void SetPreviousThermalVal();
    void GetTotalConservedQuantities( int init );
    
    void SetThermalVal( const std::array<int, 3> &i );
    void GetThermalVal( const std::array<double, 5> &U,
                        std::array<double, 4> &u,
                        double &e, double &rhob, double &p, double &temp );
    
    void SetBoundary();
    double SolveV(double Uabs, double U0);
    double Direction(double Ux, double Uabs );
    
    double GetU_0( double e, double p, double u0){
        return (this->*U_0)( e, p, u0);
    }
    double GetU_i( double e, double p, double u0, double ux){
        return (this->*U_i)( e, p, u0, ux);
    }
    double GetU_c( double rhob, double u0 ){
        return (this->*U_c)( rhob, u0 );
    }
    double GetZeroComp( double U0, double U3, double eta){
        return (this->*ZeroComp)( U0, U3, eta);
    }
    double GetThirdComp( double U0, double U3, double eta){
        return (this->*ThirdComp)( U0, U3, eta);
    }
    
    double SetEfromU( double Uabs, double U0, double vtilde ){
        return (this->*GetEfromU)( Uabs, U0, vtilde);
    }
    
    double abs_vector3(double x,double y, double z){
            return sqrt( x*x+y*y+z*z );
    }
    
    std::string AskDirection3dim( int direction ){
        return (this->*AskDirection3)(direction);
    }
    
    std::string AskDirection( int d5 ){
        return (this->*AskDirection5)(d5);
    }
    
};

#endif 
