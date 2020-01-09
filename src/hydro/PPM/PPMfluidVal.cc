#include "PPMfluidVal.h"
#include "JetScapeLogger.h"


using namespace Jetscape;
using namespace std;

FluidValuables::FluidValuables(std::shared_ptr<EOS> eos_in,
                               std::shared_ptr<Coordinates> coord_in,
                               InitData &DATA_in,
                               SCGrid &arena_in
                               ):DATA(DATA_in), arena(arena_in){
    
    JSINFO << "<-[PPM] Creating FluidValuables ->";
    
    grid_nx = DATA.nx;
    grid_ny = DATA.ny;
    grid_neta = DATA.neta;
    
    eos = eos_in;
    coord = coord_in;
    
    max_times = 50;
    err = 1.0e-6;
    
    miss_energy = 0.0;
    U0standard = eos->E(0.04);
    
}

FluidValuables::~FluidValuables(){
    JSINFO << "<-[PPM] Deleting FluidValuables ->";
}// destructor

void FluidValuables::SetCartesian(){
    U_0 = &FluidValuables::U_tCart;
    U_i = &FluidValuables::U_xCart;
    U_c = &FluidValuables::U_currentCart;
    ZeroComp = &FluidValuables::TComp;
    ThirdComp = &FluidValuables::ZComp;
    Difference = &FluidValuables::DiffCart;
    GetEfromU = &FluidValuables::EfromUCart;
}

void FluidValuables::SetTauEta(){
    U_0 = &FluidValuables::U_tau;
    U_i = &FluidValuables::U_x;
    U_c = &FluidValuables::U_current;
    ZeroComp = &FluidValuables::TauComp;
    ThirdComp = &FluidValuables::EtaComp;
    Difference = &FluidValuables::Diff;
    GetEfromU = &FluidValuables::EfromU;
}


void FluidValuables::SetInitialProfile(){
    JSINFO << "<-[PPM] Set Initial Fluid Valuables ->";
    for (int ix = 0; ix < grid_nx; ix++) {
        for (int iy = 0; iy < grid_ny; iy++) {
            for (int ieta = 0; ieta < grid_neta; ieta++) {
                
                arena(ix,iy,ieta).U[0]
                = (this->*U_0)( arena(ix,iy,ieta).epsilon,
                               arena(ix,iy,ieta).p,
                               arena(ix,iy,ieta).u[0] );
                
                for( int d=0; d<3; d++){
                    arena(ix,iy,ieta).U[d+1]
                    = (this->*U_i)( arena(ix,iy,ieta).epsilon,
                                   arena(ix,iy,ieta).p,
                                   arena(ix,iy,ieta).u[0],
                                   arena(ix,iy,ieta).u[d+1] );
                }
                
                arena(ix,iy,ieta).U[4]
                = (this->*U_c)( arena(ix,iy,ieta).rhob,
                               arena(ix,iy,ieta).u[0] );
                
                arena(ix,iy,ieta).U_prev = arena(ix,iy,ieta).U;
                
            }
        }
    }
}

void FluidValuables::GetTotalConservedQuantities( int init ){
    JSINFO << "<-[PPM] Show Conserved Quantities ->";
    double TotalE = 0.0;
    double TotalPx = 0.0;
    double TotalPy = 0.0;
    double TotalPz = 0.0;
    
    for (int ieta = 0; ieta < grid_neta; ieta++) {
        double eta = coord->GetEta(ieta);
        for (int ix = 0; ix < grid_nx; ix++) {
            for (int iy = 0; iy < grid_ny; iy++) {
                
                TotalE +=
                (this->*ZeroComp)(arena(ix,iy,ieta).U[0],
                                  arena(ix,iy,ieta).U[3],
                                  eta)*coord->dV;
                TotalPx += arena(ix,iy,ieta).U[1]*coord->dV;
                TotalPy += arena(ix,iy,ieta).U[2]*coord->dV;
                
                TotalPz +=
                (this->*ThirdComp)(arena(ix,iy,ieta).U[0],
                                   arena(ix,iy,ieta).U[3],
                                   eta)*coord->dV;
                
            }
        }
    }
    
    JSINFO
    << "<-[PPM]  tau (t) = " << coord->tau*hbarc << " fm/c, "
    << "TotalE = " << TotalE << " GeV, "
    << "TotalPx = " << TotalPx << " GeV/c, "
    << "TotalPy = " << TotalPy << " GeV/c, "
    << "TotalPz = " << TotalPz << " GeV/c. ->";
    
    
    if( init == 1 ){
        initTotalE = TotalE;
        initTotalPx = TotalPx;
        initTotalPy = TotalPy;
        initTotalPz = TotalPz;
    }
}

void FluidValuables::SetPreviousThermalVal(){
    JSINFO << "<-[PPM] Show Conserved Quantities ->" ;
    
    double TotalE = 0.0;
    double TotalPx = 0.0;
    double TotalPy = 0.0;
    double TotalPz = 0.0;
    
    for (int ieta = 0; ieta < grid_neta; ieta++) {
        double eta = coord->GetEta(ieta);
        for (int ix = 0; ix < grid_nx; ix++) {
            for (int iy = 0; iy < grid_ny; iy++) {
                
                //arena(ix,iy,ieta).epsilon_prev = arena(ix,iy,ieta).epsilon;
                arena(ix,iy,ieta).T_prev = arena(ix,iy,ieta).T;
                arena(ix,iy,ieta).U_prev = arena(ix,iy,ieta).U;
                //arena(ix,iy,ieta).u_prev = arena(ix,iy,ieta).u;
                
                TotalE +=
                (this->*ZeroComp)(arena(ix,iy,ieta).U[0],
                                  arena(ix,iy,ieta).U[3],
                                  eta)*coord->dV;
                TotalPx += arena(ix,iy,ieta).U[1]*coord->dV;
                TotalPy += arena(ix,iy,ieta).U[2]*coord->dV;
                TotalPz +=
                (this->*ThirdComp)(arena(ix,iy,ieta).U[0],
                                   arena(ix,iy,ieta).U[3],
                                   eta)*coord->dV;
                
            }
        }
    }
    
    JSINFO << "<-[PPM] tau (t) = " << coord->tau*hbarc << " fm/c ->";
    
    JSINFO
    << "<-[PPM] TotalE = " << TotalE << " GeV, "
    << "TotalPx = " << TotalPx << " GeV/c, "
    << "TotalPy = " << TotalPy << " GeV/c, "
    << "TotalPz = " << TotalPz << " GeV/c. ->";
    
    JSINFO
    << "<-[PPM] DeltaE = " << TotalE - initTotalE << " GeV, "
    << "DeltaPx = " << TotalPx - initTotalPx << " GeV/c, "
    << "DeltaPy = " << TotalPy - initTotalPy<< " GeV/c, "
    << "DeltaPz = " << TotalPz - initTotalPz<< " GeV/c. ->";
    
    if(DATA.profileType == 1){
        JSINFO << "<-[PPM] Bjorken Test ->";
        int ixc = int(grid_nx/2.0);
        int iyc = int(grid_ny/2.0);
        int ietac = int(grid_neta/2.0);
        double t_num
        = arena(ixc,iyc,ietac).T;
        double t_ana
        = (DATA.T0)*pow((coord->tau0/coord->tau),(1.0/3.0) );
        JSINFO
        << "<-[PPM] T_center = " << t_num
        << " GeV, T_analitic = " << t_ana
        << " GeV, deviation: " << 100*fabs(t_num - t_ana)/t_ana
        << " % ->";
    }else if(DATA.profileType == 2){
        int ixc = int(grid_nx/2.0);
        int iyc = int(grid_ny/2.0);
        int ietac = int(grid_neta/2.0);
        double t_num
        = arena(ixc,iyc,ietac).T;
        JSINFO
        << "<-[PPM] Dynamical Brick Test, T_center = " << t_num << " GeV ->";
    }else if(DATA.profileType == 3){
        int ixc = int(grid_nx/2.0);
        int iyc = int(grid_ny/2.0);
        int ietac = int(grid_neta/2.0);
        double t_num
        = arena(ixc,iyc,ietac).T;
        JSINFO << "<-[PPM] Dynamical Gaussian  Test, T_center = " << t_num << " GeV ->";
    }
    
}



double FluidValuables::U_tau( double e, double p, double u0){
    return coord->tau * ( u0*u0*(e+p) - p );
}

double FluidValuables::U_x( double e, double p, double u0, double ux){
    return coord->tau * u0 * ux * (e+p);
}

double FluidValuables::U_current( double rhob, double u0 ){
    return coord->tau * u0 * rhob;
}

double FluidValuables::U_tCart( double e, double p, double u0){
    return ( u0*u0*(e+p) - p );
}
double FluidValuables::U_xCart( double e, double p, double u0, double ux){
    return u0 * ux * (e+p);
}
double FluidValuables::U_currentCart( double rhob, double u0){
    return u0 * rhob;
}

double FluidValuables::TComp( double U0, double U3, double eta){
    return U0;
}
double FluidValuables::ZComp( double U0, double U3, double eta){
    return U3;
}

double FluidValuables::TauComp( double U0, double U3, double eta){
    return U0*cosh(eta) + U3*sinh(eta);
}
double FluidValuables::EtaComp( double U0, double U3, double eta){
    return U3*cosh(eta) + U0*sinh(eta);
}

double FluidValuables::Diff(double Uabs, double U0, double x){
    double e = EfromU( Uabs, U0, x);
    double p = eos->P( e );
    double w = x - ( Uabs / ( U0 + coord->tau * p) );
    return w;
}

double FluidValuables::DiffCart(double Uabs, double U0, double x){
    double e = EfromUCart( Uabs, U0, x);
    double p = eos->P( e );
    double w = x - Uabs / ( U0 + p);
    return w;
}

double FluidValuables::EfromU( double Uabs, double U0, double vtilde){
    return ( U0 - Uabs * vtilde ) / coord->tau;
}

double FluidValuables::EfromUCart( double Uabs, double U0, double v ){
    return ( U0 - Uabs * v );
}

void FluidValuables::SetBoundary(){
    //eta
    for (int ix = 3; ix < grid_nx-3; ix++) {
        for (int iy = 3; iy < grid_ny-3; iy++) {
            std::array<int, 3> i_og_l = {ix,iy,3};
            std::array<int, 3> i_og_h = {ix,iy,(grid_neta-4)};
            std::array<int, 3> i_copy = {ix,iy,0};
            for(int i = 0; i<3; i++){
                i_copy[2] = i;
                CopyArena( i_copy, i_og_l );
                i_copy[2] = grid_neta-i-1;
                CopyArena( i_copy, i_og_h );
            }
        }
    }
    //x
    for (int iy = 3; iy < grid_ny-3; iy++) {
        for (int ieta = 3; ieta < grid_neta-3; ieta++) {
            std::array<int, 3> i_og_l = {3,iy,ieta};
            std::array<int, 3> i_og_h = {(grid_nx-4),iy,ieta};
            std::array<int, 3> i_copy = {0,iy,ieta};
            for(int i = 0; i<3; i++){
                i_copy[0] = i;
                CopyArena( i_copy, i_og_l );
                i_copy[0] = grid_nx-i-1;
                CopyArena( i_copy, i_og_h );
            }
        }
    }
    //y
    for (int ix = 3; ix < grid_nx-3; ix++) {
        for (int ieta = 3; ieta < grid_neta-3; ieta++) {
            std::array<int, 3> i_og_l = {ix,3,ieta};
            std::array<int, 3> i_og_h = {ix,(grid_ny-4),ieta};
            std::array<int, 3> i_copy = {ix,0,ieta};
            for(int i = 0; i<3; i++){
                i_copy[1] = i;
                CopyArena( i_copy, i_og_l );
                i_copy[1] = grid_ny-i-1;
                CopyArena( i_copy, i_og_h );
            }
        }
    }
}

void FluidValuables::CopyArena( std::array<int, 3> i_copy, std::array<int, 3> i_og){
    arena(i_copy[0], i_copy[1], i_copy[2]) = arena(i_og[0], i_og[1], i_og[2]);
}

double FluidValuables::SolveV(double Uabs, double U0){

    double vmin = 0.0;
    double vmax = 1.0;
    double v;

    for(int i = 0; i < max_times; i++){
        v = ( vmin + vmax )/2.0;
        double diff = (this->*Difference)( Uabs, U0, v ); // B = v - f(v)
        if( fabs(diff) < err ){
            break;
        }else if( diff > 0.0 ){
            vmax = v;
        }else{
            vmin = v;
        }
    }
    return v;
}

double FluidValuables::Direction(double Ux, double Uabs ){
    return ( Uabs == 0. ? 0. : Ux / Uabs);
}

void FluidValuables::CalcThermalVal( const std::array<double, 5> &U, const double Uabs,
                                     std::array<double, 4> &u,
                                     double &e, double &rhob, double &p, double &temp ){
    
    double vabs =  Uabs<1e-80 ? 0.
    :SolveV( Uabs, U[0] );

    e = (this->*GetEfromU)( Uabs, U[0], vabs );
    rhob = 0.0;
    p = eos->P(e);
    temp = eos->T(e);

    u[0] = 1.0/sqrt(1.0-vabs*vabs);
    for( int d4=1; d4<4; d4++){
        u[d4] = u[0] * vabs * Direction( U[d4], Uabs );
    }
    
}

void FluidValuables::GetThermalVal( const std::array<double, 5> &U,
                                    std::array<double, 4> &u,
                                    double &e, double &rhob, double &p, double &temp ){

    double uabs2 = 0.0;
    for( int d=0; d<3; d++){
        uabs2 += U[d+1]*U[d+1];
    }
    double Uabs = sqrt(uabs2);
    
    CalcThermalVal( U, Uabs, u, e, rhob, p, temp );

}


void FluidValuables::SetThermalVal( const std::array<int, 3> &i ){
    
    double det = arena(i[0],i[1],i[2]).U[0]*arena(i[0],i[1],i[2]).U[0];
    double uabs2 = 0.0;
    for( int d=0; d<3; d++){
        uabs2 += arena(i[0],i[1],i[2]).U[d+1]*arena(i[0],i[1],i[2]).U[d+1];
    }
    
    det -= uabs2;
    
    if( arena(i[0],i[1],i[2]).U[0] <= 0.0 ){
        
        for( int d5=0; d5<5; d5++){
            arena(i[0],i[1],i[2]).U[d5] = 0.0;
        }
        
        arena(i[0],i[1],i[2]).epsilon = 0.0;
        arena(i[0],i[1],i[2]).rhob = 0.0;
        arena(i[0],i[1],i[2]).p = 0.0;
        arena(i[0],i[1],i[2]).T = 0.0;
        
        arena(i[0],i[1],i[2]).u[0] = 1.0;
        arena(i[0],i[1],i[2]).u[1] = 0.0;
        arena(i[0],i[1],i[2]).u[2] = 0.0;
        arena(i[0],i[1],i[2]).u[3] = 0.0;
        
    }else{

        double Uabs = sqrt(uabs2);

        if( det <= 0.0 ){
            //            double corf =
            //            //0.9999
            //            * arena(i[0],i[1],i[2]).U[0]
            //            /Uabs;
            //            for( int d5=1; d5<4; d5++){
            //                arena(i[0],i[1],i[2]).U[d5]
            //                *= corf;
            //            }
            //            Uabs
            //            = 0.9999
            //            * arena(i[0],i[1],i[2]).U[0];

            double U0new = Uabs/0.999999;

            if( arena(i[0],i[1],i[2]).U[0] > U0standard &&
                fabs(U0new - arena(i[0],i[1],i[2]).U[0])/arena(i[0],i[1],i[2]).U[0] > 0.95 ){

                miss_energy += (U0new - arena(i[0],i[1],i[2]).U[0])*coord->dV;

                if( miss_energy > 3.0){
                    JSINFO << "-Error- Space Like Fluid Cell. Spacelikeness is too large. ";
                    JSINFO << "U0=" << arena(i[0],i[1],i[2]).U[0];
                    JSINFO << "Uabs=" << Uabs;
                    JSINFO << "missing energy=" << miss_energy << " GeV";
                    JSINFO << "Check initial condition or Source terms. Exit.";
                    exit(-1);
                }
            }

            arena(i[0],i[1],i[2]).U[0]
            = U0new;

        }
        
        std::array<double, 4> u_get;
        double e_get, rhob_get, p_get, temp_get;
        
        CalcThermalVal( arena(i[0],i[1],i[2]).U, Uabs,
                        u_get, e_get, rhob_get, p_get, temp_get);

        arena(i[0],i[1],i[2]).u = u_get;
        arena(i[0],i[1],i[2]).epsilon = e_get;
        arena(i[0],i[1],i[2]).rhob = rhob_get;
        arena(i[0],i[1],i[2]).p = p_get;
        arena(i[0],i[1],i[2]).T = temp_get;

    }
}

