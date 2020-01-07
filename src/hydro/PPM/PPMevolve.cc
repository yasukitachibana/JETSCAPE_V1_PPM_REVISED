#include <math.h>

#include "PPMevolve.h"
#include "JetScapeLogger.h"

using namespace Jetscape;
using namespace std;

Evolve::Evolve( int run_num_in, std::shared_ptr<EOS> eos_in, InitData &DATA_in, SCGrid &arena_in ):DATA(DATA_in), arena(arena_in){

    run_num = run_num_in;
    eos = eos_in;
    
    grid_nx = DATA.nx;
    grid_ny = DATA.ny;
    grid_neta = DATA.neta;
    
    ppm_status = ppm_not_start;
    
    JSINFO<< "<-[PPM] Preparing Hydro Run: "<< run_num <<" ->";
    InitialSetting();
    
}

Evolve::~Evolve(){
    JSINFO << "<-[PPM] Deleting Evolve ->";
}// destructor

void Evolve::InitialSetting(){
    
    coord = std::make_shared<Coordinates>( DATA );
    fval = std::make_shared<FluidValuables>( eos, coord, DATA, arena );
    
    if( DATA.profileType == 2 || DATA.profileType == 3 ){
        JSINFO << "<-[PPM] Set Cartesian Coordinates ->";
        coord->SetCartesian();
        fval->SetCartesian();
        F_0 = &Evolve::F_tCart;
        F_i = &Evolve::F_xCart;
        Flux[0] = &Evolve::FluxTrans;
        Flux[1] = &Evolve::FluxTrans;
        Flux[2] = &Evolve::FluxTrans;
    }else{
        JSINFO << "<-[PPM] Tau-Eta Coordinates ->";
        coord->SetTauEta();
        fval->SetTauEta();
        F_0 = &Evolve::F_tau;
        F_i = &Evolve::F_x;
        Flux[0] = &Evolve::FluxTrans;
        Flux[1] = &Evolve::FluxTrans;
        Flux[2] = &Evolve::FluxTransEta;
    }
    
    fval->SetInitialProfile();
    fval->GetTotalConservedQuantities( 1 );
    
    printer = std::unique_ptr<Printer>
    (new Printer( run_num, eos, coord, DATA, arena ));
    
    freezeout = std::unique_ptr<Freezeout>
    (new Freezeout( run_num, coord, DATA, arena ));
    
    if( DATA.source==1 ){
        liquefier = std::unique_ptr<Liquefier>
        (new Liquefier( run_num, coord, fval, DATA, arena ));
    }else if(DATA.source==2){
        source_gauss = std::unique_ptr<SourceGauss>
        (new SourceGauss( run_num, coord, fval, DATA, arena ));
    }


}

void Evolve::EvolveIt(){
    
    int itmax = DATA.nt;
    std::array<int, 3> evoDir = {0,1,2};
    

    //----------------------------------------------------------------
    for (int it = 0; it <= itmax; it++) {

        ppm_status = freezeout->FindFreezeoutSurface(ppm_status);
        printer->PrintProfile(ppm_status);
        if( ppm_status == ppm_finished ){
            break;
        }
        fval->SetPreviousThermalVal();
        // evolution for this time step-----

        JSINFO << "<-[PPM] ######################################################################################### ->";
        JSINFO
        << "<-[PPM] Starting time step: "<< it << "/" << itmax
        <<", tau (t): "<<coord->tau*hbarc<<"->"<<(coord->tau+coord->dtau)*hbarc
        <<" fm ->";
        ppm_status = ppm_running;
        StepDtau( it, evoDir );
    } /* it */
    //----------------------------------------------------------------
    ppm_status = ppm_finished;
    ppm_status = freezeout->FindFreezeoutSurface(ppm_status);

    JSINFO
    << "<-[PPM] Hydro Evolution Run "
    << run_num
    <<" Done. tau (t) = "
    << coord->tau*hbarc
    <<" fm/c ->";
    
}

void Evolve::StepDtau( int it, std::array<int, 3> &evoDir ){
    
    coord->CountTau(0.5);
    
    if(DATA.profileType == 2 || DATA.profileType == 3 ){
        MiniStep( 1.0, evoDir[0] );
        MiniStep( 1.0, evoDir[1] );
        MiniStep( 1.0, evoDir[2] );
        
        if( DATA.source==1 ){
            liquefier->AddSourceCart();
        }else if(DATA.source==2){
            JSINFO
            << "<-[PPM] Causal Source for Cartesiaqn coordinates is coming soon. Exit. ->";
            exit(-1);
        }

        DirExchange( it, evoDir );
    }else{
        coord->dxtilde[2] = coord->tau*coord->dx[2];
        MiniStep( 0.5, evoDir[0] );
        MiniStep( 1.0, evoDir[1] );
        MiniStep( 0.5, evoDir[0] );
        MiniStep( 1.0, evoDir[2] );
        
        if( DATA.source==1 ){
            JSINFO
            << "<-[PPM] Liquefier for Tau-Eta coordinates is coming soon. Exit. ->";
            exit(-1);
        }else if(DATA.source==2){
            source_gauss->AddSourceTauEta();
        }
        swap(evoDir[0], evoDir[1]);
    }
    
    coord->CountTau(0.5);
    
}


void Evolve::MiniStep( double strang, int direction){

    JSINFO << "<-[PPM] (" << (strang)  << "/3) Step ->";
    
    double dstrang = coord->dtau*strang;
    double dtdx = (dstrang/coord->dxtilde[direction]);
    
    std::array<int, 3> s = {0,0,0};
    std::array<int, 3> ncell;
    
    SetFluidEvoDir(direction, s, ncell);
    
    for( int k = 3; k < ncell[2]-3; k++ ){
        for( int j = 3; j < ncell[1]-3; j++ ){
            
            //0: i-2, 1:i-1, 2: i, 3: i+1, 4: i+2
            std::array<std::array<int, 3>, 5> i_evo;
            SetFluidEvoIndexSet(i_evo, 2, j, k, s);
            
            std::array<std::array<double, 5>, 5> U_mem;
            SetInitUmem(U_mem, i_evo);
            
            //U_L, R
            //0: i-1, 1: i
            std::array<std::array<double, 5>, 2> U_L;
            std::array<std::array<double, 5>, 2> U_R;
            
            //0: i-1, 1: i
            std::array<double, 2> v_L = {0.0, 0.0};
            std::array<double, 2> v_R = {0.0, 0.0};
            std::array<double, 2> c_L = {0.0, 0.0};
            std::array<double, 2> c_R = {0.0, 0.0};
            
            GetU_R( U_L[1],
                   U_mem[1], U_mem[2],
                   U_mem[0], U_mem[3] );
            
            GetU_R( U_R[1],
                   U_mem[2], U_mem[3],
                   U_mem[1], U_mem[4] );
            
            SetMono( U_R[1], U_L[1], U_mem[2] );
            
            double p_dummy;
            GetBorderValuables( U_L[1], direction,
                               v_L[1], c_L[1], p_dummy );
            
            GetBorderValuables( U_R[1], direction,
                               v_R[1], c_R[1], p_dummy );
            
            
            for( int i = 3; i < ncell[0]-2; i++ ){
                
                SetFluidEvoIndexSet(i_evo, i, j, k, s);
                
                SetUmem(i_evo[4],
                        U_mem, U_L, U_R, v_L, v_R, c_L, c_R);
                
                GetU_R( U_L[1],
                       U_mem[1], U_mem[2],
                       U_mem[0], U_mem[3] );
                GetU_R( U_R[1],
                       U_mem[2], U_mem[3],
                       U_mem[1], U_mem[4] );
                
                SetMono( U_R[1], U_L[1], U_mem[2] );
                
                GetBorderValuables( U_L[1], direction,
                                   v_L[1], c_L[1], p_dummy );
                
                GetBorderValuables( U_R[1], direction,
                                   v_R[1], c_R[1], p_dummy );
                
                std::array<double, 5> U_minus = {0.,0.,0.,0.,0.};
                std::array<double, 5> U_plus  = {0.,0.,0.,0.,0.};
                
                GetUbarR( U_minus, U_mem[1], U_L[0], U_R[0],
                         v_R[0], c_R[0], v_L[1], c_L[1], dtdx );
                
                GetUbarL( U_plus, U_mem[2], U_L[1], U_R[1],
                         v_L[1], c_L[1], v_R[0], c_R[0], dtdx );
                
                std::array<double, 5> flux = {0,0,0,0,0};
                GetFluxHlle(flux, U_minus, U_plus, direction);
                
                //0: i-2, 1:i-1, 2: i, 3: i+1, 4: i+2
                (this->*Flux[direction])(flux, i_evo[1], i_evo[2], dtdx);
                
            }
        }
    }
    fval->SetBoundary();
}

void Evolve::SetFluidEvoDir(int direction,
                            std::array<int, 3> &s,
                            std::array<int, 3> &ncell){
    if( direction == 0 ){
        s[0]=1;
        ncell[0]=grid_nx;
        ncell[1]=grid_ny;
        ncell[2]=grid_neta;
    }else if( direction == 1 ){
        s[1]=1;
        ncell[0]=grid_ny;
        ncell[1]=grid_neta;
        ncell[2]=grid_nx;
    }else if( direction == 2 ){
        s[2]=1;
        ncell[0]=grid_neta;
        ncell[1]=grid_nx;
        ncell[2]=grid_ny;
    }
}

void Evolve::SetFluidEvoIndexSet( std::array<std::array<int, 3>, 5> &i_evo,
                                 int i, int j, int k,
                                 const std::array<int, 3> &s){
    for( int i_mem = 0; i_mem < 5; i_mem++ ){
        SetFluidEvoIndex(i_evo[i_mem], i-2+i_mem, j, k, s);
    }
}



void Evolve::SetFluidEvoIndex( std::array<int, 3> &i_evo,
                              int i, int j, int k,
                              const std::array<int, 3> &s){
    i_evo[0] = i*s[0]+j*s[2]+k*s[1];
    i_evo[1] = i*s[1]+j*s[0]+k*s[2];
    i_evo[2] = i*s[2]+j*s[1]+k*s[0];
}

void Evolve::SetInitUmem( std::array<std::array<double, 5>, 5> &U_mem,
                         const std::array<std::array<int, 3>, 5> &i_evo ){
    for(int i_mem=0; i_mem<5; i_mem++){
        U_mem[i_mem] =
        arena(i_evo[i_mem][0], i_evo[i_mem][1], i_evo[i_mem][2]).U;
    }
}


void Evolve::GetU_R( std::array<double, 5> &U_R,
                    const std::array<double, 5> &U,
                    const std::array<double, 5> &Up1,
                    const std::array<double, 5> &Um1,
                    const std::array<double, 5> &Up2 ){
    for( int d5=0; d5<5; d5++){
        double a = 0.5 * ( U[d5] + Up1[d5]);
        double b = Dmu( Up1[d5], U[d5], Um1[d5] );
        double c = Dmu( Up2[d5], Up1[d5], U[d5] );
        U_R[d5] = a + ( (b-c)/6.0 );
    }
}

double Evolve::Dmu(double Up1, double U, double Um1 ){
    double dmu;
    double difa = Up1 - U;
    double difb = U - Um1;
    double du = (Up1 - Um1)/2.0;
    
    if(difa*difb > 0.0) {
        double difaa = 2.0 * fabs(difa);
        double difba = 2.0 * fabs(difb);
        double dua = fabs(du);
        double c = min(dua, min(difaa, difba)) ;
        dmu = c * (du/dua) ;
    } else {
        dmu=0.0;
    }
    return dmu;
}

void Evolve::SetMono(std::array<double, 5> &UR,
                     std::array<double, 5> &UL,
                     std::array<double, 5> &UU ){
    for( int d5=0; d5<5; d5++){
        SetMonoComponent( UR[d5], UL[d5], UU[d5] );
    }
}

void Evolve::SetMonoComponent( double &UR, double &UL, double &UU){
    double difr = UR - UU;
    double difl = UU - UL;
    double difrl = UR - UL;
    double avrl = (UR+UL)/2.0;
    if(difr*difl  <=  0.0) {
        UR = UU;
        UL = UU;
    }else{
        if(difrl*(UU-avrl) >  difrl*difrl/6.0) {
            UL = 3.0*UU-2.0*UR;
        }else if(difrl*(UU-avrl) < -difrl*difrl/6.0 ) {
            UR = 3.0*UU-2.0*UL;
        }
    }
}

void Evolve::GetBorderValuables(std::array<double, 5> &U, int direction,
                                double &vx, double &c, double &p){
    double det = U[0]*U[0];
    for( int d=0; d<3; d++){
        det -= U[d+1] * U[d+1];
    }
    
    if( U[0] <= 0.0 || det <= 0.0){
        for( int d5=0; d5<5; d5++){
            U[d5] = 0.0;
        }
        vx = 0.0;
        c = 0.0;
        p = 0.0;
    }else{
        double Uabs = fval->abs_vector3( U[1], U[2], U[3] );
        double vabs = Uabs<1e-80 ? 0. :fval->SolveV( Uabs, U[0] );
        double e = fval->SetEfromU( Uabs, U[0], vabs );
        p = eos->P( e );
        c = eos->SV( p, e );
        vx = ( U[direction+1]==0. ? 0 :
              vabs * fval->Direction( U[direction+1], Uabs ));
    }
    
}


void Evolve::SetUmem(const std::array<int, 3> &i_evo,
                     std::array<std::array<double, 5>, 5> &U_mem,
                     std::array<std::array<double, 5>, 2> &U_L,
                     std::array<std::array<double, 5>, 2> &U_R,
                     std::array<double, 2> &v_L,
                     std::array<double, 2> &v_R,
                     std::array<double, 2> &c_L,
                     std::array<double, 2> &c_R ){
    
    for( int i_mem = 0; i_mem < 4; i_mem++ ){
        U_mem[i_mem] = U_mem[i_mem+1];
    }
    
    U_mem[4] = arena(i_evo[0], i_evo[1], i_evo[2]).U;
    
    U_L[0] = U_L[1];
    U_R[0] = U_R[1];
    
    v_L[0] = v_L[1];
    v_R[0] = v_R[1];
    
    c_L[0] = c_L[1];
    c_R[0] = c_R[1];
    
}

double Evolve::VelocitySum( double v1, double v2 ){
    return (v1 + v2) / (1.0 + (v1 * v2) );
}

void Evolve::GetUbarR( std::array<double, 5> &UbarR,
                      const std::array<double, 5> &U,
                      const std::array<double, 5> &UL,
                      const std::array<double, 5> &UR,
                      double v_r, double c_r,
                      double v_l_p, double c_l_p,
                      double dtdx){
    double v_h_r = 0.5*(v_r+v_l_p);
    double c_h_r = 0.5*(c_r+c_l_p);
    double lambda = max( 0.0,
                        max( VelocitySum(v_r, c_r),
                            VelocitySum(v_h_r, c_h_r) ) );
    lambda = fabs(lambda)*dtdx;
    
    for( int d5=0; d5<5; d5++){
        UbarR[d5] = GetUbarRComponent( U[d5], UL[d5],
                                      UR[d5], lambda );
    }
    
}

double Evolve::GetUbarRComponent(double U, double UL, double UR,
                                 double  lambda){
    double U6 = 6.0*( U-(UR+UL)/2.0 );
    double DelU = UR-UL;
    return UR-(lambda/2.0)*(DelU-(1.0-(2.0/3.0)*lambda)*U6);
}

void Evolve::GetUbarL( std::array<double, 5> &UbarL,
                      const std::array<double, 5> &U,
                      const std::array<double, 5> &UL,
                      const std::array<double, 5> &UR,
                      double v_l, double c_l,
                      double v_r_m, double c_r_m,
                      double dtdx){
    
    double v_h_l = 0.5*(v_l+v_r_m);
    double c_h_l = 0.5*(c_l+c_r_m);
    
    double lambda = max( 0.0,
                        max( VelocitySum( -v_l, c_l ),
                            VelocitySum( -v_h_l, c_h_l ) ) );
    lambda = fabs(lambda)*dtdx;
    
    for( int d5=0; d5<5; d5++){
        UbarL[d5] = GetUbarLComponent( U[d5], UL[d5],
                                      UR[d5], lambda );
    }
}

double Evolve::GetUbarLComponent(double U, double UL, double UR,
                                 double lambda){
    double U6 = 6.0 * ( U-(UR+UL)/2.0 );
    double DelU = UR - UL;
    return UL+(lambda/2.0)*(DelU+(1.0-(2.0/3.0)*lambda)*U6);
}

void Evolve::GetFluxHlle(std::array<double, 5> &flux,
                         std::array<double, 5> &U_minus,
                         std::array<double, 5> &U_plus,
                         int direction){
    
    double v_plus, v_minus;
    double c_plus, c_minus;
    double p_plus, p_minus;
    double b_plus, b_minus;
    
    GetBorderValuables(U_plus, direction,
                       v_plus,
                       c_plus,
                       p_plus);
    
    GetBorderValuables(U_minus, direction,
                       v_minus,
                       c_minus,
                       p_minus);
    
    double v_h = 0.5*(v_plus+v_minus);
    double c_h = 0.5*(c_plus+c_minus);
    
    b_plus  = GetBplus( v_h, c_h, v_plus, c_plus);
    b_minus = GetBminus( v_h, c_h, v_minus, c_minus);
    
    if( fabs(b_plus - b_minus) > 10e-30 ) {
        GetFh( flux, U_plus, U_minus,
              b_plus, b_minus, v_plus, v_minus,
              p_plus, p_minus, direction );
    }
}

double Evolve::GetBminus(double v_h, double c_h,
                         double v_minus, double c_minus){
    return min( 0.0,
               min( VelocitySum( v_h, -c_h ),
                   VelocitySum( v_minus, -c_minus ) ) );
}

double Evolve::GetBplus(double v_h, double c_h,
                        double v_plus, double c_plus){
    return max( 0.0,
               max( VelocitySum( v_h, c_h ),
                   VelocitySum( v_plus, c_plus ) ) );
}

void Evolve::GetFh(std::array<double, 5> &flux,
                   const std::array<double, 5> &U_plus,
                   const std::array<double, 5> &U_minus,
                   double b_plus, double b_minus,
                   double v_plus, double v_minus,
                   double p_plus, double p_minus,
                   int direction){
    
    std::array<double, 3> p_plus_v = {0.0, 0.0, 0.0};
    std::array<double, 3> p_minus_v = {0.0, 0.0, 0.0};
    p_plus_v[direction] = p_plus;
    p_minus_v[direction] = p_minus;
    
    flux[0] = GetFh0Component( U_plus[0], U_minus[0], b_plus, b_minus,
                              v_plus, v_minus, p_plus, p_minus );
    for( int d5=1; d5<4; d5++){
        flux[d5] = GetFhXComponent( U_plus[d5], U_minus[d5],
                                   b_plus, b_minus,
                                   v_plus, v_minus,
                                   p_plus_v[d5-1], p_minus_v[d5-1] );
    }
    
    
    flux[4] = GetFhCurrentComponent(U_plus[4], U_minus[4],
                                    b_plus, b_minus, v_plus, v_minus);
}

double Evolve::GetFhXComponent(double U_plus, double U_minus,
                               double b_plus, double b_minus,
                               double v_plus, double v_minus,
                               double p_plus, double p_minus){
    double f_plus = (this->*F_i)( U_plus, v_plus, p_plus );
    double f_minus = (this->*F_i)( U_minus, v_minus, p_minus );
    double a = b_plus * f_minus - b_minus * f_plus + b_plus * b_minus * ( U_plus - U_minus );
    return a/(b_plus - b_minus);
}

double Evolve::GetFh0Component(double U_plus, double U_minus,
                               double b_plus, double b_minus,
                               double v_plus, double v_minus,
                               double p_plus, double p_minus){
    double f_plus = (this->*F_0)( U_plus, v_plus, p_plus );
    double f_minus = (this->*F_0)( U_minus, v_minus, p_minus );
    double a = b_plus * f_minus - b_minus * f_plus + b_plus * b_minus * ( U_plus - U_minus );
    return a/(b_plus - b_minus);
}

double Evolve::GetFhCurrentComponent(double U_plus, double U_minus,
                                     double b_plus, double b_minus,
                                     double v_plus, double v_minus){
    double f_plus = Fcurrent( U_plus, v_plus);
    double f_minus = Fcurrent( U_minus, v_minus);
    double a = b_plus * f_minus - b_minus * f_plus + b_plus * b_minus * ( U_plus - U_minus );
    return a/(b_plus - b_minus);
}

double Evolve::F_x(double U, double vx, double p){
    return U * vx + coord->tau * p;
}

double Evolve::F_tau(double U, double vx, double p){
    return ( U + coord->tau * p ) * vx;
}

double Evolve::F_xCart(double U, double vx, double p){
    return U * vx + p;
}

double Evolve::F_tCart(double U, double vx, double p){
    return ( U + p ) * vx;
}

double Evolve::Fcurrent(double U, double vx){
    return U * vx;
}

void Evolve::FluxTransEta(std::array<double, 5> &flux,
                          const std::array<int, 3> &i_evo_m_1,
                          const std::array<int, 3> &i_evo,
                          double dtdx){
    
    double flux_component;
    for( int d5 = 1; d5 < 3; d5++ ){
        flux_component = dtdx * flux[d5];
        arena(i_evo_m_1[0], i_evo_m_1[1], i_evo_m_1[2]).U[d5]
        -= flux_component;
        arena(i_evo[0], i_evo[1], i_evo[2]).U[d5]
        += flux_component;
    }
    
    flux_component = dtdx * flux[4];
    
    arena(i_evo_m_1[0], i_evo_m_1[1], i_evo_m_1[2]).U[4]
    -= flux_component;
    arena(i_evo[0], i_evo[1], i_evo[2]).U[4]
    += flux_component;
    
    double eta_minus = coord->GetEta(i_evo_m_1[2]);
    double eta_plus = coord->GetEta(i_evo[2]);
    double eta_bound = 0.5*(eta_minus+eta_plus);
    
    //TauEta->Cartesian Lorentz transformation
    double flux_t = dtdx
    * (flux[0] * cosh(eta_bound) + flux[3] * sinh(eta_bound));
    double flux_z = dtdx
    * (flux[3] * cosh(eta_bound) + flux[0] * sinh(eta_bound));
    
    double flux_tau_minus = flux_t * cosh(eta_minus) - flux_z * sinh(eta_minus);
    double flux_eta_minus = - flux_t * sinh(eta_minus) + flux_z * cosh(eta_minus);
    
    double flux_tau_plus = flux_t * cosh(eta_plus) - flux_z * sinh(eta_plus);
    double flux_eta_plus = - flux_t * sinh(eta_plus) + flux_z * cosh(eta_plus);
    
    arena(i_evo_m_1[0], i_evo_m_1[1], i_evo_m_1[2]).U[0] -= flux_tau_minus;
    arena(i_evo[0], i_evo[1], i_evo[2]).U[0] += flux_tau_plus;
    
    arena(i_evo_m_1[0], i_evo_m_1[1], i_evo_m_1[2]).U[3] -= flux_eta_minus;
    arena(i_evo[0], i_evo[1], i_evo[2]).U[3] += flux_eta_plus;
    
    coord->CountTau(0.5);
    fval->SetThermalVal( i_evo_m_1 );
    coord->CountTau( - 0.5);
    
}

void Evolve::FluxTrans(std::array<double, 5> &flux,
                       const std::array<int, 3> &i_evo_m_1,
                       const std::array<int, 3> &i_evo,
                       double dtdx){
    
    double flux_component;
    for( int d5 = 0; d5 < 5; d5++ ){
        flux_component = dtdx * flux[d5];
        arena(i_evo_m_1[0], i_evo_m_1[1], i_evo_m_1[2]).U[d5]
        -= flux_component;
        arena(i_evo[0], i_evo[1], i_evo[2]).U[d5]
        += flux_component;
    }
    
    coord->CountTau( 0.5 );
    fval->SetThermalVal( i_evo_m_1 );
    coord->CountTau( -0.5 );
    
}

void Evolve::DirExchange( int it, std::array<int, 3> &evoDir ){
    if( (it+1) % 3 == 0  ){
        swap(evoDir[0], evoDir[1]);
    }else{
        int a = evoDir[0];
        int b = evoDir[1];
        int c = evoDir[2];
        evoDir[0] = b;
        evoDir[1] = c;
        evoDir[2] = a;
    }
}


//
//void Evolve::Print3dProfileBinary(Grid ***arena){
//    double t = tau * hbarc;
//    FILE *out_file;
//
//    string open_option;
//    if(evolving == 1){
//        open_option = "a";
//    }else{
//        open_option = "w";
//    }
//    open_option += "b";
//
//    JSINFO << "Printing 3D Profile (Binary Output File) at t=" << t << " fm/c";
//    //JSINFO <<  "open: " << open_option.c_str();
//
//    out_file = std::fopen((profile_filename_head+"_binary.dat").c_str(),
//                          open_option.c_str());
//
//    double hb3 = hbarc;
//    if( DATA_ptr->profileType != 2 && DATA_ptr->profileType != 3 ){
//        hb3 = 1.0;
//    }
//
//    //-------------------------------------------------
//    for (int ix = 3; ix < grid_nx-3; ix++) {
//        double x = GetX(ix)*hbarc;
//
//        for (int iy = 3; iy < grid_ny-3; iy++) {
//            double y = GetY(iy)*hbarc;
//
//            for (int ieta = 3; ieta < grid_neta-3; ieta++) {
//                double z = GetEta(ieta)*hbarc;
//
//                float cell_infos[]
//                ={
//                    static_cast<float> (t),//0
//                    static_cast<float> (x),//1
//                    static_cast<float> (y),//2
//                    static_cast<float> (z),//3
//                    static_cast<float> (arena(ix,iy,ieta).U[0]),//4
//                    static_cast<float> (arena(ix,iy,ieta).U[1]),//5
//                    static_cast<float> (arena(ix,iy,ieta).U[2]),//6
//                    static_cast<float> (arena(ix,iy,ieta).U[3]),//7
//                    static_cast<float> (arena(ix,iy,ieta).T),//8
//                    static_cast<float> (arena(ix,iy,ieta).epsilon),//9
//                    static_cast<float> (arena(ix,iy,ieta).p),//10
//                    static_cast<float> (arena(ix,iy,ieta).u[0]),//11
//                    static_cast<float> (arena(ix,iy,ieta).u[1]),//12
//                    static_cast<float> (arena(ix,iy,ieta).u[2]),//13
//                    static_cast<float> (arena(ix,iy,ieta).u[3]),//14
//                    static_cast<float> (dV)//15
//                };
//
//                //binary file
//                fwrite(cell_infos, sizeof(float), 16, out_file);
//
//            }//eta
//        }//y
//    }//x
//    fclose(out_file);
//}
//
//void Evolve::PrintRadialProfileBinary(Grid ***arena){
//    double t = tau * hbarc;
//    FILE *out_file;
//
//    string open_option;
//    if(evolving == 1){
//        open_option = "a";
//    }else{
//        open_option = "w";
//    }
//    open_option += "b";
//
//    JSINFO << "Printing Radial Profile (Binary Output File) at t=" << t << " fm/c";
//    //JSINFO <<  "open: " << open_option.c_str();
//
//    out_file = std::fopen((profile_filename_head+"_binary.dat").c_str(),
//                          open_option.c_str());
//
//    int iy = GetIy(0.0);
//    int iz = GetIeta(0.0);
//
//    //-------------------------------------------------
//    for (int ix = GetIx(0.0); ix < grid_nx-3; ix++) {
//        double x = GetX(ix)*hbarc;
//
//        float cell_infos[]
//        ={
//            static_cast<float> (t),//0
//            static_cast<float> (x),//1
//            static_cast<float> (arena(ix,iy,iz).epsilon/hbarc/hbarc/hbarc),//2
//            static_cast<float> (arena(ix,iy,iz).T),//3
//            static_cast<float> (eos->S( arena(ix,iy,iz).T )/hbarc/hbarc/hbarc),//4
//            static_cast<float> (arena(ix,iy,iz).p/hbarc/hbarc/hbarc),//5
//            static_cast<float> (arena(ix,iy,iz).u[0]),//6
//            static_cast<float> (arena(ix,iy,iz).u[1]),//7
//            static_cast<float> (arena(ix,iy,iz).u[2]),//8
//            static_cast<float> (arena(ix,iy,iz).u[3]),//9
//        };
//
//        //binary file
//        fwrite(cell_infos, sizeof(float), 10, out_file);
//
//    }//x
//    fclose(out_file);
//}
//
//
//void Evolve::Print3dProfile(Grid ***arena){
//
//    double t = tau * hbarc;
//
//    JSINFO << "Printing 3D Profile at t=" << t << " fm/c";
//    //profile_ofs.open(profile_output_file, std::ios::out);
//
//    profile_ofs.open(profile_filename_head+".txt", ios::out);
//
//    double hb3 = hbarc;
//    if( DATA_ptr->profileType != 2 && DATA_ptr->profileType != 3 ){
//        hb3 = 1.0;
//    }
//
//    int iyc = int(grid_ny/2.0);
//    //-------------------------------------------------
//    for (int ix = 3; ix < grid_nx-3; ix++) {
//        double x = GetX(ix)*hbarc;
//
//        for (int iy = 3; iy < grid_ny-3; iy++) {
//            double y = GetY(iy)*hbarc;
//
//        double z = 0.0;
//        int ieta = GetIeta(z);
//            // for (int ieta = 3; ieta < grid_neta-3; ieta++) {
//            //     double z = GetEta(ieta)*hb3;
//
//                //if(iy==iyc){
//                profile_ofs << std::defaultfloat;
//                profile_ofs
//                << t << " "
//                << x << " " << y << " "<< z << " "
//                << arena(ix,iy,ieta).U[0] << " "
//                << arena(ix,iy,ieta).U[1] << " "
//                << arena(ix,iy,ieta).U[2] << " "
//                << arena(ix,iy,ieta).U[3] << " "
//                << dV << " "
//                << arena(ix,iy,ieta).T << " "
//                << arena(ix,iy,ieta).epsilon << " "
//                << arena(ix,iy,ieta).p << " "
//                << arena(ix,iy,ieta).u[0] << " "
//                << arena(ix,iy,ieta).u[1] << " "
//                << arena(ix,iy,ieta).u[2] << " "
//                << arena(ix,iy,ieta).u[3] << "\n";
//                //}
//
//        //}//eta
//
//        }//y
////        profile_ofs
////        <<"\n";
//    }//x
//
//    profile_ofs << std::flush;
//    profile_ofs.close();
//}

 
