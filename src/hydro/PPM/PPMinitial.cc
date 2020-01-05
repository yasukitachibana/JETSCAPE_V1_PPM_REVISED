#include "PPMinitial.h"
#include <vector>
#include <iostream>

#include "JetScapeLogger.h"

using namespace std;
using namespace Jetscape;

Initial::Initial(std::shared_ptr<EOS> eos_in, InitData &DATA_in, SCGrid &arena_in ):DATA(DATA_in), arena(arena_in){
    eos = eos_in;
}

// destructor
Initial::~Initial() {
    JSINFO << "<-[PPM] Deleting Initial ->";
}

void Initial::InitArena( int run_num ) {
    
    if( run_num == 0 ){
        GenerateArena();
    }
    switch( DATA.profileType){
        case 0: InitProfileFromPreeq(); break;
        case 1: InitProfileBjorken(); break;
        case 2: InitProfileBjorken(); break;
        case 3: InitProfile3DGaussian(); break;
        case 4: InitProfile2DFromFile();
    }
    
}

void Initial::GenerateArena() {
    
    JSINFO << "<-[PPM] First Run, Generating Fluid Cell ->";
    
    if( DATA.profileType == 0 ){
        Reshape();
    }
        
    AddEdges();
    arena = SCGrid(DATA.nx, DATA.ny, DATA.neta);

}

void Initial::Reshape(){
    JSINFO << "<-[PPM] Reshaping Fluid Cell Configuration ->";
    DATA.delta_x = initial_dx;
    DATA.delta_y = initial_dx;
    if( fabs(initial_dz) < 0.000001 || initial_nz == 1){
        DATA.delta_eta = DATA.delta_tau;
    }else{
        DATA.delta_eta = initial_dz;
    }

    const int nx = static_cast<int>(sqrt(initial_energy_density.size()/initial_nz));
    const int ny = nx;

    extra_nx = 0;
    extra_ny = 0;
    extra_neta = 0;
    
    if( DATA.addCell == 0 ){
        DATA.nx = nx;
        DATA.ny = ny;
        DATA.neta = initial_nz;
    }else{
        
        if( DATA.nx > nx ){
            extra_nx = (DATA.nx - nx);
            if( extra_nx%2 != 0){
                extra_nx++;
                DATA.nx++;
            }
        }else{
            DATA.nx = nx;
        }
        if( DATA.ny > ny ){
            extra_ny = (DATA.ny - ny) ;
            if( extra_ny%2 != 0){
                extra_ny++;
                DATA.ny++;
            }
        }else{
            DATA.ny = ny;
        }
        if( DATA.neta > initial_nz ){
            extra_neta = (DATA.neta - initial_nz);
            if( extra_neta%2 != 0){
                extra_neta++;
                DATA.neta++;
            }
        }else{
            DATA.neta = initial_nz;
        }
        
    }
    
    DATA.x_size = DATA.delta_x*(DATA.nx-1);
    DATA.y_size = DATA.delta_y*(DATA.ny-1);
    DATA.eta_size = DATA.delta_eta*(DATA.neta-1);
    
}

void Initial::AddEdges( ){

    //Active Cell in PPM
    JSINFO << "<-[PPM] ACTIVE CELLS IN PPM ->";
    JSINFO << "neta(nz) = " << DATA.neta << ", nx = " << DATA.nx << ", ny = " << DATA.ny;
    JSINFO << "deta(dz)=" << DATA.delta_eta << " (fm), dx=" << DATA.delta_x << " fm, dy=" << DATA.delta_y << " fm";
    JSINFO << "x_size = "     << DATA.x_size << " fm, y_size = "   << DATA.y_size << " fm, eta_size(z_size) = " << DATA.eta_size << " (fm)";
    JSINFO << "--------------------------";
    //--
    DATA.nx += 6;
    DATA.ny += 6;
    DATA.neta += 6;
    DATA.x_size = DATA.delta_x*(DATA.nx-1);
    DATA.y_size = DATA.delta_y*(DATA.ny-1);
    DATA.eta_size = DATA.delta_eta*(DATA.neta-1);
    //--
    //All Cell in PPM
    JSINFO << "<-[PPM] ALL CELLS IN PPM ->";
    JSINFO << "neta(nz) = " << DATA.neta << ", nx = " << DATA.nx << ", ny = " << DATA.ny;
    JSINFO << "deta(dz)=" << DATA.delta_eta << " (fm), dx=" << DATA.delta_x << " fm, dy=" << DATA.delta_y << " fm";
    JSINFO << "x_size = "     << DATA.x_size << " fm, y_size = "   << DATA.y_size << " fm, eta_size(z_size) = " << DATA.eta_size << " (fm)";
    JSINFO << "--------------------------";

}

void Initial::InitProfileBjorken(){
    JSINFO << "<-[PPM] Generating Bjorken(tau-eta)/Brick(Cartesian) Initial Condition ->";

    const double temp = DATA.T0;
    const double epsilon = eos->E(temp);
    const double p = eos->P(epsilon);
    const double rhob = 0.0;
    
    for (int ieta = 0; ieta < DATA.neta; ieta++) {
        for (int ix = 0; ix < DATA.nx; ix++) {
            for (int iy = 0; iy < DATA.ny; iy++) {
                
                arena(ix,iy,ieta).epsilon = epsilon;
                arena(ix,iy,ieta).epsilon_prev = epsilon;
                arena(ix,iy,ieta).rhob = rhob;
                arena(ix,iy,ieta).p = p;
                arena(ix,iy,ieta).T = temp;
                arena(ix,iy,ieta).T_prev = temp;

                arena(ix,iy,ieta).u[0] = 1.0;
                arena(ix,iy,ieta).u[1] = 0.0;
                arena(ix,iy,ieta).u[2] = 0.0;
                arena(ix,iy,ieta).u[3] = 0.0;
                
                arena(ix,iy,ieta).u_prev[0] = 1.0;
                arena(ix,iy,ieta).u_prev[1] = 0.0;
                arena(ix,iy,ieta).u_prev[2] = 0.0;
                arena(ix,iy,ieta).u_prev[3] = 0.0;
                
            }
        }
    }
}

void Initial::InitProfile3DGaussian(){
    JSINFO << "<-[PPM] Generating 3-D Gaussian (Cartesian) Initial Condition ->";
    
    const double temp0 = DATA.T0;
    const double epsilon0 = eos->E(temp0);
    const double sig_gauss = 1.5;
    const double sig2 = sig_gauss*sig_gauss;
    const double r_max2 = 16.0*sig2;

    for (int ix = 3; ix < DATA.nx-3;  ix++) {
        const double x = GetX(ix);
        for (int iy = 3; iy < DATA.ny-3; iy++) {
            const double y = GetY(iy);
            for (int ieta = 3; ieta < DATA.neta-3; ieta++) {
                const double z = GetEta(ieta);
                
                const double r2 = x*x + y*y + z*z;
                
                
                if( r2 <  r_max2 ){
                    
                    double epsilon = epsilon0*exp(-r2/(2.0*sig2));
                    double temp = eos->T(epsilon);
                    double p = eos->P(epsilon);
                    
                    arena(ix,iy,ieta).epsilon = epsilon;
                    arena(ix,iy,ieta).epsilon_prev = epsilon;
                    arena(ix,iy,ieta).p = p;
                    arena(ix,iy,ieta).T = temp;
                    arena(ix,iy,ieta).T_prev = temp;
                    
                    //                    JSINFO
                    //                    << "|| r: "<<sqrt(r2)
                    //                    << "| x: "<<x
                    //                    << ", y: "<<y
                    //                    << ", z: "<<z
                    //                    << "| enegy density: "<<epsilon;
                }else{
                    
                    arena(ix,iy,ieta).epsilon = 0.0;
                    arena(ix,iy,ieta).epsilon_prev = 0.0;
                    arena(ix,iy,ieta).p = 0.0;
                    arena(ix,iy,ieta).T = 0.0;
                    arena(ix,iy,ieta).T_prev = 0.0;
                    
                }
                
                arena(ix,iy,ieta).rhob = 0.0;
                
                arena(ix,iy,ieta).u[0] = 1.0;
                arena(ix,iy,ieta).u[1] = 0.0;
                arena(ix,iy,ieta).u[2] = 0.0;
                arena(ix,iy,ieta).u[3] = 0.0;
                
                arena(ix,iy,ieta).u_prev[0] = 1.0;
                arena(ix,iy,ieta).u_prev[1] = 0.0;
                arena(ix,iy,ieta).u_prev[2] = 0.0;
                arena(ix,iy,ieta).u_prev[3] = 0.0;
                
                
                //                JSINFO
                //                << "| x: "<<x
                //                << ", y: "<<y
                //                << ", z: "<<z;
//
//                JSINFO
//                << "| "<<ix
//                << ", "<<iy
//                << ", "<<ieta;

                
                
            }
        }
    }
    
    SetBoundary( 0, 0, 0 );
}

void Initial::InitProfile2DFromFile(){
    JSINFO << "<-[PPM] Generating Initial Condition From File (2D) ->";
    JSINFO << "<-[PPM] Input profile filename: '"+ DATA.profile_input_file+"' ->";
    JSINFO << "<-[PPM] Load File Assuming nx = ny = 201, dx = dy = 0.1 fm for Input File ->";
    JSINFO << "<-[PPM] Initial tau = " << DATA.tau0 << "fm/c (= " << DATA.tau0/hbarc << "GeV^-1) ->";
    
    if( DATA.delta_x < 0.1 && DATA.delta_y < 0.1 ){
        JSWARN << "<-[PPM] dx, dy must be larger than 0.1 fm. ->";
        exit(-1);
    }
    
    
    std::array<std::array<double, n_x_input>, n_x_input> e_trans;
    LoadInitProfile2D( e_trans );
        
    for (int ieta = 3; ieta < DATA.neta-3; ieta++) {
        const double eta = GetEta(ieta);
        
        
        //             //rootS   C    alpha  sigmaNN   nucleus  eta_flat  eta_gauss
        //      RHIC = { 200,  15.0, 0.18,  4.194,      Au,      1.3,      2.1     };//values from [1204.5814]
        //      LHC  = { 2760, 41.4, 0.08,  6.136,      Pb,      1.9,      3.2     };//values from [1204.5814]
        // LHC_5TeV  = { 5020, 47.5, 0.08,  7.0,      Pb,      2.0,      3.1     };
        
        //double weight = 1.0;
        double eta_flat  = 2.0;
        double sigma_eta = 3.1;
        double e_cm = 5020;
        double weight = ProfileEtaFlatGauss( eta, eta_flat, sigma_eta, e_cm );
        
        
        for (int ix = 3; ix < DATA.nx-3; ix++) {
            const double x = GetX(ix);
            for (int iy = 3; iy < DATA.ny-3; iy++) {
                const double y = GetY(iy);
                
                if( (fabs(x) - size_input ) < DBL_EPSILON &&
                   (fabs(y) - size_input ) < DBL_EPSILON ){
                    
                    int ix_data = GetIx2DInput(x);
                    int iy_data = GetIy2DInput(y);

                    double epsilon = e_trans[ix_data][iy_data]*weight*(hbarc*hbarc*hbarc);
                    double temp = eos->T(epsilon);
                    double p = eos->P(epsilon);
                    arena(ix,iy,ieta).epsilon = epsilon;
                    arena(ix,iy,ieta).epsilon_prev = epsilon;
                    arena(ix,iy,ieta).p = p;
                    arena(ix,iy,ieta).T = temp;
                    arena(ix,iy,ieta).T_prev = temp;
                    
                }else{
                    
                    arena(ix,iy,ieta).epsilon = 0.0;
                    arena(ix,iy,ieta).epsilon_prev = 0.0;
                    arena(ix,iy,ieta).p = 0.0;
                    arena(ix,iy,ieta).T = 0.0;
                    arena(ix,iy,ieta).T_prev = 0.0;
                    
                }
                
                arena(ix,iy,ieta).rhob = 0.0;
                arena(ix,iy,ieta).u[0] = 1.0;
                arena(ix,iy,ieta).u[1] = 0.0;
                arena(ix,iy,ieta).u[2] = 0.0;
                arena(ix,iy,ieta).u[3] = 0.0;
                
                arena(ix,iy,ieta).u_prev[0] = 1.0;
                arena(ix,iy,ieta).u_prev[1] = 0.0;
                arena(ix,iy,ieta).u_prev[2] = 0.0;
                arena(ix,iy,ieta).u_prev[3] = 0.0;
                
            }
        }
        
    }
    
    SetBoundary( 0, 0, 0 );
    
}



void Initial::LoadInitProfile2D(std::array<std::array<double, n_x_input>, n_x_input> &e_trans ){
    
    std::ifstream ifs;
    std::stringstream str_stream;
    ifs.open(DATA.profile_input_file.c_str()); //open the input file
    str_stream << ifs.rdbuf(); //read the file
    if (ifs.is_open()){
        JSINFO << "<-[PPM] FOUND '" << DATA.profile_input_file << "' ->";
        ifs.close();
    }else{
        ifs.close();
        JSINFO << "<-[PPM] NOT FOUND '" << DATA.profile_input_file << "' ->";
        exit(-1);
    }

    std::string line;

    int n_total = 0;
    int n_y = 0;
    int n_x = 0;
    double dummy0, x, y, ed, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8;

    while ( getline( str_stream, line ) ){

        if( line.find("n_eta") != std::string::npos && n_total == 0){
            //JSINFO << "First line: " << line;
            continue;
        }

        if( line.empty() && n_total != 0 ){
            n_x++;
            n_y = 0;
            continue;
        }

        if( ! line.empty() ){

            sscanf( line.data(),
                   "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &dummy0, &x, &y, &ed, &dummy2, &dummy3, &dummy4,
                   &dummy5, &dummy6, &dummy7, &dummy8 );

            e_trans[n_x][n_y] = ed;

            n_y++;
            n_total++;
        }
    }

    JSINFO << "<-[PPM] LOADED. Total Data Line: " << n_total << " ->";
}

int Initial::GetIx2DInput(double x){
    return int( x/dx_input + 0.5*n_x_input );
}

int Initial::GetIy2DInput(double y){
    return int( y/dx_input + 0.5*n_x_input );
}

double Initial::ProfileEtaFlatGauss(double eta, double eta_flat, double sigma_eta, double e_cm){

  double weight=1.0;

  const double eta_abs = fabs(eta);
  const double eta_out = eta_abs - eta_flat;
  const double mass_nucleon = 0.939;
  const double y_beam = acosh(e_cm/mass_nucleon);

  if( eta_abs >= y_beam ){
    weight=0.0;
  }else if( eta_out >= 0.0 ){
    weight = exp( - ( eta_out*eta_out ) / ( sigma_eta*sigma_eta ) );
  }
  return weight;

}

void Initial::SetBoundary( int extra_nx, int extra_ny, int extra_neta ){
    //eta
    for (int ix = 3+(extra_nx/2); ix < DATA.nx-(extra_nx/2)-3; ix++) {
        for (int iy = 3+(extra_ny/2); iy < DATA.ny-(extra_ny/2)-3; iy++) {
            
            std::array<int, 3> i_og_l = {ix,iy,3+(extra_neta/2)};
            std::array<int, 3> i_og_h = {ix,iy,(DATA.neta-(extra_neta/2)-4)};
            std::array<int, 3> i_copy = {ix,iy,0};
            
            for(int i = 0; i<3+(extra_neta/2); i++){
                i_copy[2] = i;
                CopyThermo( i_copy, i_og_l );
                i_copy[2] = DATA.neta-i-1;
                CopyThermo( i_copy, i_og_h );
            }
            
        }
    }
    
    //x
    for (int iy = 3+(extra_ny/2); iy < DATA.ny-(extra_ny/2)-3; iy++) {
        for (int ieta = 0; ieta < DATA.neta; ieta++) {
            
            std::array<int, 3> i_og_l = {3+(extra_nx/2),iy,ieta};
            std::array<int, 3> i_og_h = {(DATA.nx-(extra_nx/2)-4),iy,ieta};
            std::array<int, 3> i_copy = {0,iy,ieta};
            
            for(int i = 0; i<3+(extra_nx/2); i++){
                i_copy[0] = i;
                CopyThermo( i_copy, i_og_l );
                i_copy[0] = DATA.nx-i-1;
                CopyThermo( i_copy, i_og_h );
            }
            
        }
    }
    
    //y
    for (int ix = 0; ix < DATA.nx; ix++) {
        for (int ieta = 0; ieta < DATA.neta; ieta++) {
            
            std::array<int, 3> i_og_l = {ix,3+(extra_ny/2),ieta};
            std::array<int, 3> i_og_h = {ix,(DATA.ny-(extra_ny/2)-4),ieta};
            std::array<int, 3> i_copy = {ix,0,ieta};
            
            for(int i = 0; i<3+(extra_ny/2); i++){
                i_copy[1] = i;
                CopyThermo( i_copy, i_og_l );
                i_copy[1] = DATA.ny-i-1;
                CopyThermo( i_copy, i_og_h );
            }
            
        }
    }

}

void Initial::CopyThermo( std::array<int, 3> i_copy, std::array<int, 3> i_og ){
    arena(i_copy[0], i_copy[1], i_copy[2]) = arena(i_og[0], i_og[1], i_og[2]);
}

double Initial::GetX(int ix){
    return ( ix - 0.5*(DATA.nx-1) ) * DATA.delta_x;
}

double Initial::GetY(int iy){
    return ( iy - 0.5*(DATA.ny-1) ) * DATA.delta_y;
}

double Initial::GetEta(int ieta){
    return ( ieta - 0.5*(DATA.neta-1) ) * DATA.delta_eta;
}


void Initial::InitProfileFromPreeq(){
    JSINFO << "<-[PPM] Generating Initial Condition from Pre-eq. Modlue ->";
    double tau_init = DATA.tau0/hbarc;
    JSINFO << "<-[PPM] Initial tau = " << DATA.tau0 << "fm/c (= " << tau_init << "GeV^-1) ->";

    for (int ieta = 3+(extra_neta/2); ieta < DATA.neta-(extra_neta/2)-3; ieta++) {
        for (int ix = 3+(extra_nx/2); ix < DATA.nx-(extra_nx/2)-3; ix++) {
            for (int iy = 3+(extra_ny/2); iy < DATA.ny-(extra_ny/2)-3; iy++) {

                const int idx
                = (ix-(extra_nx/2)-3)
                + ( (iy-(extra_ny/2)-3) + (ieta-(extra_neta/2)-3) * (DATA.ny-extra_ny-6) )
                * (DATA.nx-extra_nx-6);

                double epsilon = 0.0;
                const double rhob = 0.0;
                double p = 0.0;
                double temp = 0.0;
                epsilon = ( initial_energy_density[idx] * DATA.sFactor * hbarc*hbarc*hbarc );
                p = eos->P(epsilon);
                temp = eos->T(epsilon);

                arena(ix,iy,ieta).epsilon = epsilon;
                arena(ix,iy,ieta).epsilon_prev = epsilon;
                arena(ix,iy,ieta).rhob = rhob;
                arena(ix,iy,ieta).p = p;
                arena(ix,iy,ieta).T = temp;
                arena(ix,iy,ieta).T_prev = temp;

                arena(ix,iy,ieta).u[0] = initial_u_tau[idx];
                arena(ix,iy,ieta).u[1] = initial_u_x[idx];
                arena(ix,iy,ieta).u[2] = initial_u_y[idx];
                arena(ix,iy,ieta).u[3] = tau_init*initial_u_eta[idx];

                arena(ix,iy,ieta).u_prev[0] = initial_u_tau[idx];
                arena(ix,iy,ieta).u_prev[1] = initial_u_x[idx];
                arena(ix,iy,ieta).u_prev[2] = initial_u_y[idx];
                arena(ix,iy,ieta).u_prev[3] = tau_init*initial_u_eta[idx];

            }//y
        }//x
    }//eta
    
    clean_up_arrays();
    SetBoundary( extra_nx, extra_ny, extra_neta );
    
}


void Initial::get_preequilibrium_vectors(const double dx, const double dz,
                                         const double z_max, const int nz,
                                         vector<double> e_in,
                                         vector<double> u_tau_in, vector<double> u_x_in,
                                         vector<double> u_y_in,   vector<double> u_eta_in,
                                         vector<double> pi_00_in, vector<double> pi_01_in,
                                         vector<double> pi_02_in, vector<double> pi_03_in,
                                         vector<double> pi_11_in, vector<double> pi_12_in,
                                         vector<double> pi_13_in, vector<double> pi_22_in,
                                         vector<double> pi_23_in, vector<double> pi_33_in,
                                         vector<double> Bulk_pi_in) {
    initial_dx = dx;
    initial_dz = dz;
    initial_z_max = z_max;
    initial_nz = nz;
    initial_energy_density = e_in;
    initial_u_tau          = u_tau_in;
    initial_u_x            = u_x_in;
    initial_u_y            = u_y_in;
    initial_u_eta          = u_eta_in;
    initial_pi_00          = pi_00_in;
    initial_pi_01          = pi_01_in;
    initial_pi_02          = pi_02_in;
    initial_pi_03          = pi_03_in;
    initial_pi_11          = pi_11_in;
    initial_pi_12          = pi_12_in;
    initial_pi_13          = pi_13_in;
    initial_pi_22          = pi_22_in;
    initial_pi_23          = pi_23_in;
    initial_pi_33          = pi_33_in;
    initial_bulk_pi        = Bulk_pi_in;
}

void Initial::clean_up_arrays(){
    // clean up
    initial_energy_density.clear();
    initial_u_tau.clear();
    initial_u_x.clear();
    initial_u_y.clear();
    initial_u_eta.clear();
    initial_pi_00.clear();
    initial_pi_01.clear();
    initial_pi_02.clear();
    initial_pi_03.clear();
    initial_pi_11.clear();
    initial_pi_12.clear();
    initial_pi_13.clear();
    initial_pi_22.clear();
    initial_pi_23.clear();
    initial_pi_33.clear();
    initial_bulk_pi.clear();
}








