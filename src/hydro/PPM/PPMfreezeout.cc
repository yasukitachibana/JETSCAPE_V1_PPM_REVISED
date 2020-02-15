#include "PPMfreezeout.h"
#include "JetScapeLogger.h"


using namespace Jetscape;
using namespace std;

Freezeout::Freezeout( int run_num_in,
                     std::shared_ptr<Coordinates> coord_in,
                     std::shared_ptr<FluidValuables> fval_in,
                     std::shared_ptr<EOS> eos_in,
                     InitData &DATA_in,
                     SCGrid &arena_in
                     ):DATA(DATA_in), arena(arena_in){
    
    JSINFO << "<-[PPM] Creating Freezeout ->";
    
    run_num = run_num_in;
    coord = coord_in;
    fval = fval_in;
    eos = eos_in;
    
    temp_fo = DATA.temp_fo;
    
    grid_nx = DATA.nx;
    grid_ny = DATA.ny;
    grid_neta = DATA.neta;
    
    rapidity_window = DATA.rapidity_window;
    
    if(DATA.fo_type == 0){
        JSINFO << "<-[PPM] No Freezeout ->";
    }else{
        
        if(DATA.fo_type == 1){
            JSINFO << "<-[PPM] Freezeout Temperature: " << temp_fo <<" GeV ->";
        }else if(DATA.fo_type == 2){
            JSINFO << "<-[PPM] Isochronous Freezeout ->";
        }
        JSINFO << "<-[PPM] Freezeout for fluid cells with |eta_s| < "<< rapidity_window <<" ->";
        
        surface_filename = GenerateSurfaceFilename();
        OpenSurfaceFile();
        
    }
    
    
}

Freezeout::~Freezeout(){
}// destructor

std::string Freezeout::GenerateSurfaceFilename(){
    return DATA.fo_surface + "_run_" + to_string(run_num) + ".txt";
}

void Freezeout::OpenSurfaceFile(){
    JSINFO << "<-[PPM] Surface Info File Name: '" << surface_filename << "' ->";
    ofs_freezeout.open(surface_filename.c_str(), ios_base::out);
}

int Freezeout::FindFreezeoutSurface( int ppm_status ){
    
    if( DATA.fo_type == 1 && ppm_status == ppm_not_start ){
        //first step for isothermal
        IsochronousFreezeout( 1 );
        return ppm_status;
        
    }else if( DATA.fo_type == 1 && ppm_status == ppm_running ){
        //isothermal
        int status = FullFreezeout();
        return status;
    }else if( ppm_status == ppm_finished ){
        if( DATA.fo_type == 2 ){
            //isochronous
            IsochronousFreezeout( 0 );
        }
        ofs_freezeout.close();
        if( DATA.surface_check == 1 ){
            
            std::unique_ptr<SurfaceCheck> surf_check(new SurfaceCheck(surface_filename));
            surf_check->DoCheck();
            
        }
        return ppm_status;
    }else{
        return ppm_status;
    }
    
}

void Freezeout::IsochronousFreezeout( int threshold ){
    
    double t_cut = temp_fo;
    
    if( threshold == 1 ){
        JSINFO
        << "<-[PPM] Finding Isochronous Freezeout Surface for Fluid Cells with T < "
        << t_cut <<" GeV ->";
    }else{
        t_cut = 50000.0;
        JSINFO
        << "<-[PPM] Finding Isochronous Freezeout Surface for All Fluid Cells at tau = "<< coord->tau*hbarc <<" fm/c ->";
    }
    
    std::array<double, 4> dsigma = GetDsigma( 1 );
    
    for (int ieta = 3; ieta < grid_neta-3; ieta++) {
        double eta = coord->GetEta(ieta);
        if( fabs(eta) > rapidity_window ){continue;}
        for (int ix = 3; ix < grid_nx-3; ix++) {
            for (int iy = 3; iy < grid_ny-3; iy++) {
                
                std::array<int, 3> i_cell = { ix, iy, ieta };
                
                std::array<double, 3> x_cell =
                { coord->GetX(ix)*hbarc, coord->GetY(iy)*hbarc, eta };
                
                double temperature = arena(i_cell[0],i_cell[1],i_cell[2]).T;
                
                if( temperature > t_cut || arena(i_cell[0],i_cell[1],i_cell[2]).U[0] < DBL_MIN ){
                    continue;
                }
                
                
                double   x_surface = x_cell[0];
                double   y_surface = x_cell[1];
                double eta_surface = x_cell[2];
                double e_surface = arena(i_cell[0],i_cell[1],i_cell[2]).epsilon;
                double temp_surface = temperature;
                double rhob_surface = arena(i_cell[0],i_cell[1],i_cell[2]).rhob;
                double p_surface = arena(i_cell[0],i_cell[1],i_cell[2]).p;
                std::array<double, 4> u_surface = arena(i_cell[0],i_cell[1],i_cell[2]).u;
                
                ofs_freezeout
                << coord->tau*hbarc << " " //tau in [fm]
                << x_surface << " " //x in [fm]
                << y_surface << " " //y in [fm]
                << eta_surface << " " //eta in [fm]
                << dsigma[0] << " " //dsigma_tau in [GeV^-3]
                << 0.0 << " " //dsigma_x in [GeV^-3]
                << 0.0 << " " //dsigma_y in [GeV^-3]
                << 0.0 << " " //dsigma_eta in [GeV^-3]
                << u_surface[0] << " "// u0
                << u_surface[1] << " "// u1
                << u_surface[2] << " "// u2
                << u_surface[3] << " "// u3
                << e_surface << " "// energy density in [GeV^4]
                << temp_surface << " "// temparature in [GeV]
                << 0.0 << " " // for muB
                << (e_surface+p_surface)/temp_surface << " " // (e+p)/T
                << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " // for Wmunu
                << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "// for Wmunu
                << 0.0 << " " // for bulk viscosity
                << rhob_surface << " " // for rhoB
                << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0
                << "\n"; // for qmu
                
            }//y
        }//x
        
    }//eta
    
    ofs_freezeout << std::flush;
    
}


int Freezeout::FullFreezeout(){
    JSINFO
    << "<-[PPM] Finding Isothermal Freezeout Surface->";
    
    
    double temp_max = 0.0;
    std::array<double, 4> dsigma = GetDsigma( 1 );
    
    for (int ieta = 3; ieta < grid_neta-3; ieta++) {
        double eta = coord->GetEta(ieta);
        if( fabs(eta) > rapidity_window ){continue;}
        for (int ix = 3; ix < grid_nx-3; ix++) {
            for (int iy = 3; iy < grid_ny-3; iy++) {
                
                std::array<int, 3> i_cell = { ix, iy, ieta };
                
                std::array<double, 3> x_cell =
                { coord->GetX(ix)*hbarc, coord->GetY(iy)*hbarc, eta };
                
                BulkFreezeout( i_cell, x_cell, dsigma);
                SurfaceFreezeout( i_cell, x_cell, dsigma);
                
                if( temp_max < arena(i_cell[0],i_cell[1],i_cell[2]).T ){
                    temp_max = arena(i_cell[0],i_cell[1],i_cell[2]).T;
                }
                
            }//y
        }//x
        
    }//eta
    
    ofs_freezeout << std::flush;
    
    if( temp_max >= temp_fo ){
        JSINFO<< "<-[PPM] Maximum Temperature: " << temp_max << "GeV ->";
        return ppm_running;
    }else{
        JSINFO
        << "<-[PPM] Completion of Freezeout. tau = "<<coord->tau*hbarc<<" fm ->";
        return ppm_finished;
    }
    
}

std::array<double, 4> Freezeout::GetDsigma( int time_shift ){
    
    
    double tau_surface = coord->tau;
    
    if( time_shift == 1 ){
        tau_surface -= 0.5 * coord->dtau;
    }
    
    std::array<double, 4> dsigma =
    {                  coord->dx[0] * coord->dx[1] * (coord->tau *coord->dx[2]),
        -coord->dtau                * coord->dx[1] * (tau_surface*coord->dx[2]),
        -coord->dtau * coord->dx[0]                * (tau_surface*coord->dx[2]),
        -coord->dtau * coord->dx[0] * coord->dx[1]
    };
    return dsigma; // in GeV^{-3}
}

void Freezeout::BulkFreezeout(const std::array<int, 3> &i_cell,
                              const std::array<double, 3> &x_cell,
                              const std::array<double, 4> &dsigma ){
    
    double temp_prev = arena(i_cell[0],i_cell[1],i_cell[2]).T_prev;
    double temp_curr = arena(i_cell[0],i_cell[1],i_cell[2]).T;
    
    int sign;
    bool bulk_freezed_out = false;
    if( temp_prev > temp_fo  && temp_curr <= temp_fo) { //cooling
        sign = 1;
        bulk_freezed_out = true;
    }else if( temp_prev <= temp_fo  && temp_curr > temp_fo ){//heating
        sign = -1;
        bulk_freezed_out = true;
    }
    
    if( bulk_freezed_out ){
        
        double   x_surface = x_cell[0];
        double   y_surface = x_cell[1];
        double eta_surface = x_cell[2];
        double e_surface = arena(i_cell[0],i_cell[1],i_cell[2]).epsilon;
        double temp_surface = temp_curr;
        double rhob_surface = arena(i_cell[0],i_cell[1],i_cell[2]).rhob;
        double p_surface = arena(i_cell[0],i_cell[1],i_cell[2]).p;
        std::array<double, 4> u_surface = arena(i_cell[0],i_cell[1],i_cell[2]).u;

        ofs_freezeout
        << coord->tau*hbarc << " " //tau in [fm]
        << x_surface << " " //x in [fm]
        << y_surface << " " //y in [fm]
        << eta_surface << " " //eta
        << double(sign)*dsigma[0] << " " //dsigma_tau in [GeV^-3]
        << 0.0 << " " //dsigma_x in [GeV^-3]
        << 0.0 << " " //dsigma_y in [GeV^-3]
        << 0.0 << " " //dsigma_eta in [GeV^-3]
        << u_surface[0] << " "// u0
        << u_surface[1] << " "// u1
        << u_surface[2] << " "// u2
        << u_surface[3] << " "// u3
        << e_surface << " "// energy density in [GeV^4]
        << temp_surface << " "// temparature in [GeV]
        << 0.0 << " " // for muB
        << (e_surface+p_surface)/temp_surface << " " // (e+p)/T
        << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " // for Wmunu
        << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "// for Wmunu
        << 0.0 << " " // for bulk viscosity
        << rhob_surface << " " // for rhoB
        << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0  // for qmu
        << "\n";


    }
    
}

void Freezeout::SurfaceFreezeout(const std::array<int, 3> &i_cell,
                                 const std::array<double, 3> &x_cell,
                                 const std::array<double, 4> &dsigma ){
    
    //directions
    
    for(int d=0; d<3; d++){
        
        const int   ix_next = i_cell[0] + point_next[d][0];
        const int   iy_next = i_cell[1] + point_next[d][1];
        const int ieta_next = i_cell[2] + point_next[d][2];
        
        double temp_curr = arena(i_cell[0],i_cell[1],i_cell[2]).T;
        double temp_next = arena(ix_next,iy_next,ieta_next).T;
        
        int sign;
        bool surface_freezed_out = false;
        if( temp_curr > temp_fo  && temp_next <= temp_fo) { //cooling
            sign = 1;
            surface_freezed_out = true;
        }else if( temp_curr <= temp_fo  && temp_next > temp_fo ){//heating
            sign = -1;
            surface_freezed_out = true;
        }
        
        if( surface_freezed_out ){
            
            coord->CountTau( -0.5 );
            
            std::array<double, 5> U_surf = arena(i_cell[0],i_cell[1],i_cell[2]).U_surf[d];
            double   x_surface = 0.5*( x_cell[0] + coord->GetX(ix_next)*hbarc );
            double   y_surface = 0.5*( x_cell[1] + coord->GetY(iy_next)*hbarc );
            double eta_surface = 0.5*( x_cell[2] + coord->GetEta(ieta_next) );
            
            std::array<double, 4> u_surface;
            double e_surface, temp_surface, rhob_surface, p_surface;

            fval->GetThermalVal( U_surf, u_surface, e_surface, rhob_surface, p_surface, temp_surface);

            //Geometrical Source Term Effect
            if(d == 2){
                U_surf[0] -= 0.5 * (coord->dtau/coord->tau)*U_surf[0]*u_surface[3]*u_surface[3]/u_surface[0]/u_surface[0];
                U_surf[0] -= 0.5 * coord->dtau*p_surface*(1.0-(u_surface[3]/u_surface[0])*(u_surface[3]/u_surface[0]));
                U_surf[3] -= 0.5 * (coord->dtau/coord->tau)*U_surf[3];
                fval->GetThermalVal( U_surf, u_surface, e_surface, rhob_surface, p_surface, temp_surface);
            }
            
            

            ofs_freezeout
            << coord->tau*hbarc << " " //tau in [fm]
            << x_surface << " " //x in [fm]
            << y_surface << " " //y in [fm]
            << eta_surface << " " //eta in [fm]
            << 0.0 << " " //dsigma_tau in [GeV^-3]
            << double(sign * point_next[d][0])*dsigma[1] << " " //dsigma_x in [GeV^-3]
            << double(sign * point_next[d][1])*dsigma[2] << " " //dsigma_y in [GeV^-3]
            << double(sign * point_next[d][2])*dsigma[3] << " " //dsigma_eta in [GeV^-3]
            << u_surface[0] << " "// u0
            << u_surface[1] << " "// u1
            << u_surface[2] << " "// u2
            << u_surface[3] << " "// u3
            << e_surface << " "// energy density in [GeV^4]
            << temp_surface << " "// temparature in [GeV]
            << 0.0 << " " // for muB
            << (e_surface+p_surface)/temp_surface << " " // (e+p)/T
            << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " // for Wmunu
            << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " "// for Wmunu
            << 0.0 << " " // for bulk viscosity
            << rhob_surface << " " // for rhoB
            << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 // for qmu
            << "\n";
            
            coord->CountTau( 0.5 );
            
        }
    }
}
