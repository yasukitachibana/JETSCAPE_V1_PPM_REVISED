#include "PPMprinter.h"
#include "JetScapeLogger.h"

using namespace Jetscape;
using namespace std;

Printer::Printer( int run_num_in,
                 std::shared_ptr<EOS> eos_in,
                 std::shared_ptr<Coordinates> coord_in,
                 InitData &DATA_in,
                 SCGrid &arena_in
                 ):DATA(DATA_in), arena(arena_in){
    
    JSINFO << "<-[PPM] Creating Printer ->";
    
    run_num = run_num_in;
    coord = coord_in;
    
    if(DATA.write_output == 1){
        
        if( DATA.profileType == 2 || DATA.profileType == 3 ){
            hbc3 = hbarc;
        }
        
        grid_nx = DATA.nx;
        grid_ny = DATA.ny;
        grid_neta = DATA.neta;
        
        rapidity_window = DATA.rapidity_window;
        transverse_square = DATA.transverse_square;
        
        profile_filename = GenerateProfileFilename();
        SurvayConfiguration();
        
    }else{
        JSINFO << "<-[PPM] Printer: OFF ->";
    }
    
}

Printer::~Printer(){
}// destructor

std::string Printer::GenerateProfileFilename(){
    return DATA.profile_output + "_run_" + to_string(run_num) + ".txt";
}

void Printer::PrintProfile( int ppm_status ){
    
    if(DATA.write_output == 1){
        
        JSINFO << "<-[PPM] Printing Profile ->";
        FILE *out_file;
        
        std::string open_option;
        if( ppm_status == ppm_not_start ){
            open_option = "w";
        }else{
            open_option = "a";
        }
        open_option += "b";
        
        out_file = std::fopen(profile_filename.c_str(), open_option.c_str());
        
        double t = coord->tau * hbarc;
        
        for (int ix = 3; ix < grid_nx-3; ix++) {
            double x = coord->GetX(ix)*hbarc;
            if( fabs(x) > transverse_square ){continue;}
            for (int iy = 3; iy < grid_ny-3; iy++) {
                double y = coord->GetY(iy)*hbarc;
                if( fabs(y) > transverse_square ){continue;}
                for (int ieta = 3; ieta < grid_neta-3; ieta++) {
                    double eta = coord->GetEta(ieta)*hbc3;
                    if( fabs(eta) > rapidity_window ){continue;}
                    
                    float cell_infos[]
                    ={
                        static_cast<float> (t),//0
                        static_cast<float> (x),//1
                        static_cast<float> (y),//2
                        static_cast<float> (eta),//3
                        static_cast<float> (arena(ix,iy,ieta).u[0]),//4
                        static_cast<float> (arena(ix,iy,ieta).u[1]),//5
                        static_cast<float> (arena(ix,iy,ieta).u[2]),//6
                        static_cast<float> (arena(ix,iy,ieta).u[3]),//7
                        static_cast<float> (arena(ix,iy,ieta).T),//8
                        static_cast<float> (arena(ix,iy,ieta).epsilon),//9
                        static_cast<float> (arena(ix,iy,ieta).p),//10
                        static_cast<float> (coord->dV),//11
                        static_cast<float> (arena(ix,iy,ieta).U[0]),//12
                        static_cast<float> (arena(ix,iy,ieta).U[1]),//13
                        static_cast<float> (arena(ix,iy,ieta).U[2]),//14
                        static_cast<float> (arena(ix,iy,ieta).U[3]),//15
                    };
                    
                    //binary file
                    fwrite(cell_infos, sizeof(float), 16, out_file);
                    
                }//eta
            }//y
        }//x
        
        fclose(out_file);
        
    }
    
}


void Printer::SurvayConfiguration(){
    
    int n_x_in = 0;
    double x_min = coord->GetX(3)*hbarc;
    double x_max = coord->GetX(grid_nx-4)*hbarc;
    
    for (int ix = 3; ix < grid_nx-3; ix++) {
        double x = coord->GetX(ix)*hbarc;
        if( fabs(x) > transverse_square ){continue;}
        if( n_x_in == 0 ){ x_min = x; }
        x_max = x;
        n_x_in++;
    }
    
    int n_y_in = 0;
    double y_min = coord->GetY(3)*hbarc;
    double y_max = coord->GetY(grid_ny-4)*hbarc;
    
    for (int iy = 3; iy < grid_ny-3; iy++) {
        double y = coord->GetY(iy)*hbarc;
        if( fabs(y) > transverse_square ){continue;}
        if( n_y_in == 0 ){ y_min = y; }
        y_max = y;
        n_y_in++;
    }
    
    int n_eta_in = 0;
    double eta_min = coord->GetEta(3)*hbc3;
    double eta_max = coord->GetEta(grid_neta-4)*hbc3;
    
    for (int ieta = 3; ieta < grid_neta-3; ieta++) {
        double eta = coord->GetEta(ieta)*hbc3;
        if( fabs(eta) > rapidity_window ){continue;}
        if( n_eta_in == 0 ){ eta_min = eta; }
        eta_max = eta;
        n_eta_in++;
    }
    
    JSINFO << "<-[PPM] Printer is Ready.->";
    JSINFO << "<-[PPM] Profile File Name: '" << profile_filename << "' ->";
    JSINFO << "<-[PPM] Printed area: ->";
    JSINFO
    << "<-[PPM] x: "<<x_min<<" - "<<x_max
    << " fm, y: "<<y_min<<" - "<<y_max
    << " fm, eta(z): "<<eta_min<<" - "<<eta_max << " (fm) ->";
    JSINFO
    << "<-[PPM] nx x ny x neta(nz) = "<<n_x_in<<" x "<<n_y_in<<" x "<<n_eta_in<<" ->";
    
}
