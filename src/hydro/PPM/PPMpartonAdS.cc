
#include <math.h>

#include "PPMpartonAdS.h"
#include "JetScapeLogger.h"


using namespace Jetscape;
using namespace std;

PartonAdS::PartonAdS(string source_input_file){
    
    JSINFO << "<-[PPM] Creating PartonAdS ->";
    
    totalE = 0.0;
    totalPx=0.0;
    totalPy=0.0;
    totalPz=0.0;

    LoadPartonList(source_input_file);
}

PartonAdS::~PartonAdS(){
    Message();
    droplets.clear();
    JSINFO << "<-[PPM] Deleting PartonAdS ->";
}// destructor

void PartonAdS::LoadPartonList(string source_input_file){

    JSINFO
    << "<-[PPM] Loading Source Term Input File: "
    << source_input_file
    << " ->";
    
    ifstream ifs;
    stringstream str_stream;

    ifs.open(source_input_file.c_str()); //open the input file
    str_stream << ifs.rdbuf(); //read the file

    if (ifs.is_open()){
        ifs.close();
    }else{
        ifs.close();
        JSINFO << "<-[PPM] SourceGauss, " << source_input_file << " - Not Found ->";
        return;
    }

    std::string line;
    
    int sn, pid, status;
    double tau, x, y, eta_s; // in [fm]
    double e, px, py, pz; // in [GeV]
    double de, dpx, dpy, dpz; // in [GeV]
    
    
    while ( getline( str_stream, line ) ){
        if( line.find("#") == std::string::npos &&
            ! line.empty() ){
           
            //JSINFO << "Input: " << line;
            sscanf(line.data(),
                    "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &sn, &pid, &status,
                    &tau, &x, &y, &eta_s,
                    &e, &px, &py, &pz,
                    &de, &dpx, &dpy, &dpz );
            
            Droplet this_drop;
            
            this_drop.sn = sn;
            this_drop.pid = pid;
            this_drop.status = status;
            
            this_drop.tau = tau;
            this_drop.x = x;
            this_drop.y = y;
            this_drop.eta_s = eta_s;
            
            this_drop.e = e;
            this_drop.px = px;
            this_drop.py = py;
            this_drop.pz = pz;
            
            this_drop.de = de;
            this_drop.dpx = dpx;
            this_drop.dpy = dpy;
            this_drop.dpz = dpz;
            
            droplets.push_back(this_drop);
            
            if( status != -1){
                totalE += de;
                totalPx += dpx;
                totalPy += dpy;
                totalPz += dpz;
            }else{
                totalE -= de;
                totalPx -= dpx;
                totalPy -= dpy;
                totalPz -= dpz;
            }
            
        }
    }

    Message();

}

void PartonAdS::Message(){
    
    JSINFO
    << "<-[PPM] Source Term Loaded, Total Source Number = "
    << droplets.size()
    << " ->";

    JSINFO
    << "<-[PPM] Source to be deposited into the fluid: ->";
    JSINFO
    << "<-[PPM] total_e: " << totalE
    << " GeV, total_px: " << totalPx
    << " GeV, total_py: " << totalPy
    << " GeV, total_pz: " << totalPz
    << " GeV, total_mass^2: " << (totalE*totalE - totalPx*totalPx - totalPy*totalPy - totalPz*totalPz)
    << " GeV^2 ->";
    
}
