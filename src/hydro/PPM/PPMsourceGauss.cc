#include "PPMsourceGauss.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "tinyxml2.h"
#include "FluidDynamics.h"


using namespace Jetscape;
using namespace std;

SourceGauss::SourceGauss( int run_num_in,
                         std::shared_ptr<Coordinates> coord_in,
                         std::shared_ptr<FluidValuables> fval_in,
                         InitData &DATA_in,
                         SCGrid &arena_in
                         ):DATA(DATA_in), arena(arena_in){
    
    run_num = run_num_in;
    coord = coord_in;
    fval = fval_in;
    
    tinyxml2::XMLElement *hd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro" );
    tinyxml2::XMLElement *ppm=hd->FirstChildElement("PPM");
    tinyxml2::XMLElement *pliq=ppm->FirstChildElement("PPMsourceGauss");
    
    string s = pliq->FirstChildElement( "name" )->GetText();
    JSINFO <<"<-[PPM] "<< s << " to be initilizied ... ->";
    
    if (pliq) {

        
        std::string source_input = pliq->FirstChildElement("sourceInput")->GetText();
        pliq->FirstChildElement("tau_thermal")->QueryDoubleText(&tau_th);
        pliq->FirstChildElement("sigma_trans")->QueryDoubleText(&sigma_t);
        pliq->FirstChildElement("sigma_long")->QueryDoubleText(&sigma_l);
        
        
        std::string input_file_name = GenerateFilename(source_input);
        
        p_ads = std::unique_ptr<PartonAdS>(new PartonAdS( input_file_name ));

        partonE = 0.0, partonPx = 0.0, partonPy = 0.0, partonPz = 0.0;
        sourceE = 0.0, sourcePx = 0.0, sourcePy = 0.0, sourcePz = 0.0;

        sigma_cut = 5.75;
        dr_sub = 0.03;
        deta_sub = 0.03;
        subgrid_disc = 0.5*DATA.delta_x/dr_sub;

        r_max = sigma_cut * sigma_t;
        n_r = int(r_max/dr_sub);
        dr = r_max/n_r;

        eta_max = sigma_cut * sigma_l;
        n_eta = 2.0*int(eta_max/deta_sub);
        deta = eta_max/n_eta;

        correction = Correction();
        
        JSINFO << "<-[PPM] "<< s << " is initilizied ... ->";

    } else {
        JSWARN << " : PPM SourceGauss not properly initialized in XML file ...";
        exit(-1);
    }

}

SourceGauss::~SourceGauss(){
    JSINFO << "<-[PPM] Deleting SourceGauss ->";
}// destructor

std::string SourceGauss::GenerateFilename(std::string source_input){
    return source_input+"/source_run_"+to_string(run_num)+".txt";
}

double SourceGauss::Correction(){
    double total_weight = 0.0;
    double r = - 0.5*dr;
        for( int ir = 0; ir < n_r; ir++ ){
            r+=dr;
            double gauss_t = exp(-r*r/2.0/sigma_t/sigma_t)/2.0/M_PI/sigma_t/sigma_t;
            int n_phi = 8;
            int n_phi_sub = 8*int(r/dr/subgrid_disc);
            if( n_phi < n_phi_sub ) n_phi = n_phi_sub;
            double dphi = 2.0*PI/n_phi;
            double weight_t = r * dr * dphi * gauss_t;
            for( int iphi = 0; iphi < n_phi; iphi++ ){
                double phi = (double(iphi))*dphi;
                double eta = - 0.5*(n_eta+1)*deta;
                for( int ieta = 0; ieta < n_eta; ieta++ ){
                    eta+=deta;
                    double gauss_l = exp(-eta*eta/2.0/sigma_l/sigma_l)/sqrt(2.0*M_PI)/sigma_l;
                    double weight_l = deta * gauss_l;
                    double weight_full = weight_t*weight_l;
                    total_weight += weight_full;
            }//eta
        }//phi
    }//r
    JSINFO << "<-[PPM] source correction: " << total_weight << " ->";
    return total_weight;
}




void SourceGauss::AddSourceAtAPoint(PartonAdS::Droplet drop,
                                    double weight,
                                    double r, double phi, double delta_eta){

    double x = r*cos(phi) + drop.x;
    double y = r*sin(phi) + drop.y;
    double eta_s = delta_eta + drop.eta_s;

    int ix = coord->GetIx(x/hbarc);
    int iy = coord->GetIy(y/hbarc);
    int ieta = coord->GetIeta(eta_s);
//
////    JSINFO << "drop_x: " << drop.x;
////    JSINFO << "drop_y: " << drop.y;
////    JSINFO << "drop_eta_s: " << drop.eta_s;
////    JSINFO << "x: " << x << " " << evolve_ptr->GetX(ix)*hbarc;
////    JSINFO << "y: " << y << " " << evolve_ptr->GetY(iy)*hbarc;
////    JSINFO << "eta_s: " << eta_s << " " << evolve_ptr->GetEta(ieta);

    double dp_tau =   drop.de * cosh(eta_s) - drop.dpz * sinh(eta_s);
    double dp_eta = - drop.de * sinh(eta_s) + drop.dpz * cosh(eta_s);

    double du0dv = weight*dp_tau;
    double du1dv = weight*drop.dpx;
    double du2dv = weight*drop.dpy;
    double du3dv = weight*dp_eta;


    arena(ix,iy,ieta).U[0] += du0dv/coord->dV;
    arena(ix,iy,ieta).U[1] += du1dv/coord->dV;
    arena(ix,iy,ieta).U[2] += du2dv/coord->dV;
    arena(ix,iy,ieta).U[3] += du3dv/coord->dV;


    if( arena(ix,iy,ieta).U[0] < 0.0){
        JSINFO << "<-[PPM] Source Error";
        JSINFO << "U0=" << arena(ix,iy,ieta).U[0];
        JSINFO
        << "U0*dV ="
        << arena(ix,iy,ieta).U[0]*coord->dV;
        JSINFO << "dp_tau ="
        << dp_tau << " GeV ->";
        exit(-1);
    }
    
    std::array<int, 3> index = {ix, iy, ieta};
    
    coord->CountTau( 0.5 );
    fval->SetThermalVal( index );
    coord->CountTau( -0.5 );

    sourceE  += weight*drop.de;
    sourcePx += du1dv;
    sourcePy += du2dv;
    sourcePz += weight*drop.dpz;
}



void SourceGauss::AddSourceTauEta(){
    
    JSINFO << "<-[PPM] Gaussian Source ->";
    double tau = coord->tau * hbarc;
    double dtau = coord->dtau * hbarc;
    JSINFO
    << "<-[PPM] Adding External Source (Tau-Eta, Sub Grid) tau = "
    << tau- dtau/2.0 << "-"
    << tau + dtau/2.0 <<" fm ->";

    int i_drop = 0;

    for (auto &drop: p_ads->droplets) {

        double tau_source = drop.tau + tau_th;


        if( tau_source >= tau - 0.5*dtau &&
            tau_source <  tau + 0.5*dtau ){
            double sign = 1.0;
            if(drop.status <= -1){ //negative parton
                sign = -1.0;
            }

        partonE  += sign * drop.de;
        partonPx += sign * drop.dpx;
        partonPy += sign * drop.dpy;
        partonPz += sign * drop.dpz;
        i_drop++;


            double r = - 0.5*dr;

            for( int ir = 0; ir < n_r; ir++ ){

                r+=dr;

                double gauss_t = exp(-r*r/2.0/sigma_t/sigma_t)/2.0/M_PI/sigma_t/sigma_t;

                int n_phi = 8;
                int n_phi_sub = 8*int(r/dr/subgrid_disc);
                if( n_phi < n_phi_sub ) n_phi = n_phi_sub;
                double dphi = 2.0*PI/n_phi;

                double weight_t = r * dr * dphi * gauss_t;

                for( int iphi = 0; iphi < n_phi; iphi++ ){
                    double phi = (double(iphi))*dphi;



                    double eta = - 0.5*(n_eta+1)*deta;
                    for( int ieta = 0; ieta < n_eta; ieta++ ){

                        eta+=deta;

                        double gauss_l = exp(-eta*eta/2.0/sigma_l/sigma_l)/sqrt(2.0*M_PI)/sigma_l;

                        double weight_l = deta * gauss_l;

                        double weight_full = weight_t*weight_l/correction;
                        //tiny cheating by numerical correction

                        AddSourceAtAPoint(drop,
                                          sign * weight_full,
                                          r, phi, eta);

                    }//eta
                }//phi
            }//r
        }

    }

    JSINFO << "<-[PPM] Lost Parton Momenutm (total), "
    << "E = " << partonE << " GeV, "
    << "Px = " << partonPx << " GeV/c, "
    << "Py = " << partonPy << " GeV/c, "
    << "Pz = " << partonPz << " GeV/c. ->";

    JSINFO << "<-[PPM]         Added Source (total), "
    << "E = " << sourceE << " GeV, "
    << "Px = " << sourcePx << " GeV/c, "
    << "Py = " << sourcePy << " GeV/c, "
    << "Pz = " << sourcePz << " GeV/c. ->";
    
}

