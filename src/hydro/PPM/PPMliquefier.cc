#include "PPMliquefier.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "tinyxml2.h"
#include "FluidDynamics.h"
#include "PPMevolve.h"

using namespace Jetscape;
using namespace std;

Liquefier::Liquefier( int run_num_in,
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
    tinyxml2::XMLElement *pliq=ppm->FirstChildElement("PPMLiquefier");

    string s = pliq->FirstChildElement( "name" )->GetText();
    JSINFO <<"<-[PPM] "<< s << " to be initilizied ... ->";
    
    if (pliq) {
        source_input_file = pliq->FirstChildElement("sourceInput")->GetText();
        pliq->FirstChildElement("d_diff")->QueryDoubleText(&d_diff);
        pliq->FirstChildElement("t_damp")->QueryDoubleText(&t_damp);
        pliq->FirstChildElement("tau_relax")->QueryDoubleText(&tau_relax);
        
        p_cloud = std::unique_ptr<PartonCloud>(new PartonCloud( source_input_file ));

        dr_sub = 0.04;
        subgrid_disc = 0.5*DATA.delta_x/dr_sub;

        c_diff = sqrt(d_diff/tau_relax);
        gamma_relax = 0.5/tau_relax;

        if(c_diff > 1.0){
            JSWARN << "<-[PPM] Liquefier: Causality is broken with this set of parameters ...";
            exit(-1);
        }

        partonE = 0.0, partonPx = 0.0, partonPy = 0.0, partonPz = 0.0;
        sourceE = 0.0, sourcePx = 0.0, sourcePy = 0.0, sourcePz = 0.0;

        JSINFO << "<-[PPM] "  << s << " is initilizied ... ->";

    } else {
        JSWARN << " : PPM Liquefier not properly initialized in XML file ...";
        exit(-1);
    }
//    //ConservationCheck();
}

Liquefier::~Liquefier(){
    JSINFO << "<-[PPM] Deleting Liquefier ->";
}// destructor

void Liquefier::AddSourceCart(){

    JSINFO << "<-[PPM] Liquefier ->";
    
    double t = coord->tau * hbarc;
    double dt = coord->dtau * hbarc;
    
    JSINFO
    << "<-[PPM] Adding External Source (Cartesian, Sub Grid) t = "
    << t - dt/2.0 << "-"
    << t + dt/2.0 <<" fm ->";

    int idrop = 0;
    double source_sign = 1.0;
    
    JSINFO
    << "<-[PPM] test = " << p_cloud->drops.size();


    while(idrop < p_cloud->drops.size()){

        double delta_t = t - p_cloud->drops[idrop].x[0];

        if( p_cloud->drops[idrop].pstat == -1 ){
            source_sign = - 1.0;
        }else{
            source_sign = 1.0;
        }

        int d_status = 1;

        if( t_damp >= delta_t - 0.5*dt &&
            t_damp < delta_t + 0.5*dt){

            d_status = 0;
            total_weight = 0.0;

            double r_max = c_diff * delta_t;
            int n_r = int(r_max/dr_sub);
            double dr = r_max/n_r;

            int n_theta, n_phi;
            double dtheta, dphi;

            double r = - 0.5*dr;
            for( int ir = 0; ir < n_r; ir++ ){
                r+=dr;

                n_theta = 4;
                int n_theta_sub = 4*int(r/dr/subgrid_disc);
                if( n_theta < n_theta_sub ) n_theta = n_theta_sub;

                dtheta = PI/n_theta;

                double weight_r
                = CausalDiffusionGreenFunctionSmooth( delta_t, r);

                for( int itheta = 0; itheta < n_theta; itheta++ ){
                    double theta = (double(itheta) + 0.5)*dtheta;

                    n_phi = 8;
                    int n_phi_sub = 8*int(r*sin(theta)/dr/subgrid_disc);
                    if( n_phi < n_phi_sub ) n_phi = n_phi_sub;
                    dphi = 2.0*PI/n_phi;

                    double weight
                    = weight_r*(r*r*sin(theta))*dr*dtheta*dphi;

                    for( int iphi = 0; iphi < n_phi; iphi++ ){
                        double phi = (double(iphi))*dphi;

                        AddSourceAtAPoint(p_cloud->drops[idrop],
                                          source_sign * weight,
                                          r, theta, phi);

                    }//phi
                }//theta
            }//r

            r = r_max;

            double weight_delta
            = CausalDiffusionGreenFunctionDeltaValue( delta_t, r );

            for( int itheta = 0; itheta < n_theta; itheta++ ){
                double theta = (double(itheta) + 0.5)*dtheta;

                n_phi = 8;
                int n_phi_sub = 8*int(r*sin(theta)/dr/subgrid_disc);
                if( n_phi < n_phi_sub ) n_phi = n_phi_sub;

                dphi = 2.0*PI/n_phi;

                double weight
                = weight_delta*(r*r*sin(theta))*dtheta*dphi;

                for( int iphi = 0; iphi < n_phi; iphi++ ){
                    double phi = (double(iphi))*dphi;

                    AddSourceAtAPoint(p_cloud->drops[idrop],
                                      source_sign * weight,
                                      r, theta, phi);

                }//phi
            }//theta
        }

        if(d_status == 0){
            partonE  += source_sign * p_cloud->drops[idrop].p[0];
            partonPx += source_sign * p_cloud->drops[idrop].p[1];
            partonPy += source_sign * p_cloud->drops[idrop].p[2];
            partonPz += source_sign * p_cloud->drops[idrop].p[3];
            p_cloud->drops.erase(p_cloud->drops.begin() + idrop);
            //JSINFO << "<-[PPM] total_weight: " << total_weight << " ->";
        }else{
            idrop++;
        }

    }

    JSINFO << "<-[PPM] Lost  Parton (total), "
    << "E = " << partonE << " GeV, "
    << "Px = " << partonPx << " GeV/c, "
    << "Py = " << partonPy << " GeV/c, "
    << "Pz = " << partonPz << " GeV/c. ->";

    JSINFO << "<-[PPM] Added Source (total), "
    << "E = " << sourceE << " GeV, "
    << "Px = " << sourcePx << " GeV/c, "
    << "Py = " << sourcePy << " GeV/c, "
    << "Pz = " << sourcePz << " GeV/c. ->";

    JSINFO
    << "<-[PPM] Remaining Sources: "
    << p_cloud->drops.size()
    << " ->";

}

void Liquefier::AddSourceAtAPoint(PartonCloud::Droplet drop,
                                  double weight,
                                  double r, double theta, double phi){

    total_weight += weight;
    
    double x
    = r*sin(theta)*cos(phi) + drop.x[1];;
    double y
    = r*sin(theta)*sin(phi) + drop.x[2];
    double z
    = r*cos(theta) + drop.x[3];
    
    int ix = coord->GetIx(x/hbarc);
    int iy = coord->GetIy(y/hbarc);
    int ieta = coord->GetIeta(z/hbarc);
    
    std::array<double, 4> sourceXdt = {0,0,0,0};
    for( int d4=0; d4<4; d4++){
        sourceXdt[d4]
        = drop.p[d4]
        * weight;
        arena(ix,iy,ieta).U[d4] +=  sourceXdt[d4]/coord->dV;
    }

    if( arena(ix,iy,ieta).U[0] < 0){
        JSINFO << "<-[PPM] Source Error";
        JSINFO << "U0=" << arena(ix,iy,ieta).U[0];
        JSINFO
        << "U0*dV ="
        << arena(ix,iy,ieta).U[0]*coord->dV;
        JSINFO << "p0 ="
        << drop.p[0] << " GeV ->";
        exit(-1);
    }
    
    std::array<int, 3> index = {ix, iy, ieta};
    
    coord->CountTau( 0.5 );
    fval->SetThermalVal( index );
    coord->CountTau( -0.5 );
    
    sourceE  += sourceXdt[0];
    sourcePx += sourceXdt[1];
    sourcePy += sourceXdt[2];
    sourcePz += sourceXdt[3];
}

double Liquefier::CausalDiffusionGreenFunctionSmooth(double delta_t, double delta_r){
    
    if( delta_r < (c_diff * delta_t) ){
        
        double u = sqrt( c_diff*c_diff * delta_t*delta_t - delta_r*delta_r );
        double x = gamma_relax * u / c_diff; // unitless

        double i1 = gsl_sf_bessel_I1(x);
        double i2 = gsl_sf_bessel_In(2,x);
        
        return
        ( (exp(-gamma_relax*delta_t))/(20.*PI) )
        *(2.*gamma_relax*gamma_relax/c_diff)
        *(i1/(c_diff*u) + 4.*delta_t *i2/u/u);
        
    }else{
        return 0.0;
    }
    
}


double Liquefier::CausalDiffusionGreenFunctionDeltaValue(double delta_t, double delta_r){
    
    return
    ( (exp(-gamma_relax*delta_t))/(20.*PI) )
    *(8. - 3.*exp(-gamma_relax*delta_t)
      + 2.*gamma_relax*delta_t
      +4.*gamma_relax*gamma_relax*delta_t*delta_t )
    /delta_r/delta_r;
    
}


////- function for debug
//void Liquefier::ConservationCheck(){
//    
//    JSINFO << "Checking the conservation in Causal Diffusin Green Function";
//    ofstream ofs;
//    ofs.open( "test.txt" );
//
//    double dt = 0.3;
//    double dr = 0.1;
//    double t = 0.0;
//    for(int it=0; t < 40; it++){
//        t = it*dt;
//
//        double integrated_value = 0.0;
//        double integrated_smooth = 0.0;
//        double integrated_spike = 0.0;
//        for(double r = 0.5 * dr; r < 3.*t; r+=dr){
//            double all = 4. * PI * r * r * dr * CausalDiffusionGreenFunction( t, r);
//            double smooth = 4. * PI * r * r * dr * CausalDiffusionGreenFunctionSmooth( t, r );
//            double spike = 4. * PI * r * r * dr * CausalDiffusionGreenFunctionDeltaFinite( t, r );
//            integrated_value += all;
//            integrated_smooth += smooth;
//            integrated_spike += spike;
//        }
//
//        ofs
//        << setprecision(10) << t << " "
//        << setprecision(10) << integrated_value <<" "
//        << setprecision(10) << integrated_smooth <<" "
//        << setprecision(10) << integrated_spike << endl;
//    }
//    
//    ofs << flush;
//    ofs.close();
// 
//    exit(1);
//
//}
