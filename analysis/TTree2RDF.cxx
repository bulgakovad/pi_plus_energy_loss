// to run, use:
//g++ TTree2RDF.cxx -o executable `root-config --cflags --libs`
// ./executable

#include <iostream>
#include "plots.cxx"
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <fstream>
#include <TLine.h>
#include <TLegend.h>
#include <string>
#include <TMath.h>
#include <cmath>
#include <chrono>
#include "TSystem.h"

#include <ROOT/RVec.hxx>



int isData = 0;  // 1 for real data, 0 for MC
bool isBigStatistics = true;
bool toFarm = false;

std::string farm_out = (toFarm == true) ? "/farm_out/" : "/";
std::string test_stat = (isBigStatistics == true) ? "" : "_TEST";


// ROOT file

//std::string root_file_path = "../data/clasdis_rga_fall18_inbendingBIG.root";

std::string root_file_path = (isBigStatistics == true) ? "../data/multipion_1trig_el_no_other_el.root" : "../data/pi_plus_toy.root";


//OUTPUT folder

//const std::string OUTPUT_FOLDER = "analysis_out_clasdis_rga_fall18_inbending_BIG" + farm_out ;
const std::string OUTPUT_FOLDER = "multipion_analysis_out_pi_plus_toy" + test_stat + farm_out ;



ROOT::RDataFrame convert_ttrees_to_rdataframe(const std::string &root_file_path) {
    TFile *file = TFile::Open(root_file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file " << root_file_path << std::endl;
        return ROOT::RDataFrame(0);
    }

    std::vector<std::string> keys;
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (std::string(key->GetClassName()) == "TTree") {
            keys.push_back(key->GetName());
        }
    }

    if (keys.empty()) {
        std::cerr << "No TTrees found in the ROOT file." << std::endl;
        return ROOT::RDataFrame(0);
    }

    std::string tree_name = keys[0];
    std::cout << "Processing TTree: " << tree_name << std::endl;

    ROOT::RDataFrame rdf(tree_name, root_file_path);
    file->Close();
    return rdf;
}

// Use ROOT::RDF::RNode instead of RDataFrame& to fix type mismatch

int main() {

    ROOT::EnableImplicitMT(); // Enable multi-threading


    auto start = std::chrono::high_resolution_clock::now(); // STRAT

    // Load ROOT file and convert TTrees to RDataFrame
    auto rdf = convert_ttrees_to_rdataframe(root_file_path);
    if (rdf.GetColumnNames().empty()) {
        std::cerr << "Error: Could not create RDataFrame." << std::endl;
        return 1;
    }
    //Create folder if it does not exist
    std::cout << "Output folder name = " << OUTPUT_FOLDER<<std::endl;
    gSystem->mkdir(OUTPUT_FOLDER.c_str(), true);
    std::cout << "ROOT file path = " << root_file_path << std::endl;


    // Define necessary variables in RDataFrame
    auto init_rdf = rdf//.Filter(p_piplus_gen > 0.0)
                        //.Filter("n_rec_piplus == 1")
                        .Define("delta_p", "p_piplus_rec - p_piplus_gen")
                        .Define("piplus_rec_4_momentum", "TLorentzVector(px_piplus_rec, py_piplus_rec, pz_piplus_rec, 0.13957039)") // const number is mass of pi+
                        .Define("piplus_gen_4_momentum", "TLorentzVector(px_piplus_gen, py_piplus_gen, pz_piplus_gen, 0.13957039)") // const number is mass of pi+
                        .Define("Phi_rec", "piplus_rec_4_momentum.Phi()*TMath::RadToDeg()")
                        .Define("Phi_gen", "piplus_gen_4_momentum.Phi()*TMath::RadToDeg()")
                        .Define("Theta_rec", "piplus_rec_4_momentum.Theta()*TMath::RadToDeg()")
                        .Define("Theta_gen", "piplus_gen_4_momentum.Theta()*TMath::RadToDeg()")
                        .Define("Theta_piplus_DC", "TMath::ATan(sqrt(x1_piplus*x1_piplus + y1_piplus*y1_piplus)/z1_piplus)*TMath::RadToDeg()")
                        .Define("Phi_piplus_DC", "TMath::ATan2(y1_piplus, x1_piplus)*TMath::RadToDeg()")
                        .Define("detector", [](int status) {
                            return status < 4000 ? std::string("FD") : (status < 8000 ? std::string("CD") : std::string("NA"));
                        }, {"status_piplus"})
                        .Define("electron_rec_4_momentum", "TLorentzVector(px_electron_rec, py_electron_rec, pz_electron_rec, 0.0)")
                        .Define("electron_gen_4_momentum", "TLorentzVector(px_electron_gen, py_electron_gen, pz_electron_gen, 0.0)")
                        .Define("Phi_electron_rec", "electron_rec_4_momentum.Phi()*TMath::RadToDeg()")
                        .Define("Phi_electron_gen", "electron_gen_4_momentum.Phi()*TMath::RadToDeg()")
                        .Define("Theta_electron_rec", "electron_rec_4_momentum.Theta()*TMath::RadToDeg()")
                        .Define("Theta_electron_gen", "electron_gen_4_momentum.Theta()*TMath::RadToDeg()")
                        .Define("E_piplus_rec", "sqrt(px_piplus_rec*px_piplus_rec + py_piplus_rec*py_piplus_rec + pz_piplus_rec*pz_piplus_rec + 0.13957039*0.13957039)") // const number is mass of pi+
                        .Define("E_piplus_gen", "sqrt(px_piplus_gen*px_piplus_gen + py_piplus_gen*py_piplus_gen + pz_piplus_gen*pz_piplus_gen + 0.13957039*0.13957039)") // const number is mass of pi+
                        .Define("delta_E", "E_piplus_rec - E_piplus_gen")
                        .Define("dp_norm", "delta_p / p_piplus_rec")
                        .Define("sector_pi_plus", "sector_piplus") 
                        .Define("DC_fiducial_cut_electron", "detector == \"FD\" && edge1_electron > 5.0 && edge2_electron > 5.0 && edge3_electron > 10.0")
                        .Define("DC_fiducial_cut_piplus",  "detector == \"FD\" && edge1_piplus > 10.0 && edge2_piplus > 10.0 && edge3_piplus > 10.0")
                        .Define("Vz_cut", "vz_piplus > -10 && vz_piplus < 2")
                        .Define("Vxy_cut", "abs(vx_piplus) < 2.0 && abs(vy_piplus) < 2.0");





    // Print column names 
   //std::cout << "Columns in RDataFrame:" << std::endl;
   //for (const auto &col : init_rdf.GetColumnNames()) { 
   //std::cout << col << std::endl; }
   
    //init_rdf.Filter("n_rec_piplus > 1 ").Display("status_piplus_all",50)->Print();
    //int n_electrons = init_rdf.Filter("n_rec_piplus == 1 && detector == \"FD\" && Vz_cut == true && Vxy_cut == true").Count().GetValue();
    
    //Theta_VS_momentum_FD_CD(init_rdf, OUTPUT_FOLDER);
    //Phi_VS_momentum_FD_CD(init_rdf, OUTPUT_FOLDER);
    //Phi_VS_Theta_FD_CD(init_rdf, OUTPUT_FOLDER);
    //delta_P_VS_P_rec_FD_CD(init_rdf, OUTPUT_FOLDER);
    //plot_delta_P_VS_P_rec_FD(init_rdf, OUTPUT_FOLDER);
    //plot_delta_P(init_rdf, OUTPUT_FOLDER);
    //plot_momenta_components(init_rdf, OUTPUT_FOLDER);

    //Theta_VS_momentum_FD_gen_cuts(init_rdf, OUTPUT_FOLDER, 0, 10);
    //Theta_VS_momentum_FD_gen_cuts(init_rdf, OUTPUT_FOLDER, 50,55,4,5);
    //Theta_VS_momentum_FD_gen_cuts(init_rdf, OUTPUT_FOLDER, 0,30,4,5);



    
    ThetaToPBins cfg;

    // for FD
    //cfg[{5,10}]  = {3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0};
    //cfg[{10,15}] = {1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8};
    //cfg[{15,20}] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.2, 5.4, 5.6,5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4};    
    cfg[{20,21}] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0};
    //cfg[{25,30}] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8};
    //cfg[{30,35}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0,7.2, 7.4};
    //cfg[{35,36}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2};
    //cfg[{36,37}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2};
    //cfg[{37,38}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8};
    //cfg[{38,39}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8};
    //cfg[{39,40}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8};
    //cfg[{40,41}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2};
    //cfg[{41,45}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8}; 
    //cfg[{45,60}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8};

    //cfg[{38.0,38.5}] = {1.28, 1.32 };
   

    //cfg[{30,31}] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0};



    //cfg[{0,180}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

    // for CD

    //cfg[{30,35}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{38,39}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{40,45}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{45,50}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{50,55}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{55,60}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{60,65}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{65,70}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{70,75}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{75,80}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6};
    //cfg[{80,85}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6 };
    //cfg[{85,90}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6 };
    //cfg[{90,95}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6 };
    //cfg[{95,100}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2 };
    //cfg[{100,105}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0 };
    //cfg[{105,110}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6 };
    //cfg[{110,115}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2 };
    //cfg[{115,120}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 };
    //cfg[{120,125}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6 };
    //cfg[{125,130}] = {0.2,0.4, 0.6, 0.8, 1.0, 1.2};
    //cfg[{130,135}] = {0.2,0.4, 0.6, 0.8, 1.0 };
    //cfg[{135,140}] = {0.2,0.4, 0.6 };




    //delta_P_VS_P_rec_unified_1D(init_rdf, OUTPUT_FOLDER, cfg, true, "FD","Theta_rec");
    //plot_momentum_SingleRecPion_inside_theta_dp_bin(init_rdf, OUTPUT_FOLDER, 38.0, 39.0, -0.13, -0.05, 1.8, 2.0);

    //Theta_VS_momentum_FD_CD(init_rdf, OUTPUT_FOLDER);

    //plot_theta_vs_vz_pion(init_rdf, OUTPUT_FOLDER);
    //plot_Vx_VS_Vy_piplus(init_rdf, OUTPUT_FOLDER);


    //delta_P_VS_P_rec_FD_CD(init_rdf, OUTPUT_FOLDER);
    //Theta_VS_momentum_FD_CD(init_rdf, OUTPUT_FOLDER);



    //Theta_VS_Phi_per_P_bin(init_rdf, OUTPUT_FOLDER, "p_piplus_rec","Theta_rec","Phi_rec");
    //P_VS_Phi_per_Theta_bin(init_rdf, OUTPUT_FOLDER, "p_piplus_rec","Theta_rec","Phi_rec");
    //P_VS_Theta_per_Phi_bin(init_rdf, OUTPUT_FOLDER, "p_piplus_rec","Theta_rec","Phi_rec");

    //deltaP_VS_Prec_per_Theta_bin(init_rdf, OUTPUT_FOLDER, "FD", true, "p_piplus_rec","delta_p","Theta_rec");

    //deltaP_VS_Prec_per_Phi_bin(init_rdf, OUTPUT_FOLDER, "FD", true, "p_piplus_rec","delta_p","Phi_piplus_DC");

    //deltaP_per_Phi_bin(init_rdf, OUTPUT_FOLDER, "FD", 38.0, 39.0);

   //deltaP_VS_Phi(init_rdf, OUTPUT_FOLDER, "FD", 37.0, 40.0);

    //plot_X_vs_Y_piplus_FD(init_rdf, OUTPUT_FOLDER);

    //plot_deltaP_multiRecPions_inside_theta_momentum_bin(init_rdf, OUTPUT_FOLDER,38.0, 39.0, 1.0, 1.2);

    plot_deltaP_SingleRecPion_inside_theta_momentum_bin(init_rdf, OUTPUT_FOLDER, 38.0, 39.0, 0.0, 10000);

 


   

    auto end = std::chrono::high_resolution_clock::now(); // END

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time of execution: " << elapsed.count() << " sec" << std::endl;

    return 0;
}
