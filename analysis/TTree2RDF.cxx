// to run, use:
//g++ TTree2RDF.cxx -o executable `root-config --cflags --glibs`
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



int isData = 0;  // 1 for real data, 0 for MC
bool isBigStatistics = false;
bool toFarm = false;

std::string farm_out = (toFarm == true) ? "/farm_out/" : "/";


// ROOT file

//std::string root_file_path = "../data/clasdis_rga_fall18_inbendingBIG.root";
std::string root_file_path = "../data/pi_plus_toy.root";

//OUTPUT folder

//const std::string OUTPUT_FOLDER = "analysis_out_clasdis_rga_fall18_inbending_BIG" + farm_out ;
const std::string OUTPUT_FOLDER = "analysis_out_pi_plus_toy" + farm_out ;



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
    auto start = std::chrono::high_resolution_clock::now(); // STRAT

    // Load ROOT file and convert TTrees to RDataFrame
    ROOT::EnableImplicitMT(); // Enable multi-threading
    auto rdf = convert_ttrees_to_rdataframe(root_file_path);
    if (rdf.GetColumnNames().empty()) {
        std::cerr << "Error: Could not create RDataFrame." << std::endl;
        return 1;
    }
    //Create folder if it does not exist
    std::cout << "Output folder name = " << OUTPUT_FOLDER<<std::endl;
    gSystem->mkdir(OUTPUT_FOLDER.c_str(), true);

    // Define necessary variables in RDataFrame
    auto init_rdf = rdf//.Filter(p_piplus_gen > 0.0)
                        //.Filter(p_piplus_rec > 0.0)
                        .Define("delta_p", "p_piplus_rec - p_piplus_gen")
                        .Define("piplus_rec_4_momentum", "TLorentzVector(px_piplus_rec, py_piplus_rec, pz_piplus_rec, 0.13957039)") // const number is mass of pi+
                        .Define("piplus_gen_4_momentum", "TLorentzVector(px_piplus_gen, py_piplus_gen, pz_piplus_gen, 0.13957039)") // const number is mass of pi+
                        .Define("Phi_rec", "piplus_rec_4_momentum.Phi()*TMath::RadToDeg()")
                        .Define("Phi_gen", "piplus_gen_4_momentum.Phi()*TMath::RadToDeg()")
                        .Define("Theta_rec", "piplus_rec_4_momentum.Theta()*TMath::RadToDeg()")
                        .Define("Theta_gen", "piplus_gen_4_momentum.Theta()*TMath::RadToDeg()")
                        //.Define("Theta_piplus_DC", "TMath::ATan(sqrt(x1_piplus*x1_piplus + y1_piplus*y1_piplus)/z1_piplus)*TMath::RadToDeg()")
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
                        .Define("DC_fiducial_cut_electron", "detector == \"FD\" && edge1_electron > 5.0 && edge2_electron > 5.0 && edge3_electron > 10.0")
                        .Define("DC_fiducial_cut_piplus",  "detector == \"FD\" && edge1_piplus > 2.5 && edge2_piplus > 2.5 && edge3_piplus > 9.0");





    // Print column names 
    //std::cout << "Columns in RDataFrame:" << std::endl;
    //for (const auto &col : init_rdf.GetColumnNames()) { 
    //    std::cout << col << std::endl; 
    //}

    Theta_VS_momentum_FD_CD(init_rdf, OUTPUT_FOLDER);
    Phi_VS_momentum_FD_CD(init_rdf, OUTPUT_FOLDER);
    Phi_VS_Theta_FD_CD(init_rdf, OUTPUT_FOLDER);
    delta_P_VS_P_rec_FD_CD(init_rdf, OUTPUT_FOLDER);
    plot_delta_P_VS_P_rec_FD(init_rdf, OUTPUT_FOLDER);
    plot_delta_P(init_rdf, OUTPUT_FOLDER);
    plot_momenta_components(init_rdf, OUTPUT_FOLDER);


    
    ThetaToPBins cfg;


    cfg[{5.0, 10.0}] = {4.0, 4.5, 5.0};
    cfg[{10.0, 15.0}] = {1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    cfg[{15.0, 20.0}] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    cfg[{20.0, 25.0}] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    cfg[{25.0, 30.0}] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    cfg[{30.0, 35.0}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    cfg[{35.0, 40.0}] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};



    delta_P_VS_P_rec_FD_unified_1D(init_rdf, OUTPUT_FOLDER, cfg, false);



    auto end = std::chrono::high_resolution_clock::now(); // END

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time of execution: " << elapsed.count() << " sec" << std::endl;

    return 0;
}
