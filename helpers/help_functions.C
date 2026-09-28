#include <cmath>

#include "TH1F.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"

int extractClusterId(const std::string& filepath){
	
	size_t pos = filepath.find("cluster");
	if (pos != std::string::npos){

		std::string clusterStr = filepath.substr(pos + 7, 2);
		try {
			return std::stoi(clusterStr);
		} catch (const std::exception& e){
			std::cerr << "Warning: failed to convert cluster number" << clusterStr <<std::endl;
	 	return -1;
		}
	}
	return -1;
}

 
TGraphErrors* MakeLogGraph(TH1F *h1) {
    if (!h1) return nullptr;

    TGraphErrors *graph = new TGraphErrors();
    
    TString name = TString::Format("g_ln_%s", h1->GetName());
    TString title = TString::Format(";%s;ph. e.",  h1->GetXaxis()->GetTitle());
    
    graph->SetName(name);
    graph->SetTitle(title);
 
    int point_index = 0;
 
    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double content = h1->GetBinContent(i);

        double x       = h1->GetBinCenter(i);
        double error_x = h1->GetBinWidth(i) / 2.0;

        double new_y, new_error_y;
        if (content <= 0) {
            new_y       = 0;
            new_error_y = 0;
        } else {
            double error_y = h1->GetBinError(i);
            new_y       = -std::log(content);
            new_error_y = error_y / content;
        }

        graph->SetPoint(point_index, x, new_y);
        graph->SetPointError(point_index, error_x, new_error_y);

        point_index++;
    }
 
    return graph;
}


bool EnsureDirectoryExists(const TString& dir_path) {
    // .Data() extracts const char* from TString
    if (gSystem->AccessPathName(dir_path.Data())) {
        std::cout << "[ROOT] Directory '" << dir_path << "' not found. Creating..." << std::endl;
        
        int status = gSystem->MakeDirectory(dir_path.Data());
        
        if (status == 0) {
            std::cout << "[ROOT] Directory successfully created." << std::endl;
            return true;
        } else {
            std::cerr << "[ERROR] Failed to create directory: " << dir_path << std::endl;
            return false;
        }
    }
    
    std::cout << "[ROOT] Directory '" << dir_path << "' already exists." << std::endl;
    return true;
}




