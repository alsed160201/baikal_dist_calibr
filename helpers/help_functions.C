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
        // -ln is undefined for content <= 0: skip the bin instead of adding a fake point
        if (content <= 0) continue;

        double x       = h1->GetBinCenter(i);
        double error_x = h1->GetBinWidth(i) / 2.0;

        double new_y       = -std::log(content);
        double new_error_y = h1->GetBinError(i) / content;

        graph->SetPoint(point_index, x, new_y);
        graph->SetPointError(point_index, error_x, new_error_y);

        point_index++;
    }
 
    return graph;
}


// Fraction of zero OMs per bin, 1 - sig/total, with binomial errors sqrt(f*(1-f)/T).
// Divided via sig (filled directly, errors sqrt(S)) rather than the zero-OM histogram
// (built with Add, so its errors are sqrt(T+S), not sqrt(Z)); the error of 1 - f equals
// the error of f. Bins with total = 0 are set to 0 so MakeLogGraph skips them.
TH1F* MakeZeroFraction(TH1F *sig, TH1F *total, const char* name) {
    TH1F *frac = (TH1F*) sig->Clone(name);
    frac->SetTitle("fraction of zero OM vs dist from track");
    frac->Divide(sig, total, 1.0, 1.0, "B");

    for (int i = 1; i <= frac->GetNbinsX(); ++i) {
        if (total->GetBinContent(i) <= 0) {
            frac->SetBinContent(i, 0);
            frac->SetBinError(i, 0);
            continue;
        }
        frac->SetBinContent(i, 1.0 - frac->GetBinContent(i));
    }

    return frac;
}


// Zero = Total - Sig built with Add gets sqrt(T+S) errors; replace them with Poisson sqrt(Z).
// The histogram is already normalized by nEvents, so content = Z/N and
// error = sqrt(Z)/N = sqrt(content/N). Valid only while all event weights are 1.
void SetPoissonErrors(TH1F *h, double nEvents) {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double c = h->GetBinContent(i);
        h->SetBinError(i, c > 0 ? std::sqrt(c / nEvents) : 0);
    }
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




