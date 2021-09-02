#include <iostream>
#include <string>
#include <memory>
#include <cmath>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

#include "cerang.h"

int main(void) {

	//
	// 1. Get MC data from .root file and prepare it
	//

	// Open the root file with histograms
	TFile dataFile("data-proton-1PeV.root", "read");

	// Get the histogram for s = 1 (atmospheric slice number 12)
	// Scale the histogram and set its graphical properties
	auto thetaHistPtr = std::unique_ptr<TH1>( dataFile.Get<TH1>("hTheta_12") );
	
	thetaHistPtr->Scale(1.0/thetaHistPtr->Integral("width"));
	
	thetaHistPtr->SetLineWidth(2);
	thetaHistPtr->SetLineColor(kBlack);
	
	// Get values of refractive index and shower age for slice 12
	auto refractiveIndexHistPtr = std::unique_ptr<TH1>( dataFile.Get<TH1>("refractiveIndex") );
	auto showerAgeHistPtr = std::unique_ptr<TH1>( dataFile.Get<TH1>("showerAge") );
	
	double refractiveIndex = refractiveIndexHistPtr->GetBinContent(12);
	double showerAge = showerAgeHistPtr->GetBinContent(12);
	
	// Shower energy is 1 PeV, and we write it in units of TeV
	double showerEnergyTeV = 1e3;
	
	
	
	//
	// 2. Compute the normalized distribution
	//
	
	// Number of points in the TGraph
	int nPoints = 3000;
	
	// Theta range and step size according to MC histogram and nPoints
	double minTheta = thetaHistPtr->GetBinLowEdge(1);
	double maxTheta = thetaHistPtr->GetBinLowEdge(thetaHistPtr->GetNbinsX()+1);
	double dTheta = (maxTheta - minTheta) / double(nPoints);
	
	// Create a TGraph to hold the parametrized function values
	TGraph functionGraph;
	functionGraph.SetLineWidth(2);
	functionGraph.SetLineColorAlpha(kRed,0.7);
	
	// Loop between min/max theta to fill the TGraph
	for (double theta = 0.5*dTheta; theta < maxTheta; theta += dTheta) {
		// Compute the parametrized function for the current value of theta
		// Use the values of shower age, refractive index, and shower energy defined above
		double functionValue = CherenkovPDF(theta, showerAge, refractiveIndex, showerEnergyTeV);
		                       ////////////
		
		// Fill the TGraph for the parametrized function
	 	functionGraph.SetPoint(functionGraph.GetN(), theta, functionValue);
	 	
	 	// Note there are other two functions that may be useful:
	 	// CherenkovCDF: computes the cumulative distribution
	 	// CherenkovIntegral: computes the integral of the PDF between given values of theta_min and theta_max
	}
	
	
	
	//
	// 3. Plot everything in a pdf file
	//
	
	// Create a canvas and split it into two pads
	TCanvas canvas("canvas","canvas", 1200, 600);
	
	canvas.Divide(2,1);
	
	// Left canvas: region from 0º to 5º
	canvas.cd(1);
	
	TH2D hAxisLeft("", "", 100, 0, 5*std::acos(-1)/180., 100, 0, 24.99);
	hAxisLeft.SetStats(false);
	hAxisLeft.GetXaxis()->SetTitle("Angle between shower axis and photon direction [rad]");
	hAxisLeft.GetYaxis()->SetTitle("dN_{#gamma}/d#theta");
	
	hAxisLeft.Draw("axis");
	thetaHistPtr->Draw("hist same");
	functionGraph.Draw("l same");
	
	// Right pad: region from 5º to 60º (y axis in log scale)
	canvas.cd(2);
	gPad->SetLogy(true);
	
	TH2D hAxisRight("", "", 100, 5*std::acos(-1)/180., 60*std::acos(-1)/180., 100, 3.e-3, 5);
	hAxisRight.SetStats(false);
	hAxisRight.GetXaxis()->SetTitle("Angle between shower axis and photon direction [rad]");
	hAxisRight.GetYaxis()->SetTitle("dN_{#gamma}/d#theta");
	hAxisRight.GetYaxis()->SetTitleOffset(1.0);
	
	hAxisRight.Draw("axis");
	thetaHistPtr->Draw("hist same");
	functionGraph.Draw("l same");
	
	// Put a legend in the right pad
	TLegend leg(0.55,0.8,0.89,0.89);
	leg.SetBorderSize(0);
	leg.AddEntry(thetaHistPtr.get(), "CORSIKA", "l");
	leg.AddEntry(&functionGraph, "Parametrization", "l");
	leg.Draw();
	
	// Print the canvas
	canvas.Print("demo.pdf");
	

	return 0;	
}
