#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>
#include <array>
#include <vector>

#include <TCanvas.h>
#include <TSystem.h>
#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TObject.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TGraph.h>
#include <TTree.h>

#include "cerang.h"

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

template <class ... T >
bool MakeTest(bool(*testFunction)(T ...), std::string name, T ... args)
{
  std::cout << BOLDWHITE << "Starting " << name << " test" << RESET << std::endl;

  bool result = testFunction(args ...);

  if (result == true) {
    std::cout << BOLDGREEN << "Passed ";
  } else {
    std::cout << BOLDRED << "Failed ";
  }

  std::cout << name << " test." << RESET << std::endl;

  return result; 
}

bool DoTestParameters() {

  // gSystem->Unlink("tests/parameters");
  // gSystem->mkdir("tests/parameters", true);

  double minShowerAge = 0.000001;
  double maxShowerAge = 3.0;

  double minRefIndex = 1.0000001;
  double maxRefIndex = 1.0003;

  double minETeV = 1.0e3;
  double maxETeV = 1.0e8;

  int stepsShowerAge = 10;
  int stepsRefIndex = 10;
  int stepsEnergy = 5;

  double deltaAge = (maxShowerAge - minShowerAge) / double(stepsShowerAge);
  double deltaRef = (maxRefIndex - minRefIndex) / double(stepsRefIndex);
  double factorEng = std::pow(maxETeV/minETeV, 1/double(stepsEnergy));

  for (int iAge = 0; iAge <= stepsShowerAge; iAge++) {
    double showerAge = minShowerAge + iAge * deltaAge;

    for (double iRef = 0; iRef <= stepsRefIndex; iRef++) {
      double refractiveIndex = minRefIndex + iRef * deltaRef;

      for (double iEng = 0; iEng <= stepsEnergy; iEng++) {
        double showerEnergy = minETeV * std::pow(factorEng, iEng);

        auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergy);

        std::cout
          << std::setw(13) << showerAge
          << std::setw(13) << refractiveIndex-1
          << std::setw(13) << showerEnergy
          << std::setw(5) << ""
          << std::setw(13) << parameters.at(0)
          << std::setw(13) << parameters.at(1)
          << std::setw(13) << parameters.at(2)
          << std::setw(13) << parameters.at(3)
        << std::endl;
      }
    }
  }

  return true;
}

bool DoTestSmallAngleIntegral()
{
  double showerAge = 1.0;
  double refractiveIndex = 1.0003;
  double showerEnergyTeV = 1e18/1e12;

  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV);

  double thetaEm = std::acos(1.0/refractiveIndex);

  double lowAngle = 0;

  double minHighAngle = lowAngle;
  double maxHighAngle = thetaEm;
  int stepsHighAngle = 1000;

  double deltaHighAngle = (maxHighAngle - minHighAngle) / double(stepsHighAngle);

  double integralSimpson = 0;
  double integral3per8 = 0;
  double integralTrapezoid = 0;

  for (int iHigh = 1; iHigh <= stepsHighAngle; iHigh++) {

    double highAngle = minHighAngle + iHigh * deltaHighAngle;

    double a = highAngle - deltaHighAngle;
    double b = highAngle;
    double fa = CherenkovUnnormalizedDensity(a, thetaEm, parameters);
    double fb = CherenkovUnnormalizedDensity(b, thetaEm, parameters);
    double fm = CherenkovUnnormalizedDensity((a+b)*0.5, thetaEm, parameters);
    double fma = CherenkovUnnormalizedDensity((2*a+b)/3.0, thetaEm, parameters);
    double fmb = CherenkovUnnormalizedDensity((a+2*b)/3.0, thetaEm, parameters);
    integralSimpson += (fa + 4*fm + fb) * (b-a)/6.0;
    integral3per8 += (fa + 3*fma + 3*fmb + fb) * (b-a)/8.0;
    integralTrapezoid += 0.5 * (fa + fb) * (b-a);

    double integralSeries = SmallAngleIntegral(lowAngle, highAngle, thetaEm, parameters, 10);

    std::cout
      << std::setw(3) << iHigh
      // << std::setw(15) << lowAngle
      << std::setw(15) << highAngle
      // << std::setw(15) << integralSeries
      // << std::setw(15) << integralSimpson
      // << std::setw(15) << integral3per8
      << std::setw(15) << integralSimpson/integral3per8 - 1
      << std::setw(15) << integral3per8/integralSeries - 1
      << std::setw(15) << integralTrapezoid/integralSeries - 1
    << std::endl;
  }

  return true;
}

bool DoTestLargeAngleIntegral()
{
  double showerAge = 1.0;
  double refractiveIndex = 1.0003;
  double showerEnergyTeV = 1e18/1e12;

  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV);

  double thetaEm = std::acos(1.0/refractiveIndex);

  double lowAngle = 1.1*thetaEm;

  double minHighAngle = lowAngle;
  double maxHighAngle = std::asin(1);
  int stepsHighAngle = 10000;

  double deltaHighAngle = (maxHighAngle - minHighAngle) / double(stepsHighAngle);

  double integralSimpson = 0;

  for (int iHigh = 1; iHigh <= stepsHighAngle; iHigh++) {

    double highAngle = minHighAngle + iHigh * deltaHighAngle;

    double a = highAngle - deltaHighAngle;
    double b = highAngle;
    double fa = CherenkovUnnormalizedDensity(a, thetaEm, parameters);
    double fb = CherenkovUnnormalizedDensity(b, thetaEm, parameters);
    double fm = CherenkovUnnormalizedDensity((a+b)*0.5, thetaEm, parameters);
    integralSimpson += (fa + 4*fm + fb) * (b-a)/6.0;

    double integral = LargeAngleIntegral(lowAngle, highAngle, thetaEm, parameters, 1e-3);

    std::cout
      << std::setw(3) << iHigh
      // << std::setw(15) << lowAngle
      << std::setw(15) << highAngle
      << std::setw(15) << integral
      << std::setw(15) << integralSimpson
      << std::setw(15) << integral/integralSimpson - 1
    << std::endl;
  }

  return true;
}


bool DoTestAgainstMC(std::vector<std::string> fileList = {})
{
  // Loop over input files
  for (const auto & mcFilePath : fileList) {

    // Gest base name from the current root file
    std::string basename = gSystem->BaseName(mcFilePath.c_str());

    // Check if this is, indeed, a root file
    if (basename.substr(basename.size() - 5) != ".root") {
      // If not, the test has failed due to invalid arguments
      std::cout << "Received a non .root file" << std::endl;
      std::cout << "File path is: " << mcFilePath << std::endl;
      return false;
    } else {
      // If so, remove the .root extension from the basename
      basename = basename.erase(basename.size() - 5);
    }

    // Create an empty directory to store showers from the current file
    gSystem->Unlink(("tests/mc/" + basename).c_str());
    gSystem->mkdir(("tests/mc/" + basename).c_str(), true);

    // Open file with the angular distributions of Cherenkov photons
    TFile mcFile(mcFilePath.c_str(), "read");

    // Check if file was correctly opened
    if (mcFile.IsZombie()) {
      std::cout << "Unable to open file " << mcFilePath << std::endl;
      return false;
    }

    // Print name of the current file
    std::cout << "Parsing file: " << mcFilePath << std::endl;

    // Read the event tree and set its branches to arrays of floats
    std::unique_ptr<TTree> eventTree( mcFile.Get<TTree>("tEvt") );

    std::array<float, 273> eventHeaderData;
    std::array<float, 273> eventEndData;

    eventTree->SetBranchAddress("fEvtHeader", eventHeaderData.data());
    eventTree->SetBranchAddress("fEvtEnd", eventEndData.data());

    // Get number of showers in this file
    int numberOfShowers = eventTree->GetEntries();

    // Loop over showers
    for (int iShower = 1; iShower <= numberOfShowers; iShower++) {

      // Read data from the event tree
      eventTree->GetEntry(iShower - 1);

      int runNumber = eventHeaderData.at(43);
      int eventNumber = eventHeaderData.at(1);

      std::cout
        << "\r" << std::setw(50) << "" << "\r"
        << "+ Reading shower " << iShower << " of " << numberOfShowers 
      << std::flush;

      double showerEnergyTeV = eventHeaderData.at(3) * 1e-3;

      // Get a pointer to the directory containing information from the current shower
      std::string eventDirName = "Event_" + std::to_string(runNumber) + "_" + std::to_string(eventNumber);

      std::unique_ptr<TDirectory> eventDirectory( mcFile.Get<TDirectory>(eventDirName.c_str()) );

      if (!eventDirectory) {
        continue;
      }
      
      // Read the profiles of age and refractive index for the current event
      std::unique_ptr<TH1> ageProfile( eventDirectory->Get<TH1>("Profiles/hAge") );
      std::unique_ptr<TH1> refProfile( eventDirectory->Get<TH1>("Profiles/hReffrac") );

      // This is the number of atmospheric bins in the current event
      int nBins = ageProfile->GetNbinsX();

      // Create a canvas to print the plots from this shower
      std::string canvasPath = "tests/mc/" + basename + "/" + eventDirName + ".pdf";

      TCanvas canvas("", "", 1200, 600);
      canvas.Divide(2,1);

      // Open the ouput file as an empty one
      canvas.Print((canvasPath + "[").c_str());

      // Loop over bins of atmospheric depth
      for (int iBin = 1; iBin <= nBins; iBin++) {

        //
        // Retrieve the histogram from the input file
        //

        // The name of the current histogram
        std::string histogramName = "hTheta_" + std::to_string(iBin);

        // Get the histogram with the angular distribution
        std::unique_ptr<TH1> thetaHist( eventDirectory->Get<TH1>(histogramName.c_str()) );

        // Skip empty histograms
        if (thetaHist->Integral() == 0) {
          continue;
        }

        // Scale the histogram
        thetaHist->Scale(1.0/thetaHist->Integral("width"));

        // Get number of bins
        int numberOfThetaBins = thetaHist->GetNbinsX();

        // Get shower age and refractive index
        double showerAge = ageProfile->GetBinContent(iBin);
        double refractiveIndex = refProfile->GetBinContent(iBin);
        double thetaEm = std::acos(1.0/refractiveIndex);

        // Check shower age and refractive index
        if (!(0 < showerAge && showerAge < 3)) {
          continue;
        } else if (!(1.0 < refractiveIndex && refractiveIndex < 1.0003)) {
          continue;
        }

        //
        // Create a smooth function as a TGraph
        //
        const static double pi = std::acos(-1);
        
        double minTheta = 0;
        double midTheta = pi/36.; // 5 deg.
        double maxTheta = pi/3.; // 3 deg.

        int stepsInTheta = 10000;

        double deltaThetaLeft = (midTheta - minTheta) / double(stepsInTheta);
        double deltaThetaRight = (maxTheta - midTheta) / double(stepsInTheta);

        TGraph graphLeft;
        TGraph graphRight;

        auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV);
        
        double normalization = CherenkovUnnormalizedIntegral(minTheta, maxTheta, thetaEm, parameters);

        for (int iTheta = 0; iTheta < stepsInTheta; iTheta++) {
          double thetaLeft = minTheta + (iTheta + 0.5) * deltaThetaLeft;
          double thetaRight = midTheta + (iTheta + 0.5) * deltaThetaRight;

          double functionLeft = CherenkovUnnormalizedDensity(thetaLeft, thetaEm, parameters) / normalization;
          double functionRight = CherenkovUnnormalizedDensity(thetaRight, thetaEm, parameters) / normalization;

          if (!std::isfinite(functionLeft) || functionLeft < 0) {
            std::cout << "The parametrized density is not finite!" << std::endl;
            std::cout << "+ function = " << functionLeft << std::endl;
            std::cout << "+ showerAge = " << showerAge << std::endl;
            std::cout << "+ refractiveIndex = " << refractiveIndex << std::endl;
            std::cout << "+ showerEnergyTeV = " << showerEnergyTeV << std::endl;
            std::cout << "+ theta = " << thetaLeft << std::endl;
            std::cout << "+ thetaEm = " << thetaEm << std::endl;
            std::cout << "+ nu = " << parameters[0] << std::endl;
            std::cout << "+ t1 = " << parameters[1] << std::endl;
            std::cout << "+ t2 = " << parameters[2] << std::endl;
            std::cout << "+ ep = " << parameters[3] << std::endl;
            std::cout << "+ normalization = " << normalization << std::endl;
            std::cout << std::endl;
            std::cout << "---> File: " << mcFilePath << std::endl;
            std::cout 
              << "---> Run/Event/iBin/iTheta: "
              << runNumber
              << "/" 
              << eventNumber 
              << "/"
              << iBin
              << "/"
              << iTheta
            << std::endl;
            return false;
          } else if (!std::isfinite(functionRight) || functionRight < 0) {
            std::cout << "The parametrized density is not finite!" << std::endl;
            std::cout << "+ function = " << functionRight << std::endl;
            std::cout << "+ showerAge = " << showerAge << std::endl;
            std::cout << "+ refractiveIndex = " << refractiveIndex << std::endl;
            std::cout << "+ showerEnergyTeV = " << showerEnergyTeV << std::endl;
            std::cout << "+ theta = " << thetaRight << std::endl;
            std::cout << "+ thetaEm = " << thetaEm << std::endl;
            std::cout << "+ nu = " << parameters[0] << std::endl;
            std::cout << "+ t1 = " << parameters[1] << std::endl;
            std::cout << "+ t2 = " << parameters[2] << std::endl;
            std::cout << "+ ep = " << parameters[3] << std::endl;
            std::cout << "+ normalization = " << normalization << std::endl;
            std::cout << std::endl;
            std::cout << "---> File: " << mcFilePath << std::endl;
            std::cout 
              << "---> Run/Event/iBin/iTheta: "
              << runNumber
              << "/" 
              << eventNumber 
              << "/"
              << iBin
              << "/"
              << iTheta
            << std::endl;
            return false;
          }

          graphLeft.SetPoint(iTheta, thetaLeft, functionLeft);
          graphRight.SetPoint(iTheta, thetaRight, functionRight);
        }

        //
        // Create a histogram from the parametrized function
        //
        TH1D parametrizedHistogram("", "", numberOfThetaBins, minTheta, maxTheta);

        for (int iTheta = 1; iTheta <= numberOfThetaBins; iTheta++) {
          double thetaLeft = parametrizedHistogram.GetBinLowEdge(iTheta);
          double thetaRight = parametrizedHistogram.GetBinLowEdge(iTheta+1);

          double content = CherenkovUnnormalizedIntegral(thetaLeft, thetaRight, thetaEm, parameters);
          content /= normalization;
          content /= thetaRight - thetaLeft;

          if (!std::isfinite(content) || content < 0) {
            std::cout << "The parametrized integral is not finite!" << std::endl;
            std::cout << "+ normalized integral = " << content << std::endl;
            std::cout << "+ showerAge = " << showerAge << std::endl;
            std::cout << "+ refractiveIndex = " << refractiveIndex << std::endl;
            std::cout << "+ showerEnergyTeV = " << showerEnergyTeV << std::endl;
            std::cout << "+ thetaLeft = " << thetaLeft << std::endl;
            std::cout << "+ thetaRight = " << thetaRight << std::endl;
            std::cout << "+ thetaEm = " << thetaEm << std::endl;
            std::cout << "+ nu = " << parameters[0] << std::endl;
            std::cout << "+ t1 = " << parameters[1] << std::endl;
            std::cout << "+ t2 = " << parameters[2] << std::endl;
            std::cout << "+ ep = " << parameters[3] << std::endl;
            std::cout << "+ normalization = " << normalization << std::endl;
            std::cout << std::endl;
            std::cout << "---> File: " << mcFilePath << std::endl;
            std::cout 
              << "---> Run/Event/iBin/iTheta: "
              << runNumber
              << "/" 
              << eventNumber 
              << "/"
              << iBin
              << "/"
              << iTheta
            << std::endl;
            return false;
          }

          parametrizedHistogram.SetBinContent(iBin, content);
        }

        //
        // Set graphical properties of the objects to be drawn
        //
        thetaHist->SetStats(false);
        thetaHist->SetLineColor(kBlack);
        thetaHist->SetLineWidth(2);

        graphLeft.SetLineWidth(2);
        graphLeft.SetLineColor(kRed);

        graphRight.SetLineWidth(2);
        graphRight.SetLineColor(kRed);

        parametrizedHistogram.SetLineColor(kBlue);
        parametrizedHistogram.SetLineWidth(2);

        //
        // The left plot (small angular region)
        //
        canvas.cd(1);

        TH2D hAxisLeft("hAxisLeft", "", 100, minTheta, midTheta, 100, 0, 1.05*thetaHist->GetMaximum());
        hAxisLeft.SetStats(false);
        hAxisLeft.Draw("axis");

        thetaHist->Draw("hist same");

        graphLeft.Draw("l same");

        parametrizedHistogram.Draw("hist same");

        //
        // The right plot (large angular region)
        //
        canvas.cd(2);
        gPad->SetLogy(true);

        TH2D hAxisRight("hAxisRight", "", 100, midTheta, maxTheta, 100, 0.5*thetaHist->GetMinimum(0), thetaHist->GetMaximum());
        hAxisRight.SetStats(false);
        hAxisRight.Draw("axis");

        thetaHist->Draw("hist same");

        graphRight.Draw("l same");

        parametrizedHistogram.Draw("hist same");

        // Send whatever has been drawn to the canvas file
        canvas.Print(canvasPath.c_str());
      }

      // Close the canvas file
      canvas.Print((canvasPath + "]").c_str());
    }

    std::cout << std::endl;

  }

  return true;
}




int main(int argc, char ** argv)
{
  std::string mcFilePath = "../../analise-likelihood/data/histograms/vertical_hess_QGSII_urqmd_p_1PeV.root";

  std::vector<std::string> arguments(argc-1,"");

  for (int iArg = 1; iArg < argc; iArg++) {
    arguments[iArg-1] = argv[iArg];
  }

  // MakeTest(DoTestParameters, "Parameters");
  // MakeTest(DoTestSmallAngleIntegral, "Small Angle Integral");
  // MakeTest(DoTestLargeAngleIntegral, "Large Angle Integral");
  MakeTest(DoTestAgainstMC, "Monte Carlo", arguments);

  return 0;
}