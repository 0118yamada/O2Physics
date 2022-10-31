//test2
#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HistRegistry {
//configurable for sigmacut
    Configurable<int> tpcsigmacut{"tpcsigmacut", 3, "cut of tpcsigma"};

  HistogramRegistry registry{
    "registry",
    {
      //rammda invariant mass distribution
      {"hpipmass", "", {HistType::kTH1F, {{50, 0, 5}}}}
    }
  };

 
  void process(soa::Join<aod::Tracks, aod::pidTPCPi,  aod::pidTPCPr> const& tracks) {
//get pt pi and pr & E caluculate
  double ptPix[100];
  double ptPiy[100];
  double ptPiz[100];
  
  double ptPrx[100];
  double ptPry[100];
  double ptPrz[100];
 
  double EPi[100];
  double EPr[100];

  double mPi = 0.13957061;
  double mPr = 0.938272088;

  int nPi = 0;
  int nPr = 0;

   for (auto& track: tracks) {
    //Pi
    if (abs(track.tpcNSigmaPi()) < tpcsigmacut && track.sign() == -1) {
	ptPix[nPi] = track.px();
	ptPiy[nPi] = track.py();
	ptPiz[nPi] = track.pz();
	EPi[nPi] = sqrt(pow(mPi, 2) + pow(ptPix[nPi], 2) + pow(ptPiy[nPi], 2) + pow(ptPiz[nPi], 2));
	nPi ++;
    }

    //Pr
    if (abs(track.tpcNSigmaPr()) < tpcsigmacut && track.sign() == 1) {
	ptPrx[nPr] = track.px();
        ptPry[nPr] = track.py();
        ptPrz[nPr] = track.pz();
        EPr[nPr] = sqrt(pow(mPr, 2) + pow(ptPrx[nPr], 2) + pow(ptPry[nPr], 2) + pow(ptPrz[nPr], 2));
        nPr ++;
    }
   }

//invariant mass calculate
  double mPiPr;

  for (int i = 0; i < nPi; i++){
   for (int j = 0; j < nPr; j++){

   mPiPr = sqrt(pow(EPi[i]+EPr[j], 2)
	   - (pow(ptPix[i]+ptPrx[j], 2) + pow(ptPiy[i]+ptPry[j], 2) + pow(ptPiz[i]+ptPrz[j], 2)));
   registry.get<TH1>(HIST("hpipmass"))->Fill(mPiPr);

   }
  }
 }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HistRegistry>(cfgc)
  };
}


