#ifndef _StPidHistMaker_head
#define _StPidHistMaker_head
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "TVector3.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoDstMaker;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TTree;

class StCFMult;
class TpcShiftTool;
class TriggerTool;
class MeanDcaTool;
class CentCorrTool;
class VtxShiftTool;

class StPidHistMaker : public StMaker {
	public:
		StPidHistMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName="tofMatchTree.root");
		virtual ~StPidHistMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

	private:
		StPicoDstMaker *mPicoDstMaker;
		StPicoDst      *mPicoDst;
		StPicoEvent	   *event;
		StPicoTrack    *picoTrack;

		StCFMult* mtMult;
		TpcShiftTool* mtShift;
		CentCorrTool* mtCent;
		MeanDcaTool* mtDca;
		TriggerTool* mtTrg;
		VtxShiftTool* mtVtx;

		TH3F* hPro;
		TH3F* hPbar;

		TString mOutputName;
		TFile* mFileOut;

		ClassDef(StPidHistMaker, 1)
};

ClassImp(StPidHistMaker)

#endif
