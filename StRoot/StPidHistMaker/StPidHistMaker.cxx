#include "StPidHistMaker.h"

#include <TMath.h>

#include <algorithm>
#include <fstream>
#include <vector>

#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StThreeVectorF.hh"
#include "Stiostream.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "phys_constants.h"

#include "StRoot/CentCorrTool/CentCorrTool.h"
#include "StRoot/MeanDcaTool/MeanDcaTool.h"
#include "StRoot/TpcShiftTool/TpcShiftTool.h"
#include "StRoot/TriggerTool/TriggerTool.h"
#include "StRoot/StCFMult/StCFMult.h"
#include "StRoot/VtxShiftTool/VtxShiftTool.h"

StPidHistMaker::StPidHistMaker(
	const char* name, 
	StPicoDstMaker* picoMaker,
    const char* outName
) : StMaker(name) {
	mOutputName = outName;
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
}

StPidHistMaker::~StPidHistMaker() {}

Int_t StPidHistMaker::Init() {
  	mFileOut = new TFile(mOutputName, "recreate");

	hPro = new TH3F(
		"hPro", ";y;p_{T} [GeV/c];n#sigma",
		12, -0.6, 0.6,
		17, 0.3, 2.0,
		400, -20.0, 20.0
	);
	hPbar = new TH3F(
		"hPbar", ";y;p_{T} [GeV/c];n#sigma",
		12, -0.6, 0.6,
		17, 0.3, 2.0,
		400, -20.0, 20.0
	);
	
	hProTpc = new TH3F(
		"hProTpc", ";y;p_{T} [GeV/c];n#sigma",
		16, -0.8, 0.8,
		17, 0.3, 2.0,
		400, -20.0, 20.0
	);
	hPbarTpc = new TH3F(
		"hPbarTpc", ";y;p_{T} [GeV/c];n#sigma",
		16, -0.8, 0.8,
		17, 0.3, 2.0,
		400, -20.0, 20.0
	);

	// initialize costume modules

	// mean dca tool
	mtDca = new MeanDcaTool();
	mtDca->ReadParams();

	// centrality tool
	mtCent = new CentCorrTool();
	mtCent->EnableIndianMethod(true);
	mtCent->ReadParams();

	// multiplicity and shift tool
	mtShift = new TpcShiftTool();
	mtShift->Init();
	mtMult = new StCFMult();
	mtMult->ImportShiftTool(mtShift);

	// trigger tool
	mtTrg = new TriggerTool();

	// vertext shift tool
	mtVtx = new VtxShiftTool();

	return kStOK;
}

//---------------------------------------------------------
Int_t StPidHistMaker::Finish() {
	mFileOut->cd();
	// hPro[nVz]->Write();
	// hPbar[nVz]->Write();
	// for (int i=0; i<nVz; i++) {
	// 	hPro[i]->Write();
	// 	hPbar[i]->Write();
	// }
	hPro->Write();
	hPbar->Write();
	mFileOut->Close();
	return kStOK;
}

void StPidHistMaker::Clear(Option_t* opt) {}

//---------------------------------------------------------------
Int_t StPidHistMaker::Make() {
	if (!mPicoDstMaker) {
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();
	if (!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	if (!mPicoDst) {
		return kStOK;
	}

	// Load event
	event = (StPicoEvent*)mPicoDst->event();
	if (!event) {
		cerr << "Error opening picoDst Event, skip!" << endl;
		return kStOK;
	}

	TVector3 pVtx = event->primaryVertex();
	Double_t vx = pVtx.X();
	Double_t vy = pVtx.Y();
	Double_t vz = pVtx.Z();

	if (fabs(vx) < 1.e-5 && 
		fabs(vy) < 1.e-5 &&
		fabs(vz) < 1.e-5) {
		return kStOK;
	}

	auto vr = mtVtx->GetShiftedVr(vx, vy);
	if (vr >= 1.0 || fabs(vz) > 50.0) {
		return kStOK;
	}

	Int_t runId = event->runId();
	Int_t trgid = mtTrg->GetTriggerID(event);
	if (trgid < 0) { return kStOK; }

	mtMult->make(mPicoDst);
	Int_t refMult = mtMult->mRefMult;
	Int_t tofMult = mtMult->mTofMult;
	Int_t nTofMatch = mtMult->mNTofMatch;
	Int_t nTofBeta = mtMult->mNTofBeta;

	Int_t refMult3 = mtMult->mRefMult3;
	refMult3 = mtCent->GetIndianRefMult3Corr(
		refMult, refMult3, tofMult, nTofMatch, nTofBeta,
		vz, false
	);
	if (refMult3 < 0) { return kStOK; }
	Int_t cent = mtCent->GetCentrality9(refMult3);
	if (cent < 0 || cent >= 9) { return kStOK; }

	// check DCA
	if (!mtDca->Make(mPicoDst)) { return kStOK; }
	if (mtDca->IsBadMeanDcaZEvent(mPicoDst) || mtDca->IsBadMeanDcaXYEvent(mPicoDst)) {
		return kStOK;
	}

	// track loop
  	Int_t nTracks = mPicoDst->numberOfTracks();
	const Float_t mField = event->bField();

	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
		picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
		if (!picoTrack) { continue; }

		if (!picoTrack->isPrimary()) { continue; }

		// Float_t dca = fabs(picoTrack->gDCA(vx, vy, vz));
		StPicoPhysicalHelix helix = picoTrack->helix(mField);
        Double_t dca = fabs(helix.geometricSignedDistance(pVtx));

		TVector3 pmomentum = picoTrack->pMom();
		Float_t p = pmomentum.Mag();
		Float_t pt = pmomentum.Perp();
		Float_t pz = pmomentum.Z();
		Float_t eta = pmomentum.PseudoRapidity();

		Float_t EP = sqrt(p * p + 0.938272 * 0.938272);
		Float_t YP = 0.5 * log((EP + pz) / (EP - pz));
		if (isnan(YP)) { continue; }
		if (pt < 0.3 || pt > 2.0) { continue; }
		if (fabs(YP) > 0.8) { continue; }

    	Int_t nHitsFit = picoTrack->nHitsFit();
    	Int_t nHitsdEdx = picoTrack->nHitsDedx();
    	Int_t nHitsPoss = picoTrack->nHitsMax();
    	Float_t nHitsRatio = nHitsFit*1.0 / nHitsPoss;
    	Float_t nSigProton = picoTrack->nSigmaProton();
    	Int_t charge = (Int_t)picoTrack->charge();
		nSigProton -= mtShift->GetShift(runId, pt, eta);

		if (nHitsdEdx < 5 || nHitsFit < 20 || nHitsRatio < 0.52) { continue; }
		if (dca > 1.0) { continue; }

		// set TOF flag
        Int_t tofId = picoTrack->bTofPidTraitsIndex();
        Int_t btofMatchFlag = 0;
        Double_t beta = -1.0;
        Double_t btofYLocal = -999.0;
		bool hasTof = false;
        if (tofId >= 0) {
            StPicoBTofPidTraits* tofPid = mPicoDst->btofPidTraits(tofId);
            btofMatchFlag = tofPid->btofMatchFlag();
            if (tofPid) {
                beta = tofPid->btofBeta();
				hasTof = true;
                btofYLocal = tofPid->btofYLocal();
                if (beta < 1e-4) { // recalculate time of flight
                    Double_t tof = tofPid->btof();
                    TVector3 btofHitPos = tofPid->btofHitPos();
                    const StThreeVectorF* btofHitsPosSt = new StThreeVectorF(
                        btofHitPos.X(),btofHitPos.Y(),btofHitPos.Z()
                    );
                    const StThreeVectorF* vtxPosSt = new StThreeVectorF(
                        vx, vy, vz
                    );
                    Double_t L = tofPathLength(vtxPosSt, btofHitsPosSt, helix.curvature());
                    beta = tof > 0 ? L / (tof * (C_C_LIGHT/1.e9)) : std::numeric_limits<Float_t>::quiet_NaN(); // note: quiet nan will never pass > N or < N
                }
            }
        }
        Double_t m2 = -999;
        // if (btofMatchFlag > 0 && beta > 1e-5) {
        if (btofMatchFlag > 0 && beta > 0 && fabs(btofYLocal) < 1.8) {
            m2 = p * p * (pow(1.0 / beta, 2) - 1);
        }
		if (
			(hasTof && m2 > 0.6 && m2 < 1.2) ||
			(!hasTof && pt <= 0.4) // only skip TOF quality condition in very low pt
		) {
			if (charge > 0) { hPro->Fill(YP, pt, nSigProton); }
			else { hPbar->Fill(YP, pt, nSigProton); }
		}
		if (charge > 0) { hProTpc->Fill(YP, pt, nSigProton); }
		{ hPbarTpc->Fill(YP, pt, nSigProton); }

	}  // picotracks loop end
	return kStOK;
}
