////////////////////////////////////////////////////////////////////////
// Class:       UpMuTrigger
// Module Type: filter
// File:        UpMuTrigger_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Framework/Principal/Handle.h"	
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/FindOneP.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "TTimeStamp.h"
#include "TFile.h"
#include "TTree.h"

#include "DDTUtilities/AssociationUtil.h"
#include "DDTBaseDataProducts/BaseProducts.h"
#include "DDTBaseDataProducts/Boundary.h"
#include "DDTBaseDataProducts/BoundaryList.h"
#include "DDTBaseDataProducts/HitList.h"
#include "DDTBaseDataProducts/DAQHit.h"
#include "DDTBaseDataProducts/CompareDAQHit.h"
#include "DDTBaseDataProducts/Track.h"
#include "DDTBaseDataProducts/Track3D.h"
#include "DDTBaseDataProducts/TriggerDecision.h"
#include "DDTBaseDataProducts/Track.h"	

#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TApplication.h>
#include <TNtuple.h>
#include <TMath.h>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <TFile.h>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions.hpp>

namespace novaddt {
  class UpMuTrigger;
}

class novaddt::UpMuTrigger : public art::EDFilter {
  TFile *f;
  TTree *t1;
  TTree *t2;
  TTree *t3;
  TTree *TimingCal;
  TDirectory *cdupmutrigger;

public:
  explicit UpMuTrigger(fhicl::ParameterSet const & p);
  virtual ~UpMuTrigger();
  void endJob() override;
  void beginJob() override;
  virtual bool filter(art::Event & e) override;

private:

  art::ServiceHandle<art::TFileService> tfs_;
  TNtuple* upmu;
  TH1F *fChi2;
  TH1F *fLLR;


  unsigned _prescale;
  std::string _slicelabel;
  std::string _sliceinstance;
  std::string _tracklabel;
  std::string _trackinstance;
  std::string _trackToSlice;

  std::string _detector;
  bool _containedbit;
  bool passed_tr;

  unsigned _trigger_counts, _after_prescale;

  double c0x, c0y, p0, cw, pw; // initialized later based on detector
  double _TrackLen; 
  unsigned _TrackHitsXY;
  unsigned _TrackHitsX;
  unsigned _TrackHitsY; 
  int _dX; 
  int _dY; 
  int _dZ; 
  double _Chi2;
  double _LLR; 
  double _R2X; 
  double _R2Y; 
  double _MinSlope;
  double _MaxSlope;

  //vectors filled for each track, one entry per dcm                                 
  std::vector<float> dcmhits;
  //vectors filled for each track, one entry per hit                                
  std::vector<float> dcmstime;
  std::vector<float> dcmpe;
  std::vector<float> dcmrd;
  std::vector<float> dcmpl;


  Float_t newLLR;
  Int_t eventi;
  Float_t newTrackLen;
  Float_t newTrackHitsX;
  Float_t newTrackHitsY;
  Float_t newTrackHitsXY;
  Float_t newdX;
  Float_t newdY;
  Float_t newdZ;
  Float_t newChi2;
  Float_t newR2X;
  Float_t newR2Y;
  Float_t newMinSlope;
  Float_t newMaxSlope;
  Float_t efficiency_goodtrack;
  Float_t efficiency_pass;
  Float_t newmin_TDC;
  Float_t newmax_TDC;
  Float_t newtns;

  Float_t newLLR_nocut;
  Float_t newTrackLen_nocut;
  Float_t newTrackHitsX_nocut;
  Float_t newTrackHitsY_nocut;
  Float_t newTrackHitsXY_nocut;
  Float_t newdX_nocut;
  Float_t newdY_nocut;
  Float_t newdZ_nocut;
  Float_t newChi2_nocut;
  Float_t newR2X_nocut;
  Float_t newR2Y_nocut;

  Float_t newLLR_gt;
  Float_t newTrackLen_gt;
  Float_t newTrackHitsX_gt;
  Float_t newTrackHitsY_gt;
  Float_t newTrackHitsXY_gt;
  Float_t newdX_gt;
  Float_t newdY_gt;
  Float_t newdZ_gt;
  Float_t newChi2_gt;
  Float_t newR2X_gt;
  Float_t newR2Y_gt;

  //Float_t Chi2_plot;
  //Float_t LLR_plot;
  int n_events          = 0;
  int n_tracks_output   = 0;
  int n_tracks_input    = 0;
  int n_tracks_good_input = 0;
  int n_tracks_output_true =0;
  int n_tracks_output_kill1   = 0;
  int n_tracks_output_kill2   = 0;
  int n_tracks_output_kill3   = 0;
  int n_tracks_output_kill4   = 0;
  int n_tracks_output_kill5   = 0;
  int n_tracks_output_kill6   = 0;
  int n_tracks_output_kill7   = 0;
  int n_tracks_output_kill8   = 0;
  int n_tracks_output_kill9   = 0;
  int n_tracks_output_kill10   = 0;
  int n_tracks_output_kill11   = 0;
  int n_tracks_output_kill12   = 0;
  int n_tracks_output_kill13   = 0;
  int n_tracks_output_kill14   = 0;

  double GetT(double tdc, double dist_to_readout);
  double GetRes(int ADC);
  double GetRes2(int ADC);
  void GetXYZ(DAQHit const& daqhit, TVector3 start, TVector3 end, double *xyzd);

  void LinFit(const std::vector<double>& x, const std::vector<double>& y, double *fitpar);
  void LinFit(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& ye, double *fitpar);
  void LinFitLLR(std::vector<double>& x_hit, std::vector<double>& y_hit, std::vector<double>& y_error, 
		 double& slope, double& chi2, double& P_up, double& P_dn);

  double GetPigtail(DAQHit const& daqhit) const;
  double  getDCMOff(uint16_t cell, uint8_t view);
};

novaddt::UpMuTrigger::UpMuTrigger(fhicl::ParameterSet const & p)
  : _prescale       (p.get<unsigned>("prescale"))
  , _slicelabel     (p.get<std::string>("slice_label"))
  , _sliceinstance  (p.get<std::string>("slice_instance"))
  , _tracklabel     (p.get<std::string>("track_label"))
  , _trackinstance  (p.get<std::string>("track_instance"))
  , _trackToSlice   (p.get<std::string>("TrackToSliceInstanceLabel"))
  , _detector       (p.get<std::string>("detector"))
  , _containedbit   (p.get<bool>("containedbit"))
				 
  , _TrackLen       (p.get<double>("TrackLen"))
  , _TrackHitsXY    (p.get<unsigned>("TrackHitsXY"))
  , _TrackHitsX     (p.get<unsigned>("TrackHitsX"))
  , _TrackHitsY     (p.get<unsigned>("TrackHitsY"))
  , _dX             (p.get<int>("dX"))
  , _dY             (p.get<int>("dY"))
  , _dZ             (p.get<int>("dZ"))
  , _Chi2           (p.get<double>("Chi2"))
  , _LLR            (p.get<double>("LLR"))
  , _R2X            (p.get<double>("R2X"))
  , _R2Y            (p.get<double>("R2Y"))
  , _MinSlope       (p.get<double>("MinSlope"))
  , _MaxSlope       (p.get<double>("MaxSlope"))
{
  std::cout << "\t TrackToSliceInstanceLabel:     " << _trackToSlice << std::endl;
  assert(_detector == "ndos" || _detector == "fd");

  //initialize cell and plane locations
  if (_detector == "ndos") {
    c0x = -119.5;
    c0y = -199.5;
    p0  = 4.57;
    cw = 3.97;
    pw = 6.654;
  } else if (_detector == "fd") {
    c0x = -759.5;
    c0y = -759.5;
    p0  = 4.57;
    cw = 3.97;
    pw = 6.654;
  }

  _trigger_counts=0;
  _after_prescale=0;

  produces<std::vector<TriggerDecision>>();
}

novaddt::UpMuTrigger::~UpMuTrigger()
{
  // Clean up dynamic memory and other resources here.
}

bool novaddt::UpMuTrigger::filter(art::Event & e)	
{
  n_events++;
   std::unique_ptr<std::vector<TriggerDecision>> 
     trigger_decisions(new std::vector<TriggerDecision>);
   bool result = false;
   passed_tr = false;

   art::Handle<std::vector<Track3D>> tracks;
   e.getByLabel(_tracklabel, tracks);


   art::Handle<std::vector<novaddt::HitList> > slices;
   e.getByLabel(_slicelabel, _sliceinstance, slices);

    // get the hit list for each track. Time cosmusing.
    art::FindOneP<novaddt::HitList> fohl(tracks, e, _tracklabel);
    assert(fohl.isValid());
    assert(fohl.size() == tracks->size());

      for(size_t i_track = 0; i_track < tracks->size(); ++i_track){
	//	this->clearNTuple();
	n_tracks_input++;
	eventi = Int_t(i_track);
	novaddt::Track3D track = tracks->at(i_track);
	// get the hit list for this track
	art::Ptr<novaddt::HitList> this_hit_list = fohl.at(i_track);
	
	// use these to define the trigger window
        // safest not to assume this is time sorted
	  uint64_t min_tdc = this_hit_list->front().TDC().val;
	  uint64_t max_tdc = this_hit_list->back().TDC().val;
	  for(size_t i_hit = 0; i_hit < this_hit_list->size(); ++i_hit){
	    if (this_hit_list->at(i_hit).TDC().val < min_tdc) min_tdc = this_hit_list->at(i_hit).TDC().val;
	    if (this_hit_list->at(i_hit).TDC().val > max_tdc) max_tdc = this_hit_list->at(i_hit).TDC().val;
	  } // end of loop on hits

        novaddt::TDC absolute_t0 = min_tdc;

	int StartX = track.Start().X();
	int StartY = track.Start().Y();
	int StartZ = track.Start().Z();

	int EndX = track.End().X();
	int EndY = track.End().Y();
	int EndZ = track.End().Z();
	if(StartY > EndY){
	  std::swap(StartX, EndX);
	  std::swap(StartY, EndY);
	  std::swap(StartZ, EndZ);
	}

	TVector3 length;
	length.SetXYZ(cw*(track.End().X() - track.Start().X()),
		      cw*(track.End().Y() - track.Start().Y()),
		      pw*(track.End().Z() - track.Start().Z())
		      );

	double tr_length = length.Mag();
	double exp_trav_time = tr_length/29.97;   // c=29.97 cm/ns; 1/64E6=15.625 ns

	//	if (tr_length < _TrackLen)           continue; // Select tracks longer than _TrackLen.
	//	if (TMath::Abs(EndX - StartX) < _dX) continue; // Select tracks longer than _dX cells in X.
	//	if (TMath::Abs(EndY - StartY) < _dY) continue; // Select tracks longer than _dY cells in Y.
	//	if (TMath::Abs(EndZ - StartZ) < _dZ) continue; // Select tracks longer than _dZ planes in Z.
	//	if (!track.Is3D())                  continue; // Select only 3D tracks

	int T0 = 0;
	std::vector<double> x_hit;
	std::vector<double> y_hit;
	std::vector<double> zy_hit;
	std::vector<double> zx_hit;

	// Get T0 that corresponds to the minimum cell number in Y-view
	for(size_t i_hit = 0; i_hit < this_hit_list->size(); ++i_hit){
	  uint16_t iC = this_hit_list->at(i_hit).Cell().val;
	  uint16_t iP = this_hit_list->at(i_hit).Plane().val;
	  uint8_t view = this_hit_list->at(i_hit).View().val;
	  int iTDC = static_cast<int>(this_hit_list->at(i_hit).TDC().val - absolute_t0.val);	

	  if(view == daqchannelmap::X_VIEW){
	    x_hit.push_back(iC);
	    zx_hit.push_back(iP);
	  }

	  if(view == daqchannelmap::Y_VIEW){
	    if(iC == StartY)
	      T0 = iTDC;
	    y_hit.push_back(iC);
	    zy_hit.push_back(iP);
	  }
	  
	}

	// Linear fit of the hits in XZ/YZ views. Drop tracks that has bad linear fit


	double fitpar[3];
	LinFit(x_hit, zx_hit, fitpar);
	double r2x = fitpar[2];
	LinFit(y_hit, zy_hit, fitpar);
	double r2y = fitpar[2];

	//	if(r2x < _R2X) continue;
	//	if(r2y < _R2Y) continue;

	// Vectors of expected and measured times for TGraph 
	std::vector<double> eT_vecxy;
	std::vector<double> eT_vecx;
	std::vector<double> eT_vecy;

	std::vector<double> mT_vecxy;
	std::vector<double> mT_vecx;
	std::vector<double> mT_vecy;

	// Time errors
	std::vector<double> mT_vecxy_e;
	std::vector<double> mT_vecxy_eX;
	std::vector<double> mT_vecxy_eY;

	for(size_t i_hit = 0; i_hit < this_hit_list->size(); ++i_hit){
	  uint16_t iC = this_hit_list->at(i_hit).Cell().val;
	  uint16_t iP = this_hit_list->at(i_hit).Plane().val;
	  int hit_adc = this_hit_list->at(i_hit).ADC().val;

	  if (hit_adc < 50) continue;
	
	  // Calculate time of the hit
	  DAQHit const& dqhit = this_hit_list->at(i_hit);
	  double pos[4];
	  TVector3 start = track.Start();
	  TVector3 end = track.End();
	  GetXYZ(dqhit, start, end, pos);
	  double pigtail = GetPigtail(dqhit);
	  int iTDC = static_cast<int>(this_hit_list->at(i_hit).TDC().val - absolute_t0.val);
	  double  meas_time = static_cast<double>(iTDC - T0);
	  meas_time = meas_time + 0.01*static_cast<double>((this_hit_list->at(i_hit).TDC().frac));

	  double tns = GetT(meas_time, pos[3] + pigtail);
	  // Subtract DCM offsects
 	  double DCM_off = getDCMOff(iC, this_hit_list->at(i_hit).View().val);
 	  //when using MC remove this line
	  //  tns = tns - DCM_off;

          //Fill the variable for the tree: hits size, tns, hit_adc                         
          dcmstime.push_back(tns);
          dcmpe.push_back(hit_adc);

	  if( this_hit_list->at(i_hit).View().val == daqchannelmap::X_VIEW){
	    if(iP > TMath::Max(StartZ, EndZ)) continue;
	    if(iP < TMath::Min(StartZ, EndZ)) continue;

	    if(iC > TMath::Max(StartX, EndX)) continue;
	    if(iC < TMath::Min(StartX, EndX)) continue;
	  }
	  else{
	    if(iP > TMath::Max(StartZ, EndZ)) continue;
	    if(iP < TMath::Min(StartZ, EndZ)) continue;

	    if(iC > TMath::Max(StartY, EndY)) continue;
	    if(iC < TMath::Min(StartY, EndY)) continue;
	  }

	  double exp_time = -999;
	  std::string h_view = "y";	   

	  if( this_hit_list->at(i_hit).View().val == daqchannelmap::X_VIEW){
	    h_view = "x";
	    exp_time = exp_trav_time*(iC-StartX)/(EndX-StartX);
	  }
	  else	    	    
	    exp_time = exp_trav_time*(iC-StartY)/(EndY-StartY);

	  if(false && tr_length > 1200. && tns<-1000){

	    std::cout << ", i_track: "  << std::setw(3) << i_track
		      << ",  i_hit: "   << std::setw(3) << i_hit
		      << ",  iP: "      << std::setw(3) << iP
		      << ",  iC: "      << std::setw(3) << iC
		      << ",  T0: "      << std::setw(3) << T0
		      << ",  pigtail: " << std::setw(3) << pigtail
		      << ",  DCM_off: " << std::setw(3) << DCM_off
		      << ",  ADC: "     << std::setw(5) << (int)this_hit_list->at(i_hit).ADC().val
		      << ",  TDC: "     << (uint64_t)this_hit_list->at(i_hit).TDC().val
		      << ",  meas_time: "     << meas_time
		      << ",  tns: "     << tns
		      << ",  exp_time: "     << exp_time
		      << std::endl;
	  }
	  if (pos[3]<0) continue;
	  //not sure why	  assert(pos[3]>=0);	  

	  eT_vecxy.push_back(exp_time);
	  mT_vecxy.push_back(tns);
	  mT_vecxy_e.push_back(GetRes(hit_adc)); // mT_vecxy_e.push_back(GetRes2(hit_adc));
	  mT_vecxy_eX.push_back(GetRes(hit_adc));
	  mT_vecxy_eY.push_back(GetRes(hit_adc));
	  if( this_hit_list->at(i_hit).View().val == daqchannelmap::X_VIEW){
	    eT_vecx.push_back(exp_time);
	  mT_vecx.push_back(tns);}
	    if( this_hit_list->at(i_hit).View().val == daqchannelmap::Y_VIEW){
	    mT_vecy.push_back(tns);
	  eT_vecy.push_back(exp_time);}
	} // end of loop on hits

        dcmhits.push_back(dcmstime.size());
	//	cdupmutrigger->cd();
	//	TimingCal->Fill();

	unsigned nxhit = mT_vecx.size();
	unsigned nyhit = mT_vecy.size();

	// Remove tracks with less than fTrackHits hits in X+Y, X, Y views
	//	if(nxhit < _TrackHitsX) continue;
	//	if(nyhit < _TrackHitsY) continue;
	//if((nxhit+nyhit) < _TrackHitsXY) continue;
	n_tracks_good_input++;
	double slope_xy, chi_xy, P_up_xy, P_dn_xy;
	double slope_xyX, chi_xyX, P_up_xyX, P_dn_xyX;
	double slope_xyY, chi_xyY, P_up_xyY, P_dn_xyY;

	LinFitLLR(eT_vecxy, mT_vecxy, mT_vecxy_e, slope_xy, chi_xy, P_up_xy, P_dn_xy);
	//      	double LLR = log(P_up_xy/P_dn_xy);
       	double LLR = log(P_dn_xy/P_up_xy); //for coscmic ray
	
	if(mT_vecx.size() == 0 ) continue;
        LinFitLLR(eT_vecx, mT_vecx, mT_vecxy_eX, slope_xyX, chi_xyX, P_up_xyX, P_dn_xyX);
	//        double LLRX = log(P_up_xyX/P_dn_xyX);                                 
	double LLRX = log(P_dn_xyX/P_up_xyX); //for coscmic ray 
	
	if(mT_vecy.size() == 0 ) continue;
        LinFitLLR(eT_vecy, mT_vecy, mT_vecxy_eY, slope_xyY, chi_xyY, P_up_xyY, P_dn_xyY);
        //double LLRY = log(P_up_xyY/P_dn_xyY);        
	double LLRY = log(P_dn_xyY/P_up_xyY); //for coscmic ray 
	
	if (chi_xy > 800) continue;
        newLLR_nocut = Float_t(LLR);
        newTrackLen_nocut  = Float_t(tr_length);
        newTrackHitsX_nocut = Float_t(nxhit);
        newTrackHitsY_nocut = Float_t(nyhit);
        newTrackHitsXY_nocut = Float_t(nyhit);
        newdX_nocut = Float_t(TMath::Abs(EndX - StartX));
        newdY_nocut = Float_t(TMath::Abs(EndY - StartY));
        newdZ_nocut = Float_t(TMath::Abs(EndZ - StartZ));
        newChi2_nocut = Float_t(chi_xy);
        newR2X_nocut = Float_t(r2x);
        newR2Y_nocut = Float_t(r2y);
        //cdupmutrigger->cd();
	//t3->Fill();

	/*        n_tracks_output_kill1++;
        if (tr_length < _TrackLen)           continue;
        n_tracks_output_kill2++;
        if (TMath::Abs(EndX - StartX) < _dX) continue; 
        n_tracks_output_kill3++;
        if (TMath::Abs(EndY - StartY) < _dY) continue;
        n_tracks_output_kill4++;
        if (TMath::Abs(EndZ - StartZ) < _dZ) continue;
        n_tracks_output_kill5++;
        if (!track.Is3D())                  continue;
        n_tracks_output_kill6++;
        if(r2x < _R2X) continue;
        n_tracks_output_kill7++;
        if(r2y < _R2Y) continue;
        n_tracks_output_kill8++;
        if(nxhit < _TrackHitsX) continue;
        n_tracks_output_kill9++;
        if(nyhit < _TrackHitsY) continue;
        n_tracks_output_kill10++;
        if((nxhit+nyhit) < _TrackHitsXY) continue;
        n_tracks_output_kill11++;
        newLLR_gt = Float_t(LLR);
        newTrackLen_gt  = Float_t(tr_length);
        newTrackHitsX_gt = Float_t(nxhit);
        newTrackHitsY_gt = Float_t(nyhit);
        newTrackHitsXY_gt = Float_t(nyhit);
        newdX_gt = Float_t(TMath::Abs(EndX - StartX));
        newdY_gt = Float_t(TMath::Abs(EndY - StartY));
        newdZ_gt = Float_t(TMath::Abs(EndZ - StartZ));
        newChi2_gt = Float_t(chi_xy);
        newR2X_gt = Float_t(r2x);
        newR2Y_gt = Float_t(r2y);
	// cdupmutrigger->cd();
	// t2->Fill();
	// Remove tracks with LLR smaller than input in the config fcl file
	if(LLR < _LLR) continue;
        n_tracks_output_kill12++;
	// Remove tracks with Chi2 fit values higher than input in the config fcl file
	if(chi_xy > _Chi2) continue;
        n_tracks_output_kill13++;
	// Remove tracks with Slope fit values lower/higher than input in the config fcl file
	if(slope_xy < _MinSlope) continue;//inverted
	if(slope_xy > _MaxSlope) continue;//inverted
	n_tracks_output_kill14++;*/
        passed_tr = true;

	/*        newLLR = Float_t(LLR);
        newTrackLen = Float_t(tr_length);
        newTrackHitsX = Float_t(nxhit);
        newTrackHitsY = Float_t(nyhit);
        newTrackHitsXY = Float_t(nyhit);
        newdX = Float_t(TMath::Abs(EndX - StartX));
        newdY = Float_t(TMath::Abs(EndY - StartY));
        newdZ = Float_t(TMath::Abs(EndZ - StartZ));
        newChi2 = Float_t(chi_xy);
        newMaxSlope = Float_t(slope_xy);
        newMinSlope = Float_t(slope_xy);
        newR2X = Float_t(r2x);
        newR2Y = Float_t(r2y);
        efficiency_pass = Int_t(n_tracks_output_kill12);
        efficiency_goodtrack = Int_t(n_tracks_good_input);
	cdupmutrigger->cd();
	t1->Fill();*/

	if(false){
	  n_tracks_output++;
	  std::cout<< std::setiosflags(std::ios::fixed)
		   << "Run: "         << std::setprecision(4) << std::setw(8) << std::left << e.id().run()
		   << ", Event: "     << std::setprecision(4) << std::setw(8) << std::left << e.id().event()
		   << ", LLR: "       << std::setprecision(4) << std::setw(5) << std::left << LLR 
		   << ", Chi2: "      << std::setprecision(4) << std::setw(5) << std::left << chi_xy
		   << ", slope: "     << std::setprecision(4) << std::setw(5) << std::left << slope_xy
		   << ", tr_length: " << std::setprecision(3) << std::setw(3) << std::left << tr_length
		   << ", dX: "        << std::setprecision(4) << std::setw(4) << std::left << TMath::Abs(EndX - StartX)
		   << ", dY: "        << std::setprecision(4) << std::setw(4) << std::left << TMath::Abs(EndY - StartY)
		   << ", dZ: "        << std::setprecision(4) << std::setw(4) << std::left << TMath::Abs(EndZ - StartZ)
		   << ", 3d: "        << std::setprecision(1) << std::setw(1) << std::left << track.Is3D()
		   << ", nxhit: "     << std::setprecision(4) << std::setw(3) << std::left << nxhit
		   << ", nyhit: "     << std::setprecision(4) << std::setw(3) << std::left << nyhit
		   << ", nxyhit: "    << std::setprecision(4) << std::setw(3) << std::left << nxhit+nyhit
		   << ", r2x: "       << std::setprecision(4) << std::setw(5) << std::left << r2x
		   << ", r2y: "       << std::setprecision(4) << std::setw(5) << std::left << r2y
		   << std::endl;
	}

       
	_trigger_counts++;
        //fChi2->Fill(chi_xy);
        //fLLR->Fill(LLR);

        //float track_entries[2]=
        //  {  (float)LLR,
        //     (float)chi_xy   };
        //fupmu->Fill(track_entries);

	if (_trigger_counts%_prescale == _prescale-1) {
	  _after_prescale++;
	  
	  if(_containedbit)
	    trigger_decisions->emplace_back
	      (min_tdc, max_tdc - min_tdc, daqdataformats::trigID::TRIG_ID_DATA_CONTAINED, _prescale);
	  else
	     trigger_decisions->emplace_back
	      (min_tdc, max_tdc - min_tdc, daqdataformats::trigID::TRIG_ID_DATA_UPMU, _prescale);

	  result = true;
	  n_tracks_output_true++;
	  //	  break; // don't need to loop any more on this slice
	} // end if (prescale)
          float track_entries[22]=
	    { (float)e.id().run(), (float)e.id().subRun(),
	      (float)e.id().event() , (float)i_track,
	      (float)(nxhit+nyhit),
	      (float)LLR,
	      (float)chi_xy, (float)slope_xy,
              (float)LLRX,
              (float)chi_xyX, (float)slope_xyX,
              (float)LLRY,
              (float)chi_xyY, (float)slope_xyY,	     
 (float)r2x, (float)r2y,
	      (float)TMath::Abs(EndX - StartX), (float)TMath::Abs(EndY - StartY), (float)TMath::Abs(EndZ - StartZ),
	      (float)nxhit,(float)nyhit, tr_length
	    };
	  upmu->Fill(track_entries);
      } //end loop over tracks

  e.put(std::move(trigger_decisions));
  return result;

} // end Filter()
//----------------------------------------------------------------------------------
void novaddt::UpMuTrigger::beginJob()
{
  upmu = tfs_->make<TNtuple>("upmu_ntuple", "Track Ntuple","Run:SubRun:Event:TrackID:Nhits:LLR:Chi2:Slope:LLRX:Chi2X:SlopeX:LLRY:Chi2Y:SlopeY:R2X:R2Y:dX:dY:dZ:TrackHitsX:TrackHitsY:Length");

  //  f = new TFile("ht_MC_cosmic_newgain.root","recreate");
  //TDirectory  cdupmutrigger = f->mkdir("upmutrigger");
  // cdupmutrigger->cd();
  // t1 = new TTree("t1","/upmutrigger");
  // t2 = new TTree("t2","/upmutrigger");
  /* t3 = new TTree("t3","/upmutrigger");
  TimingCal = new TTree("TimingCal","/upmutrigger");
  t1->Branch("newLLR",&newLLR,"newLLR/F");
  t1->Branch("eventi",&eventi,"eventi/I");
  t1->Branch("newTrackLen",&newTrackLen,"newTrackLen/F");
  t1->Branch("newTrackHitsX",&newTrackHitsX,"newTrackHitsX/F");
  t1->Branch("newTrackHitsY",&newTrackHitsY,"newTrackHitsY/F");
  t1->Branch("newdX",&newdX,"newdX/F");
  t1->Branch("newdY",&newdY,"newdY/F");
  t1->Branch("newdZ",&newdZ,"newdZ/F");
  t1->Branch("newChi2",&newChi2,"newChi2/F");
  t1->Branch("newR2X",&newR2X,"newR2X/F");
  t1->Branch("newR2Y",&newR2Y,"newR2Y/F");
  t1->Branch("efficiency_goodtrack",&efficiency_goodtrack,"efficiency_goodtrack/F");
  t1->Branch("efficiency_pass",&efficiency_pass,"efficiency_pass/F");
  t1->Branch("passed_tr",&passed_tr,"passed_tr/O");

  t2->Branch("newLLR_gt",&newLLR_gt,"newLLR_gt/F");
  t2->Branch("newTrackLen_gt",&newTrackLen_gt,"newTrackLen_gt/F");
  t2->Branch("newTrackHitsX_gt",&newTrackHitsX_gt,"newTrackHitsX_gt/F");
  t2->Branch("newTrackHitsY_gt",&newTrackHitsY_gt,"newTrackHitsY_gt/F");
  t2->Branch("newdX_gt",&newdX_gt,"newdX_gt/F");
  t2->Branch("newdY_gt",&newdY_gt,"newdY_gt/F");
  t2->Branch("newdZ_gt",&newdZ_gt,"newdZ_gt/F");
  t2->Branch("newChi2_gt",&newChi2_gt,"newChi2_gt/F");
  t2->Branch("newR2X_gt",&newR2X_gt,"newR2X_gt/F");
  t2->Branch("newR2Y_gt",&newR2Y_gt,"newR2Y_gt/F");

  t3->Branch("newLLR_nocut",&newLLR_nocut,"newLLR_nocut/F");
  t3->Branch("newTrackLen_nocut",&newTrackLen_nocut,"newTrackLen_nocut/F");
  t3->Branch("newTrackHitsX_nocut",&newTrackHitsX_nocut,"newTrackHitsX_nocut/F");
  t3->Branch("newTrackHitsY_nocut",&newTrackHitsY_nocut,"newTrackHitsY_nocut/F");
  t3->Branch("newdX_nocut",&newdX_nocut,"newdX_nocut/F");
  t3->Branch("newdY_nocut",&newdY_nocut,"newdY_nocut/F");
  t3->Branch("newdZ_nocut",&newdZ_nocut,"newdZ_nocut/F");
  t3->Branch("newChi2_nocut",&newChi2_nocut,"newChi2_nocut/F");
  t3->Branch("newR2X_nocut",&newR2X_nocut,"newR2X_nocut/F");
  t3->Branch("newR2Y_nocut",&newR2Y_nocut,"newR2Y_nocut/F");
*/

  //  TimingCal->Branch("dcmstime",&dcmstime);
  // TimingCal->Branch("dcmpe",&dcmpe);
  // TimingCal->Branch("dcmhits",&dcmhits);
 }
//----------------------------------------------------------------------------------
void novaddt::UpMuTrigger::endJob()
{
  // f->cd();
  // t1->Write();
  // t2->Write();
  // t3->Write();
  // TimingCal->Write();
  std::cout << "=== novaddt::UpmuTrigger endJob"   << std::endl;

  /*  std::cout << "\tNumber of events:             " << n_events           << std::endl;
  std::cout << "\tNumber of input track:        " << n_tracks_input     << std::endl;
  std::cout << "\tNumber of good tracks in input:"<< n_tracks_good_input << std::endl;
  std::cout << "\tNumber of tracks output:      " << n_tracks_output    << std::endl;
  std::cout << "\tNumber of tracks that passed: " << n_tracks_output_true << std::endl;
  std::cout << "\tNumber of tracks output cut 1:      " << n_tracks_output_kill1    << std::endl;
  std::cout << "\tNumber of tracks output cut 2:      " << n_tracks_output_kill2    << std::endl;
  std::cout << "\tNumber of tracks output cut 3:      " << n_tracks_output_kill3    << std::endl;
  std::cout << "\tNumber of tracks output cut 4:      " << n_tracks_output_kill4    << std::endl;
  std::cout << "\tNumber of tracks output cut 5:      " << n_tracks_output_kill5    << std::endl;
  std::cout << "\tNumber of tracks output cut 6:      " << n_tracks_output_kill6    << std::endl;
  std::cout << "\tNumber of tracks output cut 7:      " << n_tracks_output_kill7    << std::endl;
  std::cout << "\tNumber of tracks output cut 8:      " << n_tracks_output_kill8    << std::endl;
  std::cout << "\tNumber of tracks output cut 9:      " << n_tracks_output_kill9    << std::endl;
  std::cout << "\tNumber of tracks output cut 10:      " << n_tracks_output_kill10    << std::endl;
  std::cout << "\tNumber of tracks output cut 11:      " << n_tracks_output_kill11    << std::endl;
  std::cout << "\tNumber of tracks output cut 12:      " << n_tracks_output_kill12    << std::endl;
  std::cout << "\tNumber of tracks output cut 13:      " << n_tracks_output_kill13    << std::endl;
  std::cout << "\tNumber of tracks output cut 14:      " << n_tracks_output_kill14    << std::endl;
  std::cout << "\tNumber of trigger counts" <<_trigger_counts << std::endl;*/
}
//----------------------------------------------------------------------------------
void novaddt::UpMuTrigger::GetXYZ(DAQHit const& daqhit, TVector3 start, TVector3 end, double *xyzd) {
  UInt_t plane = daqhit.Plane().val;
  UInt_t cell  = daqhit.Cell().val;

  // The parameters are taken from the fit of plane vs Z
  double x = c0x + cw*cell;
  double y = x;
  double z = p0 + pw*plane;
  double d = -999.9;

  double x0 = c0x + cw*start.X(); double x1 = c0x + cw*end.X();  
  double y0 = c0y + cw*start.Y(); double y1 = c0y + cw*end.Y();    
  double z0 = p0 + pw*start.Z(); double z1 = p0 + pw*end.Z();

  if(daqhit.View().val == daqchannelmap::X_VIEW){
    y = (y1 - y0)/(z1 - z0);
    y = y*(z - z0);
    y = y + y0;
    d = 800 - y;
  }
  if(daqhit.View().val == daqchannelmap::Y_VIEW){
    x = (x1 - x0)/(z1 - z0);
    x = x*(z - z0);
    x = x + x0;
    d = 800 - x;
  }

  if(d<0) d=0;
  xyzd[0] = x;
  xyzd[1] = y;
  xyzd[2] = z;
  xyzd[3] = d;
}

double novaddt::UpMuTrigger::GetT(double tdc, double dist_to_readout) {
  // This needs to be calibrated and put in a DB table
  double speedOfFiberTransport = 15.3; // cm/ns, "first principles" calib.
  double TDC_to_ns = 15.625; // 1/64E6=15.625 ns
  
  // Differs from c/n due to zigzag
  // paths down fiber.  But, this
  // value may be way off (10%? 30%?).
  return (tdc*TDC_to_ns - dist_to_readout/speedOfFiberTransport);
}



double novaddt::UpMuTrigger::GetRes(int ADC) {

  //The fit function from Evan: f(x) = p0/(p1+x^p2) + p3;

  double p0 = 70267.4;
  double p1 = 677.305;
  double p2 = 1.48702;
  double p3 = 8.46997;
  
//  double p0 = 24004.5;
  //  double p1 = 153.621 ;
  // double p2 = 1.31519;
  // double p3 = 6.35691;

  // double p0 = 137796; //70267.4;                            
  //  double p1 = 1533.47; //677.305;                        
  //  double p2 = 2.05448; //1.48702;                   
   // double p3 = 8.39053; //8.46997; 


  double res = p0/(p1+ pow( static_cast<double>(ADC), p2) ) + p3;  
 return res;

}

double novaddt::UpMuTrigger::GetRes2(int ADC)
{
  //The fit function from Evan: f(x) = p0/(p1+x^p2) + p3;
  
  //1.273e+06
  double p0 = 1273000;
  
  //1.781e+04
  double p1 = 17810;
  
  // 
  double p2 = 2.09;
  double p3 = 8.246;
  
  double res = p0/(p1 + pow( static_cast<double>(ADC)/0.53, p2) ) + p3;
  return res;
}

double novaddt::UpMuTrigger::GetPigtail(DAQHit const& daqhit) const{
    if (_detector == "fd") {
      UInt_t cellid  = daqhit.Cell().val;
      //NOTE: Do we think 32 cells per module can ever change
      //This number can be found in DAQChannelMap but did not want to add that
      //dependancy to geometry unless needed.
      const int kCellsPerModule = 32;
      cellid = cellid % kCellsPerModule;
      // In vertical planes, high cell numbers have longer fibres.
      // In horizontal planes it's the opposite, so flip it round.
      if(daqhit.View().val == daqchannelmap::X_VIEW)
        cellid = kCellsPerModule-cellid-1;

      // This should really never happen, but just to be safe...
      if(cellid < 0 || cellid >= kCellsPerModule) return 100;
      // Email from Tom Chase 2011-04-29
      // NB: this isn't just a ~3.8cm change per cell. Every 8 cells something
      // different happens.
      const double kPigtails[kCellsPerModule] = {
        34.5738, 38.4379, 42.3020, 46.1660, 50.0301, 53.8942, 57.7582, 61.6223,
        64.7504, 68.6144, 72.4785, 76.3426, 80.2067, 84.0707, 87.9348, 91.0790,
        95.3301, 99.1941, 103.058, 106.922, 110.786, 114.650, 118.514, 122.379,
        125.507, 129.371, 133.235, 137.099, 140.963, 144.827, 148.691, 150.751
      };
      return kPigtails[cellid];
    }
    else
      return 0; //to be updated as necessary
}


void novaddt::UpMuTrigger::LinFit(const std::vector<double>& x_val, const std::vector<double>& y_val, double *fitpar){
  // http://en.wikipedia.org/wiki/Simple_linear_regression
  const auto n    = x_val.size();
  const auto x  = (std::accumulate(x_val.begin(), x_val.end(), 0.0))/n;                       // <x>
  const auto y  = (std::accumulate(y_val.begin(), y_val.end(), 0.0))/n;                       // <y>
  const auto xx = (std::inner_product(x_val.begin(), x_val.end(), x_val.begin(), 0.0))/n;     // <xx>
  const auto yy = (std::inner_product(y_val.begin(), y_val.end(), y_val.begin(), 0.0))/n;     // <yy>
  const auto xy = (std::inner_product(x_val.begin(), x_val.end(), y_val.begin(), 0.0))/n;     // <xy>

  const auto b    = (xy - x*y)/(xx - x*x);                                                    // slope
  const auto a    = y - b*x;                                                                  // intercept
  const auto r    = (xy - x*y)/sqrt((xx - x*x)*(yy - y*y));                                   // Rxy - coeffcient of determination
  fitpar[0] = a;
  fitpar[1] = b;
  fitpar[2] = r*r;

}

void novaddt::UpMuTrigger::LinFit(const std::vector<double>& x_val, const std::vector<double>& y_val, const std::vector<double>& y_err, double *fitpar){
  // http://en.wikipedia.org/wiki/Simple_linear_regression
  int n    = x_val.size();
  double x = 0; 
  double y = 0; 
  double xx = 0;
  double yy = 0;
  double xy = 0;
  double ee = 0;

  for ( int i=0; i<n; i++ ){

    x = x + x_val[i]/y_err[i]/y_err[i];
    y = y + y_val[i]/y_err[i]/y_err[i];

    xx = xx + x_val[i]*x_val[i]/y_err[i]/y_err[i];
    yy = yy + y_val[i]*y_val[i]/y_err[i]/y_err[i];
    xy = xy + x_val[i]*y_val[i]/y_err[i]/y_err[i];
    ee = ee + 1./y_err[i]/y_err[i];
  }

  const auto b    = (xy*ee - x*y)/(xx*ee - x*x);                                              // slope
  const auto a    = (xy - b*xx)/x;                                                            // intercept
  const auto r    = (xy - x*y)/sqrt((xx - x*x)*(yy - y*y));                                   // Rxy - coeffcient of determination
  fitpar[0] = a;
  fitpar[1] = b;
  fitpar[2] = r*r;


}

void novaddt::UpMuTrigger::LinFitLLR(std::vector<double>& x_hit, std::vector<double>& y_hit, std::vector<double>& y_error, 
				 double& slope, double& chi2, double& P_up, double& P_dn) {

  int n = x_hit.size();

  double fitpar[3];
  LinFit(x_hit, y_hit, y_error, fitpar);	
  double a = fitpar[0];
  double b = fitpar[1];


  int totAllowOutlier = static_cast<int>(0.1*n);
  const double sigma_cut = 5;
  int iOutlier = 0;
  std::vector<double> x_filt, y_filt, ye_filt, x_filt_e;
  for (int i = 0; i < n; i++) {
    double y_fit = a + b*x_hit[i];
    if ( fabs(y_hit[i] - y_fit) < sigma_cut*y_error[i] || iOutlier>totAllowOutlier) 
      { // > 5 * sigma. Accept 10% of 5sigma outliers. These numbers should be optimized.
	x_filt.push_back(x_hit[i]);
	y_filt.push_back(y_hit[i]);
	ye_filt.push_back(y_error[i]);
	x_filt_e.push_back(0);
      }    
    else
      iOutlier++;

  }

  LinFit(x_filt, y_filt, ye_filt, fitpar);	
  a = fitpar[0];
  b = fitpar[1];
  n = x_filt.size();
  if (n < 5){
    slope = 0;
    chi2  = 999;
    P_dn   = 1e-30;
    P_dn   = 1e-30;
    return;
  }


  slope = fitpar[1];
  
  // chi2
  chi2=0.0;
  for (int i=0; i<n; i++) {
    double y_expected = a + b*x_filt.at(i);
    chi2+=pow((y_filt.at(i)-y_expected) / ye_filt.at(i), 2.0);
  }

  chi2/=static_cast<double>(n-2);      // NDF = N points - 2 parameters in the fit


  // Calculate up/down intercepts
  double one_over_ee = 0.0; 
  double x_over_ee   = 0.0;
  double y_over_ee   = 0.0;

  for (int i=0; i<n; i++) {
    double e = ye_filt.at(i);
    one_over_ee+=1.0/e/e;
    x_over_ee+=x_filt.at(i)/e/e;
    y_over_ee+=y_filt.at(i)/e/e;
  }

  // up/down intercepts defined below
  double up_inter   = (y_over_ee-x_over_ee)/one_over_ee;
  double down_inter = (y_over_ee+x_over_ee)/one_over_ee;


  // Calculate up/down chi2  			
  double up_chi2=0.0, down_chi2=0.0;
  for (int i=0;i<n;i++) {

    double e = ye_filt.at(i);
    double up_expected   = up_inter + x_filt.at(i);
    double down_expected = down_inter - x_filt.at(i);
    up_chi2   += pow((y_filt.at(i)-up_expected) / e, 2.0);
    down_chi2 += pow((y_filt.at(i)-down_expected) / e, 2.0);
  }
  
  // if statement is here for safety to prevent gamma_q failing

    double prob_up   = boost::math::gamma_q(static_cast<double>(n-2)/2.0,up_chi2/2.0);
    double prob_down = boost::math::gamma_q(static_cast<double>(n-2)/2.0,down_chi2/2.0);

    // Use Root TMath if you prefer. Root and boost libraries give identical results.
    // double prob_up = TMath::Prob(up_chi2, n-2);
    // double prob_down = TMath::Prob(down_chi2, n-2);

    if (prob_up<1e-30) prob_up = 1e-30;
    if (prob_down<1e-30) prob_down = 1e-30;

    P_up = prob_up;
    P_dn = prob_down;

    if(false){ 

    TCanvas *can = new TCanvas("can", "can", 800, 600);
    // Debugging part. Will produce the time distribution fits using Root fitter and analytical fuction
    TGraphErrors *eTmT_grx = new TGraphErrors(x_filt.size(), &x_filt[0], &y_filt[0], &x_filt_e[0], &ye_filt[0]);
    TF1 *SlopeFit = new TF1("SlopeFit", "pol1", -10, 150);
    TF1 *UpFit = new TF1("UpFit", "pol1", -10, 150);
    TF1 *UpFit1 = new TF1("UpFit1", "pol1", -10, 150);
    TF1 *DnFit2 = new TF1("DnFit2", "pol1", -10, 150);
    TF1 *DnFit = new TF1("DnFit", "pol1", -10, 150);
    UpFit->SetLineColor(4);
    DnFit->SetLineColor(4);
    DnFit->SetLineStyle(3);
    UpFit->SetParameter(0, up_inter);
    UpFit->SetParameter(1, 1);

    DnFit2->SetParameter(0, fitpar[0]);
    DnFit2->SetParameter(1, fitpar[1]);
    DnFit2->SetLineColor(kGreen);
    DnFit2->SetLineStyle(3);


    DnFit->SetParameter(0, down_inter);
    DnFit->SetParameter(1, -1);
    


    SlopeFit->SetParameter(0, a);
    SlopeFit->SetParameter(1, b);
    SlopeFit->SetLineColor(3);
    SlopeFit->SetLineStyle(3);
    can->Draw();
    can->cd();

    eTmT_grx->GetXaxis()->SetRangeUser(-10, 100);
    eTmT_grx->Draw("AP");
    SlopeFit->Draw("same");
    //    UpFit->Draw("same");
    DnFit2->Draw("same");
    DnFit->Draw("same");

    UpFit1->SetLineColor(kRed);
    UpFit1->SetLineWidth(1);
    //    UpFit1->SetLineStyle(4);

    eTmT_grx->Fit("UpFit1", "BR", "", -10, 150);

    TPaveText *pt = new TPaveText(.75, .9, 0.9, 1.0, "brNDC");    
    pt->SetBorderSize(1);
    pt->AddText(Form("P_up = %.1e",   P_up            ));
    pt->AddText(Form("P_dn = %.1e",   P_dn            ));
    pt->AddText(Form("LLR =  %.1f",   log(P_up/P_dn)  ));  
    pt->AddText(Form("Chi2 = %.2f",   chi2  ));   
    pt->AddText(Form("NDF =  %d",   n-2  )); 
    
 
    pt->Draw();
    can->Update();
    can->SaveAs(Form("plots/Chi2_%f.pdf", chi2));
    //    can->Delete();
  }





}

double novaddt::UpMuTrigger::getDCMOff(uint16_t cell, uint8_t view){

  // Cell number starts from 1 

  if (_detector == "ndos") {
    if (view == daqchannelmap::Y_VIEW)
      if (cell > 31 && cell < 64) return -26;
      else if (cell > 63) return -40;
      else return 0;
    else // X_VIEW
      if (cell > 7 && cell < 32) return -10;
      else if (cell > 31) return 5;
      else return 0;
  }
  else if (_detector == "fd")
    return (cell/64)*40.;
  else return 0;
}

DEFINE_ART_MODULE(novaddt::UpMuTrigger)
