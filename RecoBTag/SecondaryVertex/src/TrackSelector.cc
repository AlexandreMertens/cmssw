#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "RecoBTag/SecondaryVertex/interface/TrackSelector.h"

using namespace reco; 

TrackSelector::TrackSelector(const edm::ParameterSet &params) :
	minPixelHits(params.getParameter<unsigned int>("pixelHitsMin")),
	minTotalHits(params.getParameter<unsigned int>("totalHitsMin")),
	minPt(params.getParameter<double>("ptMin")),
	maxNormChi2(params.getParameter<double>("normChi2Max")),
	maxJetDeltaR(params.getParameter<double>("jetDeltaRMax")),
	maxDistToAxis(params.getParameter<double>("maxDistToAxis")),
	maxDecayLen(params.getParameter<double>("maxDecayLen")),
	sip2dValMin(params.getParameter<double>("sip2dValMin")),
	sip2dValMax(params.getParameter<double>("sip2dValMax")),
	sip2dSigMin(params.getParameter<double>("sip2dSigMin")),
	sip2dSigMax(params.getParameter<double>("sip2dSigMax")),
	sip3dValMin(params.getParameter<double>("sip3dValMin")),
	sip3dValMax(params.getParameter<double>("sip3dValMax")),
	sip3dSigMin(params.getParameter<double>("sip3dSigMin")),
	sip3dSigMax(params.getParameter<double>("sip3dSigMax")),
	useMvaSelection_(params.existsAs<bool>("useMvaSelection") ?  params.getParameter<bool>("useMvaSelection") : false),
	useVariableJTA_(params.existsAs<bool>("useVariableJTA") ?  params.getParameter<bool>("useVariableJTA") : false)
{
	std::string qualityClass =
			params.getParameter<std::string>("qualityClass");
	if (qualityClass == "any" || qualityClass == "Any" ||
	    qualityClass == "ANY" || qualityClass == "") {
		selectQuality = false;
		quality = reco::TrackBase::undefQuality;
	} else {
		selectQuality = true;
		quality = reco::TrackBase::qualityByName(qualityClass);
	}
	if (useVariableJTA_){
	  varJTApars = {
	    params.getParameter<double>("a_dR"),
	    params.getParameter<double>("b_dR"),
	    params.getParameter<double>("a_pT"),
	    params.getParameter<double>("b_pT"),
	    params.getParameter<double>("min_pT"),  
	    params.getParameter<double>("max_pT"),
	    params.getParameter<double>("min_pT_dRcut"),  
	    params.getParameter<double>("max_pT_dRcut"),
	    params.getParameter<double>("max_pT_trackPTcut") };
	}
        
        std::cout<<"do mva sel: "<<useMvaSelection_<<std::endl;
      
        if (useMvaSelection_){
           trackSelBDTVarMin = params.getParameter<double>("trackSelBDTVarMin");
           weightFile_ = params.existsAs<edm::FileInPath>("weightFile") ? params.getParameter<edm::FileInPath>("weightFile") : edm::FileInPath("RecoBTag/SecondaryVertex/data/TMVAClassification_BDT.weights.xml.gz");
           //weightFile_(params.existsAs<edm::FileInPath>("weightFile") ? params.getParameter<edm::FileInPath>("weightFile") : edm::FileInPath()) 
           std::cout<<"mva min value: "<<trackSelBDTVarMin<<std::endl;
           std::cout<<"xml file     : "<<weightFile_<<std::endl;
           // initialize MVA evaluators
           evaluator_MVA_.reset( new TMVAEvaluator() );
           std::cout<<"reset mva eval "<<std::endl;
           std::vector<std::string> variables({"Track_dz", "Track_length", "Track_dist", "Track_IP2D", "Track_pt", "Track_chi2", "Track_nHitPixel", "Track_nHitAll"});
           std::vector<std::string> spectators;
           evaluator_MVA_->initialize("V:Color:!Silent:Error", "BDT", weightFile_.fullPath(), variables, spectators);
           std::cout<<"initialized mva "<<std::endl;

        }

}

bool
TrackSelector::operator () (const Track &track,
                            const btag::TrackIPData &ipData,
                            const Jet &jet,
                            const GlobalPoint &pv) const
{

 
  bool jtaPassed = false;
  if (useVariableJTA_) {
    jtaPassed = TrackIPTagInfo::passVariableJTA( varJTApars,
					jet.pt(),track.pt(),
					reco::deltaR(jet.momentum(),track.momentum()));
  }
  else  jtaPassed = true;

  return track.pt() >= minPt &&
    reco::deltaR2(jet.momentum(),
		       track.momentum()) < maxJetDeltaR*maxJetDeltaR &&
    jtaPassed &&
    trackSelection(track, ipData, jet, pv);
}

bool
TrackSelector::operator () (const CandidatePtr &track,
                            const btag::TrackIPData &ipData,
                            const Jet &jet,
                            const GlobalPoint &pv) const
{

 
  bool jtaPassed = false;
  if (useVariableJTA_) {
    jtaPassed = TrackIPTagInfo::passVariableJTA( varJTApars,
					jet.pt(),track->pt(),
					reco::deltaR(jet.momentum(),track->momentum()));
  }
  else  jtaPassed = true;

  return track->pt() >= minPt &&
    reco::deltaR2(jet.momentum(),
		       track->momentum()) < maxJetDeltaR*maxJetDeltaR &&
    jtaPassed &&
    trackSelection(*reco::btag::toTrack(track), ipData, jet, pv);
}

bool
TrackSelector::trackSelection(const Track &track,
                              const btag::TrackIPData &ipData,
                              const Jet &jet,
                              const GlobalPoint &pv) const
{
   // MVA stuff
  reco::TrackBase::Point p (pv.x(), pv.y(), pv.z());
  std::map<std::string,float> variables;

  variables["Track_dz"] = track.dz(p);
  variables["Track_length"] = (ipData.closestToJetAxis - pv).mag();
  variables["Track_dist"] = std::abs(ipData.distanceToJetAxis.value());
  variables["Track_IP2D"] = ipData.ip2d.value();
  variables["Track_pt"] = track.pt();
  variables["Track_chi2"] = track.normalizedChi2();
  variables["Track_nHitPixel"] = track.hitPattern().numberOfValidPixelHits();
  variables["Track_nHitAll"] = track.hitPattern().numberOfValidHits();

  std::cout<<"got here, before evaluate"<<std::endl;  
  double mvaValue = evaluator_MVA_->evaluate(variables); 
 
  std::cout<<mvaValue<<std::endl;
 
  return (!selectQuality || track.quality(quality)) &&
    (minPixelHits <= 0 ||
     track.hitPattern().numberOfValidPixelHits() >= (int)minPixelHits) &&
    (minTotalHits <= 0 ||
     track.hitPattern().numberOfValidHits() >= (int)minTotalHits) &&
    track.normalizedChi2() < maxNormChi2 &&
    std::abs(ipData.distanceToJetAxis.value()) <= maxDistToAxis &&
    (ipData.closestToJetAxis - pv).mag() <= maxDecayLen &&
    ipData.ip2d.value()        >= sip2dValMin &&
    ipData.ip2d.value()        <= sip2dValMax &&
    ipData.ip2d.significance() >= sip2dSigMin &&
    ipData.ip2d.significance() <= sip2dSigMax &&
    ipData.ip3d.value()        >= sip3dValMin &&
    ipData.ip3d.value()        <= sip3dValMax &&
    ipData.ip3d.significance() >= sip3dSigMin &&
    ipData.ip3d.significance() <= sip3dSigMax;
}
