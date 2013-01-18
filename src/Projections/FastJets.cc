// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Math/MathUtils.hh"

namespace Rivet {


  void FastJets::_init1(JetAlgName alg, double rparameter, double seed_threshold) {
    setName("FastJets");
    MSG_DEBUG("R parameter = " << rparameter);
    MSG_DEBUG("Seed threshold = " << seed_threshold);
    if (alg == KT) {
      _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == CAM) {
      _jdef = fastjet::JetDefinition(fastjet::cambridge_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == ANTIKT) {
      _jdef = fastjet::JetDefinition(fastjet::antikt_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == DURHAM) {
      _jdef = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::E_scheme);
    } else {
      // Plugins:
      if (alg == SISCONE) {
        const double OVERLAP_THRESHOLD = 0.75;
        _plugin.reset(new fastjet::SISConePlugin(rparameter, OVERLAP_THRESHOLD));
      } else if (alg == PXCONE) {
        string msg = "PxCone currently not supported, since FastJet doesn't install it by default. ";
        msg += "Please notify the Rivet authors if this behaviour should be changed.";
        throw Error(msg);
        //_plugin.reset(new fastjet::PxConePlugin(rparameter));
      } else if (alg == ATLASCONE) {
        const double OVERLAP_THRESHOLD = 0.5;
        _plugin.reset(new fastjet::ATLASConePlugin(rparameter, seed_threshold, OVERLAP_THRESHOLD));
      } else if (alg == CMSCONE) {
        _plugin.reset(new fastjet::CMSIterativeConePlugin(rparameter, seed_threshold));
      } else if (alg == CDFJETCLU) {
        const double OVERLAP_THRESHOLD = 0.75;
        _plugin.reset(new fastjet::CDFJetCluPlugin(rparameter, OVERLAP_THRESHOLD, seed_threshold));
      } else if (alg == CDFMIDPOINT) {
        const double OVERLAP_THRESHOLD = 0.5;
        _plugin.reset(new fastjet::CDFMidPointPlugin(rparameter, OVERLAP_THRESHOLD, seed_threshold));
      } else if (alg == D0ILCONE) {
        const double min_jet_Et = 6.0;
        _plugin.reset(new fastjet::D0RunIIConePlugin(rparameter, min_jet_Et));
      } else if (alg == JADE) {
        _plugin.reset(new fastjet::JadePlugin());
      } else if (alg == TRACKJET) {
        _plugin.reset(new fastjet::TrackJetPlugin(rparameter));
      }
      _jdef = fastjet::JetDefinition(_plugin.get());
    }
  }

  void FastJets::_init2(fastjet::JetAlgorithm type,
                        fastjet::RecombinationScheme recom, double rparameter) {
    setName("FastJets");
    _jdef = fastjet::JetDefinition(type, rparameter, recom);
  }

  void FastJets::_init3(fastjet::JetDefinition::Plugin* plugin) {
    setName("FastJets");
    /// @todo Should we be copying the plugin?
    _plugin.reset(plugin);
    _jdef = fastjet::JetDefinition(_plugin.get());
  }



  int FastJets::compare(const Projection& p) const {
    const FastJets& other = dynamic_cast<const FastJets&>(p);
    return \
      (_useInvisibles ? mkNamedPCmp(other, "FS") : mkNamedPCmp(other, "VFS")) ||
      cmp(_jdef.jet_algorithm(), other._jdef.jet_algorithm()) ||
      cmp(_jdef.recombination_scheme(), other._jdef.recombination_scheme()) ||
      cmp(_jdef.plugin(), other._jdef.plugin()) ||
      cmp(_jdef.R(), other._jdef.R()) ||
      cmp(_adef, other._adef);
  }



  void FastJets::project(const Event& e) {
    ParticleVector particles;
    if (_useInvisibles) {
      particles = applyProjection<FinalState>(e, "FS").particles();
    } else {
      particles = applyProjection<FinalState>(e, "VFS").particles();
    }
    calc(particles);
  }


  void FastJets::calc(const ParticleVector& ps) {
    _particles.clear();
    vector<fastjet::PseudoJet> vecs;
    // Store 4 vector data about each particle into vecs
    int counter = 1;
    foreach (const Particle& p, ps) {
      const FourMomentum fv = p.momentum();
      fastjet::PseudoJet pJet(fv.px(), fv.py(), fv.pz(), fv.E());
      pJet.set_user_index(counter);
      vecs.push_back(pJet);
      _particles[counter] = p;
      ++counter;
    }
    MSG_DEBUG("Running FastJet ClusterSequence construction");

    // Choose CSeq as basic or area-calculating depending on whether _adef pointer is non-null.
    if (_adef == 0) {
      _cseq.reset(new fastjet::ClusterSequence(vecs, _jdef));
    } else {
      _cseq.reset(new fastjet::ClusterSequenceArea(vecs, _jdef, *_adef));
    }
  }


  Jets FastJets::_pseudojetsToJets(const PseudoJets& pjets) const {
    Jets rtn;
    foreach (const fastjet::PseudoJet& pj, pjets) {
      assert(clusterSeq());
      const PseudoJets parts = clusterSeq()->constituents(pj);
      vector<Particle> constituents;
      constituents.reserve(parts.size());
      foreach (const fastjet::PseudoJet& p, parts) {
        map<int, Particle>::const_iterator found = _particles.find(p.user_index());
        assert(found != _particles.end());
        constituents.push_back(found->second);
      }
      FourMomentum pjet(pj.E(), pj.px(), pj.py(), pj.pz());
      Jet j(constituents, pjet);
      rtn.push_back(j);
    }
    /// @todo Cache?
    return rtn;
  }


  void FastJets::reset() {
    _yscales.clear();
    _particles.clear();
    /// @todo _cseq = fastjet::ClusterSequence();
  }


  size_t FastJets::numJets(double ptmin) const {
    if (_cseq.get() != 0) {
      return _cseq->inclusive_jets(ptmin).size();
    } else {
      return 0;
    }
  }


  Jets FastJets::_jets(double ptmin) const {
    Jets rtn = _pseudojetsToJets(pseudoJets(ptmin));
    return rtn;
  }

  /// D===1/R_{12} \Sum_{i\in J} p_{Ti}/p_{TJ} R_i
  double FastJets::Dipolarity(const fastjet::PseudoJet &j, const double dcut=0.5) const
  {
    assert(clusterSeq());
    double D(0);
    PseudoJets subJets=j.validated_cs()->exclusive_subjets(j,dcut);
    unsigned int nSubJets=j.validated_cs()->n_exclusive_subjets(j,dcut);
    double denom=j.pt()*pow(subJets.at(0).delta_R(subJets.at(1)),2);
    if(nSubJets>2)  {
      //for now take the first two 
      D+=subJets.at(0).pt()*(pow(subJets.at(0).eta(),2) + pow(subJets.at(0).phi(),2));
      D+=subJets.at(1).pt()*(pow(subJets.at(1).eta(),2) + pow(subJets.at(1).phi(),2));
      return D/denom;
    }
    return D;
  }

  ///\vec{t} ===\Sum_{i\in J} |r_i|p_{Ti}/p_{TJ}\vec{r_i}
  std::pair<double,double> FastJets::JetPull(const fastjet::PseudoJet &j, const double ptmin) const {
    assert(clusterSeq());
    const PseudoJets parts = clusterSeq()->constituents(j);
    const double jetRap = j.rapidity(), jetPhi = j.phi();
    double ty=0, tphi=0, tmag=0, ttheta=0, dphi=0;
    foreach (const fastjet::PseudoJet& p, parts) {
      dphi = mapAngleMPiToPi(p.phi()-jetPhi); //don't generate a large pull for jets at 2pi
      if(p.pt() > ptmin) { //pt always > 0, if the user hasn't defined a cut, this will always pass
	double ptTimesRmag=sqrt(pow(p.rapidity()-jetRap,2) + pow(dphi,2))*p.pt();//use dphi
	ty+=ptTimesRmag*(p.rapidity()-jetRap);
	tphi+=ptTimesRmag*(dphi);//use dphi
      }
    }
    //now parametrize \vec{t}=|t|(cos(\theta_t),sin(\theta_t))
    tmag=sqrt(pow(ty/j.pt(),2) + pow(tphi/j.pt(),2));
    if(tmag>0) {
      ttheta=atan2(tphi,ty);
    }
    else
      MSG_ERROR("Couldn't calculate theta_t, tphi:"<<tphi<<" ty:"<<ty);
    
    return std::pair<double,double>(tmag,ttheta);
  }

  ///Q===\Sum_{i\in J} q_i*p_{Ti}^k/p_{TJ} 
  double FastJets::JetCharge(const fastjet::PseudoJet &j, const double k, const double ptmin) const {
    assert(clusterSeq());
    const PseudoJets parts = clusterSeq()->constituents(j);
    double q(0);
    foreach (const fastjet::PseudoJet& p, parts) {
      map<int, Particle>::const_iterator found = _particles.find(p.user_index());
      assert(found != _particles.end());
      if(p.pt() > ptmin) //pt always > 0, if the user hasn't defined a cut, this will always pass
	q += PID::charge(found->second) * pow(p.pt(),k);       
    }
    return q/j.pt();
  }

  // Jets FastJets::jetsByPt(double ptmin) const {
  //   return _pseudojetsToJets(pseudoJetsByPt(ptmin));
  // }


  // Jets FastJets::jetsByE(double ptmin) const {
  //   return _pseudojetsToJets(pseudoJetsByE(ptmin));
  // }


  // Jets FastJets::jetsByRapidity(double ptmin) const {
  //   return _pseudojetsToJets(pseudoJetsByRapidity(ptmin));
  // }


  PseudoJets FastJets::pseudoJets(double ptmin) const {
    if (_cseq.get() != 0) {
      return _cseq->inclusive_jets(ptmin);
    } else {
      return PseudoJets();
    }
  }


  vector<double> FastJets::ySubJet(const fastjet::PseudoJet& jet) const {
    assert(clusterSeq());
    fastjet::ClusterSequence subjet_cseq(clusterSeq()->constituents(jet), _jdef);
    vector<double> yMergeVals;
    for (int i = 1; i < 4; ++i) {
      // Multiply the dmerge value by R^2 so that it corresponds to a
      // relative k_T (fastjet has 1/R^2 in the d_ij distance by default)
      const double ktmerge = subjet_cseq.exclusive_dmerge(i) * _jdef.R()*_jdef.R();
      yMergeVals.push_back(ktmerge/jet.perp2());
    }
    _yscales.insert(make_pair( jet.cluster_hist_index(), yMergeVals ));
    return yMergeVals;
  }



  fastjet::PseudoJet FastJets::splitJet(fastjet::PseudoJet jet, double& last_R) const {
    // Sanity cuts
    if (jet.E() <= 0 || _cseq->constituents(jet).size() <= 1) {
      return jet;
    }

    // Build a new cluster sequence just using the consituents of this jet.
    assert(clusterSeq());
    fastjet::ClusterSequence cs(clusterSeq()->constituents(jet), _jdef);

    // Get the jet back again
    fastjet::PseudoJet remadeJet = cs.inclusive_jets()[0];
    MSG_DEBUG("Jet2:" << remadeJet.m() << "," << remadeJet.e());

    fastjet::PseudoJet parent1, parent2;
    fastjet::PseudoJet split(0.0, 0.0, 0.0, 0.0);
    while (cs.has_parents(remadeJet, parent1, parent2)) {
      MSG_DEBUG("Parents:" << parent1.m() << "," << parent2.m());
      if (parent1.m2() < parent2.m2()) {
        fastjet::PseudoJet tmp;
        tmp = parent1; parent1 = parent2; parent2 = tmp;
      }

      double ktdist = parent1.kt_distance(parent2);
      double rtycut2 = 0.3*0.3;
      if (parent1.m() < ((2.0*remadeJet.m())/3.0) && ktdist > rtycut2*remadeJet.m2()) {
        break;
      } else {
        remadeJet = parent1;
      }
    }

    last_R = 0.5 * sqrt(parent1.squared_distance(parent2));
    split.reset(remadeJet.px(), remadeJet.py(), remadeJet.pz(), remadeJet.E());
    return split;
  }



  fastjet::PseudoJet FastJets::filterJet(fastjet::PseudoJet jet,
                                         double& stingy_R, const double def_R) const {
    assert(clusterSeq());

    if (jet.E() <= 0.0 || clusterSeq()->constituents(jet).size() == 0) {
      return jet;
    }
    if (stingy_R == 0.0) {
      stingy_R = def_R;
    }

    stingy_R = def_R < stingy_R ? def_R : stingy_R;
    fastjet::JetDefinition stingy_jet_def(fastjet::cambridge_algorithm, stingy_R);

    //FlavourRecombiner recom;
    //stingy_jet_def.set_recombiner(&recom);
    fastjet::ClusterSequence scs(clusterSeq()->constituents(jet), stingy_jet_def);
    std::vector<fastjet::PseudoJet> stingy_jets = sorted_by_pt(scs.inclusive_jets());

    fastjet::PseudoJet reconst_jet(0.0, 0.0, 0.0, 0.0);

    for (unsigned isj = 0; isj < std::min(3U, (unsigned int) stingy_jets.size()); ++isj) {
      reconst_jet += stingy_jets[isj];
    }
    return reconst_jet;
  }
  fastjet::JetAlgorithm FastJets::setJetAlgorithm(JetAlgName subJetAlgorithm) const
  {
    //Do we want to support all enums? This is only a subset...
    switch(subJetAlgorithm)
      {
      case KT:
	return fastjet::kt_algorithm;
      case ANTIKT:
	return fastjet::antikt_algorithm;
      case CAM:
	return fastjet::cambridge_algorithm;
      case DURHAM:
	return fastjet::ee_kt_algorithm;
      default: 
	cout << "No plugin jet algorithms accepted for Filter! Prepare to die!" << endl;
	return fastjet::undefined_jet_algorithm;
      }
  }
  
  fastjet::PseudoJet FastJets::Filter(fastjet::PseudoJet jet, JetAlgName subjet_def, 
				      int hardest,double subjet_R=0.3) const {
    assert(clusterSeq());
    //sanity check on the jet
    if (jet.E() <= 0.0 || clusterSeq()->constituents(jet).size() == 0) {
      return jet;
    }
    fastjet::Filter filter(fastjet::JetDefinition(setJetAlgorithm(subjet_def), subjet_R), fastjet::SelectorNHardest(hardest));
    return filter(jet);
  }

  fastjet::PseudoJet FastJets::Trimmer(fastjet::PseudoJet jet, JetAlgName subjet_def, 
				       double percentage, double subjet_R=0.3) const {
    assert(clusterSeq());
    //sanity check on the jet
    if (jet.E() <= 0.0 || clusterSeq()->constituents(jet).size() == 0) {
      return jet;
    }
    fastjet::Filter filter(fastjet::JetDefinition(setJetAlgorithm(subjet_def), subjet_R), fastjet::SelectorPtFractionMin(percentage));
    return filter(jet);
  }
  fastjet::PseudoJet FastJets::Pruner(fastjet::PseudoJet jet, JetAlgName subjet_def, 
						double zcut=0.1, double Rcut_factor=0.5) const {
    //sanity check on the jet
    assert(clusterSeq());
    if (jet.E() <= 0.0 || clusterSeq()->constituents(jet).size() == 0) {
      return jet;
    }
    //NOTE: Pruner sets the jet algorithm R to the maximum allowed, so setting it
    //to 1 here only to follow jetalg syntax
    //finally declare the pruner and apply it to the jet
    fastjet::Pruner pruner(fastjet::JetDefinition(setJetAlgorithm(subjet_def), 1), zcut, Rcut_factor);
    return pruner(jet);
  }
  PseudoJets FastJets::GetAxes(unsigned int n_jets, 
		PseudoJets& inputJets, JetAlgName subjet_def, double subR) const {

    assert(clusterSeq());
    //sanity check
    if (inputJets.size() < n_jets) {
	std::cout << "Not enough input particles." << endl;
	return inputJets;
	}
    //get subjets, return
    fastjet::ClusterSequence sub_clust_seq(inputJets, fastjet::JetDefinition(setJetAlgorithm(subjet_def), subR));
    return sub_clust_seq.exclusive_jets((signed)n_jets);
  }

  double FastJets::TauValue(double beta, double jet_rad, 
	PseudoJets& particles, PseudoJets& axes) const {

    double tauNum = 0.0;
    double tauDen = 0.0;

    if(particles.size() == 0)return 0.0;
   
    for (unsigned int i = 0; i < particles.size(); i++) {

     	 // find minimum distance (set R large to begin)

      double minR = 10000.0;
      for (unsigned int j = 0; j < axes.size(); j++) {
         double tempR = sqrt(particles[i].squared_distance(axes[j]));
         if (tempR < minR) minR = tempR;
      }

	//calculate nominator and denominator
      
      tauNum += particles[i].perp() * pow(minR,beta);
      tauDen += particles[i].perp() * pow(jet_rad,beta);
    }
	
	//return N-subjettiness

    return tauNum/tauDen;
  }

  void FastJets::UpdateAxes(double beta,
	PseudoJets& particles, PseudoJets& axes) const {

    vector<int> belongsto;
    //no reason not to use foreach here
    for (unsigned int i = 0; i < particles.size(); i++) {

     	 // find minimum distance axis

      int assign = 0;
      double minR = 10000.0;
      for (unsigned int j = 0; j < axes.size(); j++) {
         double tempR = sqrt(particles[i].squared_distance(axes[j]));
         if (tempR < minR) {
            minR = tempR;
            assign = j;
	  }
       }
     belongsto.push_back(assign);
    }

	// iterative step
    
    double deltaR2, distphi;
    vector<double> ynom, phinom, den;
    ynom.resize(axes.size());
    phinom.resize(axes.size());
    den.resize(axes.size());

    for (unsigned int i = 0; i < particles.size(); i++) {
      distphi = particles[i].phi() - axes[belongsto[i]].phi();
      deltaR2 = particles[i].squared_distance(axes[belongsto[i]]);

      if(FuzzyEquals(deltaR2, 0.))continue;
      
      if (abs(distphi) <= M_PI) phinom.at(belongsto[i]) += particles[i].perp() * particles[i].phi() * pow(deltaR2, (beta-2)/2);
      else if ( distphi > M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (-2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);
      else if ( distphi < M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (+2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);

      ynom.at(belongsto[i]) += particles[i].perp() * particles[i].rap() * pow(deltaR2, (beta-2)/2);
      
      den.at(belongsto[i]) += particles[i].perp() * pow(deltaR2, (beta-2)/2);
    }
      
    // reset to new axes

    for (unsigned int j = 0; j < axes.size(); j++) {

      if (FuzzyEquals(den[j], 0.)) axes.at(j) = axes[j];
      else {
	axes.at(j).reset_momentum_PtYPhiM(axes[j].perp(), ynom[j] / den[j], fmod( 2*M_PI + (phinom[j] / den[j]), 2*M_PI ), axes[j].perp()/2);
      }
    }
  }

  bool FastJets::FuzzyEquals(double a, double b) const {
    return fabs(a - b) < 0.00000001;
  }

  double FastJets::KeyColToRight(int p, vector<ACFpeak> peaks, vector<double> ASF_erf) const {

    int higherpeak = -1;
    double height = peaks[p].height;
    double keycol = height;

    for (unsigned int k = p+1; k < peaks.size(); k++){

      if (peaks[k].height > height) { higherpeak = k; break; }
      }

    if (higherpeak != -1){

      int startindex = peaks[p].index;
      int endindex = peaks[higherpeak].index;
      for (int j = startindex+1; j < endindex; j++){

        if (ASF_erf[j] < keycol) keycol = ASF_erf[j];
        }
      }
    else keycol = 0.;

    return keycol;
  }

  double FastJets::KeyColToLeft(int p, vector<ACFpeak> peaks, vector<double> ASF_erf) const {

    int higherpeak = -1;
    double height = peaks[p].height;
    double keycol = height;

    for (int k = p-1; k >= 0; k--){

      if (peaks[k].height > height) { higherpeak = k; break; }
      }

    if (higherpeak != -1){

      int endindex = peaks[p].index;
      int startindex = peaks[higherpeak].index;
      for (int j = startindex+1; j < endindex; j++){

        if (ASF_erf[j] < keycol) keycol = ASF_erf[j];
        }
      }
    else keycol = 0.;

    return keycol;
  }
  
  vector<ACFpeak> FastJets::ASFPeaks(PseudoJets& particles,
		unsigned int most_prominent = 0, double minprominence = 0.) const {

    unsigned int meshsize = 500;
    double sigma = 0.06;

    vector<ACFpeak> peaks;

	//sanity check
    if(particles.size() < 2) {
      cout << "Not enough particles in jet for ACF." << endl;
      return peaks;
      }

	//pair all particles up
    vector<ACFparticlepair> pairs;
    ACFparticlepair dummy;
    for(unsigned int k = 0; k < particles.size(); k++) {
      for(unsigned int j = 0; j < k; j++) {
	double phidist = abs(particles[k].phi_std() - particles[j].phi_std());
	if( phidist > M_PI ) phidist = 2 * M_PI - phidist;
        //dummy.deltaR = sqrt( pow(particles[k].pseudorapidity() - particles[j].pseudorapidity(), 2) +
	//	pow(phidist, 2) );
	dummy.deltaR = sqrt(particles[k].plain_distance(particles[j]));
	dummy.weight = particles[k].perp() * particles[j].perp() * dummy.deltaR * dummy.deltaR;
	pairs.push_back(dummy);
	}
      }

	//sort by delta R
    sort(pairs.begin(), pairs.end(), ppsortfunction());

    double Rmax = pairs[pairs.size() - 1].deltaR;

    double rVal = 0.;
    double xArg = 0.;
    vector<double> ACF(meshsize), ASF_gauss(meshsize), Rvals(meshsize), erf_denom(meshsize), gauss_peak(meshsize);
    ACF[0] = 0.;
    ASF_gauss[0] = 0.;
    Rvals[0] = 0.;
    erf_denom[0] = 0.;
    gauss_peak[0] = 0.;
    double fVal = 0.;
    double eVal = 0.;
    double gVal = 0.;

	//mesh loop
    for (unsigned int k = 1; k < meshsize; k++) {

      rVal = (double)k*Rmax/(meshsize-1);
      Rvals[k] = rVal;

	//reset
      fVal = 0.;
      eVal = 0.;
      gVal = 0.;

	//Loop on pairs.
      for (unsigned int j = 0; j < pairs.size(); j++){
	//ACF-Add pairs within mesh's deltaR.
	if (pairs[j].deltaR <= rVal) fVal += pairs[j].weight;

	//Smoothing function argument
	xArg = (double)(rVal-pairs[j].deltaR)/sigma;

	//ASF Error Function Denominator: Add pairs weighted by Erf values.
	eVal += (double)pairs[j].weight*0.5*(1.0+erf(xArg));

	//ASF Gaussian Numerator: Add pairs weight by Gaussian values.
	gVal += (double)pairs[j].weight * (exp( -pow ((abs(xArg)),2.0)));

	}//end pair loop
      ACF[k] = fVal;
      erf_denom[k] = eVal;
      gauss_peak[k] = gVal;
      ASF_gauss[k] = (double)gVal*(1/sqrt(M_PI))*rVal/sigma; //Normalized Gaussian value

      }//end mesh loop

    vector<double> ASF_erf(meshsize), ASF(meshsize);
    ASF[0] = 0.;
    ASF[meshsize-1] = 0.;
    ASF_erf[0] = 0.;

	//Total jet mass
    double jetmass = ACF[meshsize-1];   

	//Second mesh loop
    for (unsigned int k = 1; k < meshsize; k++){

	//Compute simple (spiked) ASF
      if (ACF[k-1] == 0.) {ASF[k] = 0.;}
      else if (k < meshsize - 2) 
        {ASF[k] = (log(ACF[k+1]) - log(jetmass*ACF[k-1]))/(log(Rvals[k+1]) - log(Rvals[k-1]));}

      //Compute gaussian (smoothed) ASF
      ASF_erf[k] = (FuzzyEquals(erf_denom[k],0.)) ? 0. : ASF_gauss[k]/erf_denom[k];
      //if(erf_denom[k] == 0.0)ASF_erf[k] = 0.0;
      //else ASF_erf[k] = ASF_gauss[k]/erf_denom[k];

      //Normalize ACF
      ACF[k] /= jetmass;

      }//end mesh loop

    ACFpeak myPeak;

    double lefth =  0.;
    double height = 0.;
    double righth = 0.;

    for (unsigned int k = 1; k < meshsize-1; k++) {

      lefth  = ASF_erf[k-1];
      height = ASF_erf[k];
      righth = ASF_erf[k+1];

	//Found a peak?
      if (lefth < height && height > righth) {

	myPeak.height = height;
	myPeak.Rval   = Rvals[k];
	myPeak.index  = k;

	//Peaks are stored according to
	//index (equilvalent to R values)
	//from lowest to highest
	peaks.push_back(myPeak);
	}

      }

	//Partial mass of peak
    for (unsigned int p = 0; p < peaks.size(); p++)peaks[p].partialmass = (double)sqrt(gauss_peak[peaks[p].index]);

	//peaks[p].partialmass = (double)sqrt(sqrt(M_PI)*sigma*peaks[p].height*ACF[peaks[p].index]*jetmass/peaks[p].Rval);

	//Prominence of peak
    for (unsigned int k = 0; k < peaks.size(); k++){
      double height = peaks[k].height;
      double leftdescent = height - KeyColToLeft(k, peaks, ASF_erf);
      double rightdescent = height - KeyColToRight(k, peaks, ASF_erf);
      if (leftdescent < rightdescent) peaks[k].prominence = leftdescent;
      else peaks[k].prominence = rightdescent;
      }

	//return all peaks
    if(most_prominent == 0 && FuzzyEquals(minprominence, 0.) ) {
      return peaks;
      }

	//return all peaks over given prominence
    else if(most_prominent == 0) {
      vector<ACFpeak> dummyp;
      for(unsigned int i = 0; i < peaks.size(); i++) {
	if(peaks[i].prominence > minprominence) {
	  dummyp.push_back(peaks[i]);
	  }
	}
      return dummyp;
      }

	//else return the N most_prominent ones (might return less than N peaks if there aren't enough over minprominence)
    else {
      vector<ACFpeak> dummyp;
      for(unsigned int i = 0; i < peaks.size(); i++) {
	if(peaks[i].prominence > minprominence) {
	  dummyp.push_back(peaks[i]);
	  }
	}

      if(dummyp.size() > most_prominent) {
        vector<ACFpeak> newdummy;
	for(unsigned int j = 0; j < most_prominent; j++) {
	  int assign = -1;
	  double prom = 0.;
          for(unsigned int k = 0; k < dummyp.size(); k++) {
	    if(dummyp[k].prominence > prom){assign = k; prom = dummyp[k].prominence;}
	    }
	  if(assign != -1)newdummy.push_back(dummyp[assign]);
	  dummyp[assign].prominence = 0.;
	  }
	return newdummy;
        }


      else return dummyp;
      }
  }

}
