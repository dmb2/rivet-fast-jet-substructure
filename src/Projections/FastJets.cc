// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

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

  ///Q===\Sum_{i\in J} q_i*p_{Ti}^k/p_{TJ} 
  double FastJets::JetCharge(const fastjet::PseudoJet &j, const double k=0.5, const double ptmin=-1*GeV) const {
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
    //finally declare the filter and apply it to the jet
    fastjet::Pruner pruner(fastjet::JetDefinition(setJetAlgorithm(subjet_def), 1), zcut, Rcut_factor);
    return pruner(jet);
  }
  PseudoJets FastJets::GetAxes( int n_jets, 
		PseudoJets inputJets, JetAlgName subjet_def, double subR) const {

    assert(clusterSeq());
    //sanity check
    if ((signed) inputJets.size() < n_jets) { // Dirty cast...
	std::cout << "Not enough input particles." << endl;
	return inputJets;
	}
    //get subjets, return
    fastjet::ClusterSequence sub_clust_seq(inputJets, fastjet::JetDefinition(setJetAlgorithm(subjet_def), subR));
    return sub_clust_seq.exclusive_jets(n_jets);
  }

  double FastJets::TauValue(double beta, double jet_rad, 
	PseudoJets particles, PseudoJets axes) const {

    double tauNum = 0.0;
    double tauDen = 0.0;
   
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

  vector<fastjet::PseudoJet> FastJets::UpdateAxes(double beta,
	PseudoJets particles, PseudoJets axes) const {

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
    
    double PI = 3.14159265359;
    double deltaR2, distphi = 0.0;
    vector<double> ynom, phinom, den;
    ynom.resize(axes.size());
    phinom.resize(axes.size());
    den.resize(axes.size());

    for (unsigned int i = 0; i < particles.size(); i++) {
      distphi = particles[i].phi() - axes[belongsto[i]].phi();
      deltaR2 = particles[i].squared_distance(axes[belongsto[i]]);
      
      if (abs(distphi) <= PI) phinom.at(belongsto[i]) += particles[i].perp() * particles[i].phi() * pow(deltaR2, (beta-2)/2);
      else if ( distphi > PI) phinom.at(belongsto[i]) += particles[i].perp() * (-2 * PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);
      else if ( distphi < PI) phinom.at(belongsto[i]) += particles[i].perp() * (+2 * PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);

      ynom.at(belongsto[i]) += particles[i].perp() * particles[i].rap() * pow(deltaR2, (beta-2)/2);
      
      den.at(belongsto[i]) += particles[i].perp() * pow(deltaR2, (beta-2)/2);
    }
      
    // declare new axes
    vector<fastjet::PseudoJet> newaxes;
    newaxes.resize(axes.size());
    for (unsigned int j = 0; j < axes.size(); j++) {

      if (den[j] == 0) newaxes.at(j) = axes[j];
      else {
	newaxes.at(j).reset_momentum_PtYPhiM(axes[j].perp(), ynom[j] / den[j], fmod( 2*PI + (phinom[j] / den[j]), 2*PI ), axes[j].perp()/2);
      }
    }
    return newaxes;
  }
}
