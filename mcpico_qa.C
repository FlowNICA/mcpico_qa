using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PxPyPzE4D<double>>;

void mcpico_qa(std::string infilelist, std::string outfile, std::string cm_energy="2.4", std::string Cms="true")
{
  const bool isCms = (bool)std::stoi(Cms);
  const double sNN = std::stod( cm_energy ); // in GeV
  const double M = 0.938; // in GeV/c^2
  const double T = sNN*sNN/(2.*M) - 2.*M;
  const double E = T + M;
  const double P = sqrt( E*E - M*M );
  const double Y_BEAM = 0.25 * log( (E + P) / (E - P) );
  const double GAMMA_CM = cosh(Y_BEAM);
  const double BETA_CM =  tanh(Y_BEAM);

  std::cout << "sqrtSnn = " << sNN << " GeV; T = " << T << "A GeV; Y_BEAM = " << Y_BEAM << std::endl;

  vector <RResultPtr<::TH1D >> hists;
  vector <RResultPtr<::TH2D >> hists2d;
  vector <RResultPtr<::TH3D >> hists3d;
  vector <RResultPtr<::THnD >> histsnd;
  vector <RResultPtr<::TProfile >> profs;
  vector <RResultPtr<::TProfile2D >> profs2d;

  // Centrality for 2.5 GeV
  std::vector<float> cent_bins;
  std::vector<float> cent_b;

  // Xe+W @ 2.87 GeV
  // cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 95.};
  // cent_b = { 1.43172, 2.87406, 4.05374, 5.0333, 5.86304, 6.58249, 7.22197, 7.80409, 8.34523, 8.8571, 9.34824, 9.8255, 10.2956, 10.7666, 11.2494, 11.7595, 12.3179, 12.9533, 13.7034};
  // Xe+Xe @ 2.87 GeV
  cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 95.};
  cent_b = { 1.39543, 2.70341, 3.76519, 4.64771, 5.40276, 6.06901, 6.67396, 7.23605, 7.76659, 8.27183, 8.75493, 9.21803, 9.66421, 10.0995, 10.5351, 10.9889, 11.488, 12.0706, 12.7879};

  TStopwatch timer;
  timer.Start();
  std::string treename = "mctree";
  TFileCollection collection( "collection", "", infilelist.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
    .Alias( "b", "bimp" )
    .Filter( "b<20.")
    .Define("cent", [cent_bins,cent_b](float b){
      float cent = -1.;
      if (cent_b.size() == 0) return (float)-1.;
      if (cent_bins.size() == 0) return (float)-1.;
      if (b <= cent_b.at(0))
        cent = cent_bins.at(0);
      for (int i=1; i<cent_b.size(); ++i){
        if (b < cent_b.at(i) && b >= cent_b.at(i-1))
          cent = cent_bins.at(i);
      }
      return cent;
    }, {"b"})
    .Define( "psi_rp", [](){
      std::random_device rnd_device;
      std::mt19937 generator(rnd_device());
      std::uniform_real_distribution<float> distribution(-1.*Pi(),Pi()); // distribution in range [-pi, pi]
      return distribution(generator);
    }, {} )
    .Define("particles", [GAMMA_CM,BETA_CM,isCms](RVec<float> &_px, RVec<float> &_py, RVec<float> &_pz, RVec<float> &_en, float _psi){
      RVec<fourVector> p;
      int Np = _px.size();
      for (int i=0; i<Np; i++){
        auto px_0 = _px.at(i);
        auto py_0 = _py.at(i);
        auto pt = sqrt(px_0*px_0+py_0*py_0);
        auto pz = (isCms) ? _pz.at(i) : GAMMA_CM*(_pz.at(i) - BETA_CM*_en.at(i));
        auto mass = sqrt( _en.at(i)*_en.at(i) - pt*pt - _pz.at(i)*_pz.at(i) );
        auto en = (isCms) ? _en.at(i) : GAMMA_CM*(_en.at(i) - BETA_CM*_pz.at(i)); //sqrt( mass*mass + pt*pt + pz*pz );
        auto phi = atan2(py_0, px_0);
        phi += _psi;
        auto px_1 = pt*cos(phi);
        auto py_1 = pt*sin(phi);
        p.push_back( {px_1, py_1, pz, en} );
      }
      return p;
    }, {"momx", "momy", "momz", "ene", "psi_rp"})
    .Define("target", [Y_BEAM,sNN](RVec<fourVector> _p){
      RVec<fourVector> result;
      for (auto &p:_p){
        if (abs(p.Rapidity()+Y_BEAM)<0.1 && p.Pt()<0.2 && abs(p.E() - sNN/2.)<0.2){ // target candidate
          result.push_back(p);
        }
      }
      return result;
    }, {"particles"})
    .Define("projectile", [Y_BEAM,sNN](RVec<fourVector> _p){
      RVec<fourVector> result;
      for (auto &p:_p){
        if (abs(p.Rapidity()-Y_BEAM)<0.1 && p.Pt()<0.2 && abs(p.E() - sNN/2.)<0.2){ // target candidate
          result.push_back(p);
        }
      }
      return result;
    }, {"particles"})
    .Define("phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"particles"})
    .Define("eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"particles"})
    .Define("pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; }, {"particles"})
    .Define("px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; }, {"particles"})
    .Define("py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; }, {"particles"})
    .Define("pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; }, {"particles"})
    .Define("en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },  {"particles"})
    .Define("y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"particles"})
    .Define("is_proton", "pdg==2212")
    .Define("is_pionM", "pdg==-211")
    .Define("is_pionP", "pdg==211")
    .Define("is_kaonM", "pdg==-321")
    .Define("is_kaonP", "pdg==321")
    .Define("protons", [](RVec<fourVector> _p, RVec<int> _pid){ RVec<fourVector> result; for (int i=0; i<_pid.size(); i++){auto pid = _pid.at(i); if (!pid) continue; result.push_back(_p.at(i));} return result; }, {"particles", "is_proton"})
    .Define("pionsM", [](RVec<fourVector> _p, RVec<int> _pid){ RVec<fourVector> result; for (int i=0; i<_pid.size(); i++){auto pid = _pid.at(i); if (!pid) continue; result.push_back(_p.at(i));} return result; }, {"particles", "is_pionM"})
    .Define("pionsP", [](RVec<fourVector> _p, RVec<int> _pid){ RVec<fourVector> result; for (int i=0; i<_pid.size(); i++){auto pid = _pid.at(i); if (!pid) continue; result.push_back(_p.at(i));} return result; }, {"particles", "is_pionP"})
    .Define("kaonsM", [](RVec<fourVector> _p, RVec<int> _pid){ RVec<fourVector> result; for (int i=0; i<_pid.size(); i++){auto pid = _pid.at(i); if (!pid) continue; result.push_back(_p.at(i));} return result; }, {"particles", "is_kaonM"})
    .Define("kaonsP", [](RVec<fourVector> _p, RVec<int> _pid){ RVec<fourVector> result; for (int i=0; i<_pid.size(); i++){auto pid = _pid.at(i); if (!pid) continue; result.push_back(_p.at(i));} return result; }, {"particles", "is_kaonP"})
    .Define("prot_phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"protons"})
    .Define("prot_eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"protons"})
    .Define("prot_pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; },      {"protons"})
    .Define("prot_px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; },      {"protons"})
    .Define("prot_py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; },      {"protons"})
    .Define("prot_pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; },      {"protons"})
    .Define("prot_en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },       {"protons"})
    .Define("prot_y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"protons"})
    .Define("pionM_phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"pionsM"})
    .Define("pionM_eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"pionsM"})
    .Define("pionM_pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; },      {"pionsM"})
    .Define("pionM_px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; },      {"pionsM"})
    .Define("pionM_py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; },      {"pionsM"})
    .Define("pionM_pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; },      {"pionsM"})
    .Define("pionM_en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },       {"pionsM"})
    .Define("pionM_y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"pionsM"})
    .Define("pionP_phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"pionsP"})
    .Define("pionP_eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"pionsP"})
    .Define("pionP_pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; },      {"pionsP"})
    .Define("pionP_px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; },      {"pionsP"})
    .Define("pionP_py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; },      {"pionsP"})
    .Define("pionP_pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; },      {"pionsP"})
    .Define("pionP_en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },       {"pionsP"})
    .Define("pionP_y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"pionsP"})
    .Define("kaonM_phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"kaonsM"})
    .Define("kaonM_eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"kaonsM"})
    .Define("kaonM_pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; },      {"kaonsM"})
    .Define("kaonM_px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; },      {"kaonsM"})
    .Define("kaonM_py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; },      {"kaonsM"})
    .Define("kaonM_pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; },      {"kaonsM"})
    .Define("kaonM_en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },       {"kaonsM"})
    .Define("kaonM_y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"kaonsM"})
    .Define("kaonP_phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"kaonsP"})
    .Define("kaonP_eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"kaonsP"})
    .Define("kaonP_pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; },      {"kaonsP"})
    .Define("kaonP_px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; },      {"kaonsP"})
    .Define("kaonP_py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; },      {"kaonsP"})
    .Define("kaonP_pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; },      {"kaonsP"})
    .Define("kaonP_en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },       {"kaonsP"})
    .Define("kaonP_y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"kaonsP"})
    .Define("protons_target", [Y_BEAM,sNN](RVec<fourVector> _p){
      RVec<fourVector> result;
      for (auto &p:_p)
        {
          if (abs(p.Rapidity()+Y_BEAM)<0.1 && p.Pt()<0.4 && abs(p.E() - sNN/2.)<0.2) // target candidate
            result.push_back(p);
        } return result;
    }, {"protons"})
    .Define("protons_projectile", [Y_BEAM,sNN](RVec<fourVector> _p){
      RVec<fourVector> result;
      for (auto &p:_p)
        {
          if (abs(p.Rapidity()-Y_BEAM)<0.1 && p.Pt()<0.4 && abs(p.E() - sNN/2.)<0.2) // projectile candidate
            result.push_back(p);
        } return result;
    }, {"protons"})
    .Define("prot_targ_phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"protons_target"})
    .Define("prot_targ_eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"protons_target"})
    .Define("prot_targ_pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; },      {"protons_target"})
    .Define("prot_targ_px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; },      {"protons_target"})
    .Define("prot_targ_py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; },      {"protons_target"})
    .Define("prot_targ_pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; },      {"protons_target"})
    .Define("prot_targ_en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },       {"protons_target"})
    .Define("prot_targ_y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"protons_target"})
    .Define("prot_proj_phi", [](RVec<fourVector> _p){ vector<float> phi; for (auto &p:_p){phi.push_back((float)p.Phi());} return phi; }, {"protons_projectile"})
    .Define("prot_proj_eta", [](RVec<fourVector> _p){ vector<float> eta; for (auto &p:_p){eta.push_back((float)p.Eta());} return eta; }, {"protons_projectile"})
    .Define("prot_proj_pT", [](RVec<fourVector> _p){ vector<float> pt; for (auto &p:_p){pt.push_back((float)p.Pt());} return pt; },      {"protons_projectile"})
    .Define("prot_proj_px", [](RVec<fourVector> _p){ vector<float> px; for (auto &p:_p){px.push_back((float)p.Px());} return px; },      {"protons_projectile"})
    .Define("prot_proj_py", [](RVec<fourVector> _p){ vector<float> py; for (auto &p:_p){py.push_back((float)p.Py());} return py; },      {"protons_projectile"})
    .Define("prot_proj_pz", [](RVec<fourVector> _p){ vector<float> pz; for (auto &p:_p){pz.push_back((float)p.Pz());} return pz; },      {"protons_projectile"})
    .Define("prot_proj_en", [](RVec<fourVector> _p){ vector<float> e;  for (auto &p:_p){e.push_back((float)p.E());}   return e; },       {"protons_projectile"})
    .Define("prot_proj_y", [Y_BEAM,isCms](RVec<fourVector> _p){
      vector<float> y;
      for (auto &p:_p){
        y.push_back((float)p.Rapidity());
        //if (isCms) y.push_back((float)p.Rapidity());
        //else y.push_back((float)(p.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"protons_projectile"})
    .Define("N_prot_targ", [](RVec<fourVector> _p){ return (float)_p.size(); }, {"protons_target"})
    .Define("N_prot_proj", [](RVec<fourVector> _p){ return (float)_p.size(); }, {"protons_projectile"})
    .Define("N_prot_targ_proj_ratio", [](float _ntarg, float _nproj){ return (float)(_ntarg/_nproj); }, {"N_prot_targ", "N_prot_proj"})
    ;

  TFile fOut(outfile.c_str(),"recreate");

  // Make lists of histograms for QA
  hists.push_back(dd.Histo1D({"h1_b","Impact parameter;b (fm)",200,0.,20.}, "b"));
  
  hists.push_back(dd.Histo1D({"h1_hadr_pT",";p_{T} (GeV/c)",500,0.,5.}, "pT"));
  hists.push_back(dd.Histo1D({"h1_prot_pT",";p_{T} (GeV/c)",500,0.,5.}, "prot_pT"));
  hists.push_back(dd.Histo1D({"h1_pionM_pT",";p_{T} (GeV/c)",500,0.,5.}, "pionM_pT"));
  hists.push_back(dd.Histo1D({"h1_pionP_pT",";p_{T} (GeV/c)",500,0.,5.}, "pionP_pT"));

  hists.push_back(dd.Histo1D({"h1_hadr_px",";p_{x} (GeV/c)",1000,-5.,5.}, "px"));
  hists.push_back(dd.Histo1D({"h1_prot_px",";p_{x} (GeV/c)",1000,-5.,5.}, "prot_px"));
  hists.push_back(dd.Histo1D({"h1_pionM_px",";p_{x} (GeV/c)",1000,-5.,5.}, "pionM_px"));
  hists.push_back(dd.Histo1D({"h1_pionP_px",";p_{x} (GeV/c)",1000,-5.,5.}, "pionP_px"));

  hists.push_back(dd.Histo1D({"h1_hadr_py",";p_{y} (GeV/c)",1000,-5.,5.}, "py"));
  hists.push_back(dd.Histo1D({"h1_prot_py",";p_{y} (GeV/c)",1000,-5.,5.}, "prot_py"));
  hists.push_back(dd.Histo1D({"h1_pionM_py",";p_{y} (GeV/c)",1000,-5.,5.}, "pionM_py"));
  hists.push_back(dd.Histo1D({"h1_pionP_py",";p_{y} (GeV/c)",1000,-5.,5.}, "pionP_py"));

  hists.push_back(dd.Histo1D({"h1_hadr_pz",";p_{z} (GeV/c)",1000,-5.,5.}, "pz"));
  hists.push_back(dd.Histo1D({"h1_prot_pz",";p_{z} (GeV/c)",1000,-5.,5.}, "prot_pz"));
  hists.push_back(dd.Histo1D({"h1_pionM_pz",";p_{z} (GeV/c)",1000,-5.,5.}, "pionM_pz"));
  hists.push_back(dd.Histo1D({"h1_pionP_pz",";p_{z} (GeV/c)",1000,-5.,5.}, "pionP_pz"));

  hists.push_back(dd.Histo1D({"h1_hadr_y",";y",1000,-5.,5.}, "y"));
  hists.push_back(dd.Histo1D({"h1_prot_y",";y",1000,-5.,5.}, "prot_y"));
  hists.push_back(dd.Histo1D({"h1_pionM_y",";y",1000,-5.,5.}, "pionM_y"));
  hists.push_back(dd.Histo1D({"h1_pionP_y",";y",1000,-5.,5.}, "pionP_y"));

  hists.push_back(dd.Histo1D({"h1_hadr_en",";E (GeV)",500,0.,5.}, "en"));
  hists.push_back(dd.Histo1D({"h1_prot_en",";E (GeV)",500,0.,5.}, "prot_en"));
  hists.push_back(dd.Histo1D({"h1_pionM_en",";E (GeV)",500,0.,5.}, "pionM_en"));
  hists.push_back(dd.Histo1D({"h1_pionP_en",";E (GeV)",500,0.,5.}, "pionP_en"));

  hists.push_back(dd.Histo1D({"h1_prot_targ_N","Number of target protons;N_{targ}",500,0.,500.},     "N_prot_targ"));
  hists.push_back(dd.Histo1D({"h1_prot_proj_N","Number of projectile protons;N_{proj}",500,0.,500.}, "N_prot_proj"));
  hists.push_back(dd.Histo1D({"h1_prot_targ_proj_ratio","N_{targ}/N_{proj} pf protons;N_{targ}/N_{proj}",500,0.,50.}, "N_prot_targ_proj_ratio"));
  
  hists2d.push_back(dd.Histo2D({"h2_hadr_pTY","p_{T}-Y acceptance;y;p_{T} (GeV/c)",500,-5.,5.,500,0.,5.}, "y", "pT"));
  hists2d.push_back(dd.Histo2D({"h2_prot_pTY","p_{T}-Y acceptance (p);y;p_{T} (GeV/c)",500,-5.,5.,500,0.,5.}, "prot_y", "prot_pT"));
  hists2d.push_back(dd.Histo2D({"h2_pionM_pTY","p_{T}-Y acceptance (#pi^{-});y;p_{T} (GeV/c)",500,-5.,5.,500,0.,5.}, "pionM_y", "pionM_pT"));
  hists2d.push_back(dd.Histo2D({"h2_pionP_pTY","p_{T}-Y acceptance (#pi^{+});y;p_{T} (GeV/c)",500,-5.,5.,500,0.,5.}, "pionP_y", "pionP_pT"));
  hists2d.push_back(dd.Histo2D({"h2_kaonM_pTY","p_{T}-Y acceptance (#pi^{-});y;p_{T} (GeV/c)",500,-5.,5.,500,0.,5.}, "kaonM_y", "kaonM_pT"));
  hists2d.push_back(dd.Histo2D({"h2_kaonP_pTY","p_{T}-Y acceptance (#pi^{+});y;p_{T} (GeV/c)",500,-5.,5.,500,0.,5.}, "kaonP_y", "kaonP_pT"));

  hists2d.push_back(dd.Histo2D({"h2_hadr_pxY","p_{x}-Y acceptance;y;p_{x} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "y", "px"));
  hists2d.push_back(dd.Histo2D({"h2_prot_pxY","p_{x}-Y acceptance (p);y;p_{x} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "prot_y", "prot_px"));
  hists2d.push_back(dd.Histo2D({"h2_pionM_pxY","p_{x}-Y acceptance (#pi^{-});y;p_{x} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "pionM_y", "pionM_px"));
  hists2d.push_back(dd.Histo2D({"h2_pionP_pxY","p_{x}-Y acceptance (#pi^{+});y;p_{x} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "pionP_y", "pionP_px"));
  hists2d.push_back(dd.Histo2D({"h2_kaonM_pxY","p_{x}-Y acceptance (#pi^{-});y;p_{x} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "kaonM_y", "kaonM_px"));
  hists2d.push_back(dd.Histo2D({"h2_kaonP_pxY","p_{x}-Y acceptance (#pi^{+});y;p_{x} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "kaonP_y", "kaonP_px"));

  hists2d.push_back(dd.Histo2D({"h2_hadr_pyY","p_{y}-Y acceptance;y;p_{y} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "y", "py"));
  hists2d.push_back(dd.Histo2D({"h2_prot_pyY","p_{y}-Y acceptance (p);y;p_{y} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "prot_y", "prot_py"));
  hists2d.push_back(dd.Histo2D({"h2_pionM_pyY","p_{y}-Y acceptance (#pi^{-});y;p_{y} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "pionM_y", "pionM_py"));
  hists2d.push_back(dd.Histo2D({"h2_pionP_pyY","p_{y}-Y acceptance (#pi^{+});y;p_{y} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "pionP_y", "pionP_py"));
  hists2d.push_back(dd.Histo2D({"h2_kaonM_pyY","p_{y}-Y acceptance (#pi^{-});y;p_{y} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "kaonM_y", "kaonM_py"));
  hists2d.push_back(dd.Histo2D({"h2_kaonP_pyY","p_{y}-Y acceptance (#pi^{+});y;p_{y} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "kaonP_y", "kaonP_py"));

  hists2d.push_back(dd.Histo2D({"h2_hadr_pzY","p_{z}-Y acceptance;y;p_{z} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "y", "pz"));
  hists2d.push_back(dd.Histo2D({"h2_prot_pzY","p_{z}-Y acceptance (p);y;p_{z} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "prot_y", "prot_pz"));
  hists2d.push_back(dd.Histo2D({"h2_pionM_pzY","p_{z}-Y acceptance (#pi^{-});y;p_{z} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "pionM_y", "pionM_pz"));
  hists2d.push_back(dd.Histo2D({"h2_pionP_pzY","p_{z}-Y acceptance (#pi^{+});y;p_{z} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "pionP_y", "pionP_pz"));
  hists2d.push_back(dd.Histo2D({"h2_kaonM_pzY","p_{z}-Y acceptance (#pi^{-});y;p_{z} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "kaonM_y", "kaonM_pz"));
  hists2d.push_back(dd.Histo2D({"h2_kaonP_pzY","p_{z}-Y acceptance (#pi^{+});y;p_{z} (GeV/c)",500,-5.,5.,1000,-5.,5.}, "kaonP_y", "kaonP_pz"));

  hists2d.push_back(dd.Histo2D({"h2_hadr_enY","E-Y acceptance;y;E (GeV)",500,-5.,5.,500,0.,5.}, "y", "en"));
  hists2d.push_back(dd.Histo2D({"h2_prot_enY","E-Y acceptance (p);y;E (GeV)",500,-5.,5.,500,0.,5.}, "prot_y", "prot_en"));
  hists2d.push_back(dd.Histo2D({"h2_pionM_enY","E-Y acceptance (#pi^{-});y;E (GeV)",500,-5.,5.,500,0.,5.}, "pionM_y", "pionM_en"));
  hists2d.push_back(dd.Histo2D({"h2_pionP_enY","E-Y acceptance (#pi^{+});y;E (GeV)",500,-5.,5.,500,0.,5.}, "pionP_y", "pionP_en"));
  hists2d.push_back(dd.Histo2D({"h2_kaonM_enY","E-Y acceptance (#pi^{-});y;E (GeV)",500,-5.,5.,500,0.,5.}, "kaonM_y", "kaonM_en"));
  hists2d.push_back(dd.Histo2D({"h2_kaonP_enY","E-Y acceptance (#pi^{+});y;E (GeV)",500,-5.,5.,500,0.,5.}, "kaonP_y", "kaonP_en"));

  fOut.cd();
  for (auto& hist:hists)
    hist->Write();
  for (auto& hist:hists2d)
    hist->Write();
  for (auto& hist:hists3d)
    hist->Write();
  for (auto& hist:profs)
    hist->Write();
  for (auto& hist:profs2d)
    hist->Write();
  fOut.Close();

  std::cout << std::endl;
  timer.Stop();
  timer.Print();
}
