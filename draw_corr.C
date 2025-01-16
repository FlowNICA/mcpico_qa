#include "ris2.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void draw_corr()
{
  gROOT->Macro("~/Soft/ris2/macro/style.cc");

  // std::vector<std::string> files = {"~/Documents/ModelData/QnToolsData/calc_urqmd_xexe_2.87gev_cascade.root",
  //                                   "~/Documents/ModelData/QnToolsData/calc_urqmd_xexe_2.5agev_cascade.root"};

  std::vector<std::string> files = {"~/Documents/ModelData/QnToolsData/calc_urqmd_xexe_2.87gev_mf.root",
                                    "~/Documents/ModelData/QnToolsData/calc_urqmd_xexe_2.5agev_mf.root"};

  std::pair<float, float> v1_pt_cut = {1.0, 1.5};
  std::pair<float, float> v1_y_cut = {-0.5, -0.1};
  std::pair<float, float> v2_pt_cut = {1.0, 1.5};
  std::pair<float, float> v2_y_cut = {-0.5, -0.1};
  std::pair<float, float> cent_cut = {10., 40.};

  auto mc_v1_proton = ris2::Bunch<ris2::Correlation>{}
    .AddToBunch(std::string{"CMS"},
    files.at(0),
    std::vector<std::string>{
        "proton/v1.tru_proton_PLAIN.psi_rp_PLAIN.x1x1cent",
        "proton/v1.tru_proton_PLAIN.psi_rp_PLAIN.y1y1cent",
    },
    std::vector<double>{1., 1.})
    .AddToBunch(
        std::string{"LAB"},
        files.at(1),
        std::vector<std::string>{
            "proton/v1.tru_proton_PLAIN.psi_rp_PLAIN.x1x1cent",
            "proton/v1.tru_proton_PLAIN.psi_rp_PLAIN.y1y1cent",
        },
        std::vector<double>{1., 1.});
  mc_v1_proton.GetPalette().SetPalette(std::vector{
      ris2::Style().SetColor(kRed).SetMarker(-1),
      ris2::Style().SetColor(kBlue).SetMarker(-1),
  });

/*  auto mc_v1_pionP = ris2::Bunch<ris2::Correlation>{}
  .AddToBunch(std::string{"CMS"},
    files.at(0),
    std::vector<std::string>{
        "pionP/v1.tru_pionP_PLAIN.psi_rp_PLAIN.x1x1cent",
        "pionP/v1.tru_pionP_PLAIN.psi_rp_PLAIN.y1y1cent",
    },
    std::vector<double>{1., 1.})
    .AddToBunch(
        std::string{"LAB"},
        files.at(1),
        std::vector<std::string>{
            "pionP/v1.tru_pionP_PLAIN.psi_rp_PLAIN.x1x1cent",
            "pionP/v1.tru_pionP_PLAIN.psi_rp_PLAIN.y1y1cent",
        },
        std::vector<double>{1., 1.});
  mc_v1_pionP.GetPalette().SetPalette(std::vector{
      ris2::Style().SetColor(kRed).SetMarker(-1),
      ris2::Style().SetColor(kBlue).SetMarker(-1),
  });

  auto mc_v1_pionM = ris2::Bunch<ris2::Correlation>{}
  .AddToBunch(std::string{"CMS"},
    files.at(0),
    std::vector<std::string>{
        "pionM/v1.tru_pionM_PLAIN.psi_rp_PLAIN.x1x1cent",
        "pionM/v1.tru_pionM_PLAIN.psi_rp_PLAIN.y1y1cent",
    },
    std::vector<double>{1., 1.})
    .AddToBunch(
        std::string{"LAB"},
        files.at(1),
        std::vector<std::string>{
            "pionM/v1.tru_pionM_PLAIN.psi_rp_PLAIN.x1x1cent",
            "pionM/v1.tru_pionM_PLAIN.psi_rp_PLAIN.y1y1cent",
        },
        std::vector<double>{1., 1.});
  mc_v1_pionM.GetPalette().SetPalette(std::vector{
      ris2::Style().SetColor(kRed).SetMarker(-1),
      ris2::Style().SetColor(kBlue).SetMarker(-1),
  });*/

  auto mc_v2_proton = ris2::Bunch<ris2::Correlation>{}
  .AddToBunch(std::string{"CMS"},
    files.at(0),
    std::vector<std::string>{
        "proton/v2.tru_proton_PLAIN.psi_rp_PLAIN.x2x2cent",
        "proton/v2.tru_proton_PLAIN.psi_rp_PLAIN.y2y2cent",
    },
    std::vector<double>{1., 1.})
    .AddToBunch(
        std::string{"LAB"},
        files.at(1),
        std::vector<std::string>{
            "proton/v2.tru_proton_PLAIN.psi_rp_PLAIN.x2x2cent",
            "proton/v2.tru_proton_PLAIN.psi_rp_PLAIN.y2y2cent",
        },
        std::vector<double>{1., 1.});
  mc_v2_proton.GetPalette().SetPalette(std::vector{
      ris2::Style().SetColor(kRed).SetMarker(-1),
      ris2::Style().SetColor(kBlue).SetMarker(-1),
  });
/*
  auto mc_v2_pionP = ris2::Bunch<ris2::Correlation>{}
  .AddToBunch(std::string{"CMS"},
    files.at(0),
    std::vector<std::string>{
        "pionP/v2.tru_pionP_PLAIN.psi_rp_PLAIN.x2x2cent",
        "pionP/v2.tru_pionP_PLAIN.psi_rp_PLAIN.y2y2cent",
    },
    std::vector<double>{1., 1.})
    .AddToBunch(
        std::string{"LAB"},
        files.at(1),
        std::vector<std::string>{
            "pionP/v2.tru_pionP_PLAIN.psi_rp_PLAIN.x2x2cent",
            "pionP/v2.tru_pionP_PLAIN.psi_rp_PLAIN.y2y2cent",
        },
        std::vector<double>{1., 1.});
  mc_v2_pionP.GetPalette().SetPalette(std::vector{
      ris2::Style().SetColor(kRed).SetMarker(-1),
      ris2::Style().SetColor(kBlue).SetMarker(-1),
  });

  auto mc_v2_pionM = ris2::Bunch<ris2::Correlation>{}
  .AddToBunch(std::string{"CMS"},
    files.at(0),
    std::vector<std::string>{
        "pionM/v2.tru_pionM_PLAIN.psi_rp_PLAIN.x2x2cent",
        "pionM/v2.tru_pionM_PLAIN.psi_rp_PLAIN.y2y2cent",
    },
    std::vector<double>{1., 1.})
    .AddToBunch(
        std::string{"LAB"},
        files.at(1),
        std::vector<std::string>{
            "pionM/v2.tru_pionM_PLAIN.psi_rp_PLAIN.x2x2cent",
            "pionM/v2.tru_pionM_PLAIN.psi_rp_PLAIN.y2y2cent",
        },
        std::vector<double>{1., 1.});
  mc_v2_pionM.GetPalette().SetPalette(std::vector{
      ris2::Style().SetColor(kRed).SetMarker(-1),
      ris2::Style().SetColor(kBlue).SetMarker(-1),
  });
*/
  // Making projections
  auto mc_v1_proton_y = mc_v1_proton;
  mc_v1_proton_y.Perform( [cent_cut, v1_pt_cut, v1_y_cut]( auto& obj ){
    obj.Rebin( std::vector<Qn::AxisD>{
        {"cent", 1, cent_cut.first, cent_cut.second},
        {"pT", 1, v1_pt_cut.first, v1_pt_cut.second},
        {"y", 15, -1.5, 1.5 },
      }).Project(std::vector<Qn::AxisD>{{"y", 15, -1.5, 1.5 }});
  });
  auto mc_v1_proton_pT = mc_v1_proton;
  mc_v1_proton_pT.Perform( [cent_cut, v1_pt_cut, v1_y_cut]( auto& obj ){
    obj.Rebin( std::vector<Qn::AxisD>{
        {"cent", 1, cent_cut.first, cent_cut.second},
        {"y", 1, v1_y_cut.first, v1_y_cut.second},
        {"pT", {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 2.} },
      }).Project(std::vector<Qn::AxisD>{{"pT", {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 2.} }});
  });
  auto mc_v1_proton_cent = mc_v1_proton;
  mc_v1_proton_cent.Perform( [cent_cut, v1_pt_cut, v1_y_cut]( auto& obj ){
    obj.Rebin( std::vector<Qn::AxisD>{
        {"y", 1, v1_y_cut.first, v1_y_cut.second},
        {"pT", 1, v1_pt_cut.first, v1_pt_cut.second},
        {"cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.}},
      }).Project(std::vector<Qn::AxisD>{{"cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.} }});
  });

  auto mc_v2_proton_y = mc_v2_proton;
  mc_v2_proton_y.Perform( [cent_cut, v2_pt_cut, v2_y_cut]( auto& obj ){
    obj.Rebin( std::vector<Qn::AxisD>{
        {"cent", 1, cent_cut.first, cent_cut.second},
        {"pT", 1, v2_pt_cut.first, v2_pt_cut.second},
        {"y", 15, -1.5, 1.5 },
      }).Project(std::vector<Qn::AxisD>{{"y", 15, -1.5, 1.5 }});
  });
  auto mc_v2_proton_pT = mc_v2_proton;
  mc_v2_proton_pT.Perform( [cent_cut, v2_pt_cut, v2_y_cut]( auto& obj ){
    obj.Rebin( std::vector<Qn::AxisD>{
        {"cent", 1, cent_cut.first, cent_cut.second},
        {"y", 1, v2_y_cut.first, v2_y_cut.second},
        {"pT", {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 2.} },
      }).Project(std::vector<Qn::AxisD>{{"pT", {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 2.} }});
  });
  auto mc_v2_proton_cent = mc_v2_proton;
  mc_v2_proton_cent.Perform( [cent_cut, v2_pt_cut, v2_y_cut]( auto& obj ){
    obj.Rebin( std::vector<Qn::AxisD>{
        {"y", 1, v2_y_cut.first, v2_y_cut.second},
        {"pT", 1, v2_pt_cut.first, v2_pt_cut.second},
        {"cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.}},
      }).Project(std::vector<Qn::AxisD>{{"cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.} }});
  });

  // Making plots
  auto plot_v1_proton_y = ris2::Plot( {1000, 1100} );
  auto leg_v1_proton_y = mc_v1_proton_y.MakeLegend( {0.25, 0.55, 0.55, 0.85} );
  plot_v1_proton_y.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(ris2::Axis().SetTitle("y_{cm}").SetLo(-1.5).SetHi(1.5))
    .SetYAxis(ris2::Axis().SetTitle("v_{1}").SetLo(-0.75).SetHi(0.75))
    .AddToPlot( *mc_v1_proton_y )
    .AddLegend( leg_v1_proton_y )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.9})
             .SetSize(0.04)
             .SetText(std::string{"UrQMD Xe+Xe, T=2.5A GeV"}) )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.86})
             .SetSize(0.035)
             .SetText(std::string{Form("p; %1.0f-%1.0f%s; %1.1f<p_{T}<%1.1f (GeV/c)", cent_cut.first, cent_cut.second, "%", v1_pt_cut.first, v1_pt_cut.second)}) )
    .AddFunction(new TF1("hline","pol0",-1.5,1.5))
  ;
  plot_v1_proton_y.Print( "~/Documents/ModelData/mcpico_qa/draw_corr/v1_proton_y.png" );
  auto plot_v1_proton_pT = ris2::Plot( {1000, 1100} );
  auto leg_v1_proton_pT = mc_v1_proton_pT.MakeLegend( {0.7, 0.75, 0.99, 0.95} );
  plot_v1_proton_pT.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(ris2::Axis().SetTitle("p_{T}, GeV/c").SetLo(0.).SetHi(2.))
    .SetYAxis(ris2::Axis().SetTitle("v_{1}").SetLo(-0.4).SetHi(0.02))
    // .SetYAxis(ris2::Axis().SetTitle("v_{1}").SetLo(-0.15).SetHi(0.02))
    .AddToPlot( *mc_v1_proton_pT )
    .AddLegend( leg_v1_proton_pT )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.9})
             .SetSize(0.04)
             .SetText(std::string{"UrQMD Xe+Xe, T=2.5A GeV"}) )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.86})
             .SetSize(0.035)
             .SetText(std::string{Form("p; %1.0f-%1.0f%s; %1.1f<y_{CM}<%1.1f", cent_cut.first, cent_cut.second, "%", v1_y_cut.first, v1_y_cut.second)}) )
  ;
  plot_v1_proton_pT.Print( "~/Documents/ModelData/mcpico_qa/draw_corr/v1_proton_pt.png" );
  auto plot_v1_proton_cent = ris2::Plot( {1000, 1100} );
  auto leg_v1_proton_cent = mc_v1_proton_cent.MakeLegend( {0.7, 0.68, 0.99, 0.88} );
  plot_v1_proton_cent.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(ris2::Axis().SetTitle("Centrality, %").SetLo(0.).SetHi(100.))
    .SetYAxis(ris2::Axis().SetTitle("v_{1}").SetLo(-0.35).SetHi(0.28))
    // .SetYAxis(ris2::Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.1))
    .AddToPlot( *mc_v1_proton_cent )
    .AddLegend( leg_v1_proton_cent )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.9})
             .SetSize(0.04)
             .SetText(std::string{"UrQMD Xe+Xe, T=2.5A GeV"}) )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.86})
             .SetSize(0.035)
             .SetText(std::string{Form("p; %1.1f<y_{CM}<%1.1f; %1.1f<p_{T}<%1.1f (GeV/c)", v1_y_cut.first, v1_y_cut.second, v1_pt_cut.first, v1_pt_cut.second)}) )
  ;
  plot_v1_proton_cent.Print( "~/Documents/ModelData/mcpico_qa/draw_corr/v1_proton_cent.png" );

  auto plot_v2_proton_y = ris2::Plot( {1000, 1100} );
  auto leg_v2_proton_y = mc_v2_proton_y.MakeLegend( {0.38, 0.55, 0.68, 0.85} );
  plot_v2_proton_y.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(ris2::Axis().SetTitle("y_{cm}").SetLo(-1.5).SetHi(1.5))
    .SetYAxis(ris2::Axis().SetTitle("v_{2}").SetLo(-0.25).SetHi(0.3))
    // .SetYAxis(ris2::Axis().SetTitle("v_{2}").SetLo(-0.01).SetHi(0.15))
    .AddToPlot( *mc_v2_proton_y )
    .AddLegend( leg_v2_proton_y )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.9})
             .SetSize(0.04)
             .SetText(std::string{"UrQMD Xe+Xe, T=2.5A GeV"}) )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.86})
             .SetSize(0.035)
             .SetText(std::string{Form("p; %1.0f-%1.0f%s; %1.1f<p_{T}<%1.1f (GeV/c)", cent_cut.first, cent_cut.second, "%", v2_pt_cut.first, v2_pt_cut.second)}) )
    .AddFunction(new TF1("hline","pol0",-1.5,1.5))
  ;
  plot_v2_proton_y.Print( "~/Documents/ModelData/mcpico_qa/draw_corr/v2_proton_y.png" );
  auto plot_v2_proton_pT = ris2::Plot( {1000, 1100} );
  auto leg_v2_proton_pT = mc_v2_proton_pT.MakeLegend( {0.7, 0.75, 0.99, 0.95} );
  plot_v2_proton_pT.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(ris2::Axis().SetTitle("p_{T}, GeV/c").SetLo(0.).SetHi(2.))
    .SetYAxis(ris2::Axis().SetTitle("v_{2}").SetLo(-0.06).SetHi(0.016))
    // .SetYAxis(ris2::Axis().SetTitle("v_{2}").SetLo(-0.01).SetHi(0.06))
    .AddToPlot( *mc_v2_proton_pT )
    .AddLegend( leg_v2_proton_pT )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.9})
             .SetSize(0.04)
             .SetText(std::string{"UrQMD Xe+Xe, T=2.5A GeV"}) )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.86})
             .SetSize(0.035)
             .SetText(std::string{Form("p; %1.0f-%1.0f%s; %1.1f<y_{CM}<%1.1f", cent_cut.first, cent_cut.second, "%", v2_y_cut.first, v2_y_cut.second)}) )
  ;
  plot_v2_proton_pT.Print( "~/Documents/ModelData/mcpico_qa/draw_corr/v2_proton_pt.png" );
  auto plot_v2_proton_cent = ris2::Plot( {1000, 1100} );
  auto leg_v2_proton_cent = mc_v2_proton_cent.MakeLegend( {0.7, 0.68, 0.99, 0.88} );
  plot_v2_proton_cent.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(ris2::Axis().SetTitle("Centrality, %").SetLo(0.).SetHi(100.))
    .SetYAxis(ris2::Axis().SetTitle("v_{2}").SetLo(-0.09).SetHi(0.01))
    // .SetYAxis(ris2::Axis().SetTitle("v_{2}").SetLo(-0.09).SetHi(0.038))
    .AddToPlot( *mc_v2_proton_cent )
    .AddLegend( leg_v2_proton_cent )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.9})
             .SetSize(0.04)
             .SetText(std::string{"UrQMD Xe+Xe, T=2.5A GeV"}) )
    .AddText(
             ris2::Text().SetPosition(std::array<double, 2>{0.25, 0.86})
             .SetSize(0.035)
             .SetText(std::string{Form("p; %1.1f<y_{CM}<%1.1f; %1.1f<p_{T}<%1.1f (GeV/c)", v2_y_cut.first, v2_y_cut.second, v2_pt_cut.first, v2_pt_cut.second)}) )
  ;
  plot_v2_proton_cent.Print( "~/Documents/ModelData/mcpico_qa/draw_corr/v2_proton_cent.png" );
  
}
