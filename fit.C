int fit(int j = 20){

std::vector<int> runlist = {984, 987, 988, 990, 991, 992, 993, 995, 996, 1000, 1001, 1003, 1004, 1005, 1006, 1007, 1008, 1010, 1012, 1014, 1015, 1023, 1024, 1025, 1027, 1028};
//std::vector<int> runlist = {984, 987, 988};

auto planeStr = "plane"+std::to_string(j);
TH1D* sigma_histo = new TH1D(("res_"+planeStr).c_str(),("res_"+planeStr).c_str(),50,0,10);

TH1D* width_X = new TH1D("widthX","widthX",40,240,260);
TH1D* mean_X = new TH1D("meanX","meanX",10,-10,10);

TH1D* width_Y = new TH1D("widthY","widthY",40,40,60);
TH1D* mean_Y = new TH1D("meanY","meanY",10,-10,10);

TCanvas* canvas = new TCanvas("ProjCan", "Residual Fits", 3000, 750);
canvas->Divide(runlist.size(),3);
size_t i = 1;

for(auto run: runlist) {
		
		std::stringstream ss;
		ss << std::setw(4) << std::setfill('0') << run;

		auto file = "output/histograms/run00"+ss.str()+"-fitter-histo.root";

		TFile* _file0 = TFile::Open(file.c_str());
		TDirectoryFile* f1 = (TDirectoryFile*)_file0->Get("DafFitter");
		TH1D* hX = static_cast<TH1D*>(f1->Get(("pl"+std::to_string(j)+"_residualX").c_str()));
		TH1D* hY = static_cast<TH1D*>(f1->Get(("pl"+std::to_string(j)+"_residualY").c_str()));
		TFormula* form = new TFormula("box", "[0]/2*(1-erf((x-[1])/(sqrt(2)*[2])))*(x>0)+[0]/2*(1-erf(([3]-x)/(sqrt(2)*[2])))*(x<0)");

		TF1* func = new TF1("myFit","box",-1,1);
		func->SetParNames("max","right","sigma","left");
		func->SetParameters(450,0.1,0.03,-0.1);

		canvas->cd(i);
		TFitResultPtr fitX = hX->Fit("myFit","QS");

		auto paramsX = fitX->Parameters();
		auto maxX = paramsX.at(0);
		auto rightX = paramsX.at(1);
		auto sigmaX = paramsX.at(2);
		auto leftX = paramsX.at(3);
		auto nDofX = fitX->Ndf();
		auto Chi2X = fitX->Chi2();

		sigma_histo->Fill(sigmaX*1000);
		width_X->Fill((rightX-leftX)*1000);
		mean_X->Fill(1000*(rightX+leftX)/2.0);		
		//hX->Draw();

		std::cout << "Fit " << planeStr << " in X:\n" << "left|right: " << leftX << "|" << rightX << " is distance: " << rightX-leftX << "\n" << "Sigma: " << sigmaX << std::endl;
		std::cout << "Fit " << planeStr << " in X:\n" << "Chi2/Ndof = " << Chi2X << "/" << nDofX << " = " << Chi2X/nDofX << std::endl;

		TF1* func2 = new TF1("myFit2","box",-1,1);
		func2->SetParNames("max","right","sigma","left");
		func2->SetParameters(1800,0.05,sigmaX,-0.05);
		func2->SetParLimits(1,0,0.2);
		func2->SetParLimits(3,-0.2,0);
		func2->FixParameter(2,sigmaX);

		canvas->cd(i+runlist.size());
		TFitResultPtr fitY = hY->Fit("myFit2","QBS");
		
		//hY->Draw();

		auto paramsY = fitY->Parameters();
		auto maxY = paramsY.at(0);
		auto rightY = paramsY.at(1);
		auto sigmaY = paramsY.at(2);
		auto leftY = paramsY.at(3);
		auto nDofY = fitY->Ndf();
		auto Chi2Y = fitY->Chi2();

		width_Y->Fill((rightY-leftY)*1000);
		mean_Y->Fill(1000*(rightY+leftY)/2.0);		

		std::cout << "Fit " << planeStr << " in Y:\n" << "left|right: " << leftY << "|" << rightY << " is distance: " << rightY-leftY << "\n" << "Sigma: " << sigmaY << std::endl;
		std::cout << "Fit " << planeStr << " in Y:\n" << "Chi2/Ndof = " << Chi2Y << "/" << nDofY << " = " << Chi2Y/nDofY << std::endl;

		i++;
		}

auto name = "output_plane_"+std::to_string(j)+".root"; 
TFile outputFile(name.c_str(), "RECREATE");
gFile = &outputFile;

canvas->cd(2*runlist.size()+1);
sigma_histo->Draw();
sigma_histo->Write();

canvas->cd(2*runlist.size()+2);
width_X->Draw();
width_X->Write();
canvas->cd(2*runlist.size()+3);
mean_X->Draw();
mean_X->Write();

canvas->cd(2*runlist.size()+4);
width_Y->Draw();
width_Y->Write();
canvas->cd(2*runlist.size()+5);
mean_Y->Write();

outputFile.Close();

return 0;
}
