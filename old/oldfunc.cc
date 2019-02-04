// Functions not used anymore collected here, use .csv output (and MATLAB instead)
//
// mikael.mieskolainen@cern.ch, 2018


// Regge factorization plots
void CombinatoricsSuper::PlotAllFactorization() {

	printf("CombinatoricsSuper::PlotAllFactorization:: \n");

	TCanvas* c1 = new TCanvas(Form("c1_Regge"),"Regge factorization",1000,1600);
	c1->Divide(1,2);

	c1->cd(1)->SetGrid();
	c1->cd(1)->SetLogy();
   	
	const Int_t d  = sources_.at(0)->d_;
   	const UInt_t n = std::pow(2,d);

	std::vector<TGraphErrors*> gr(sources_.size(), NULL);

	// GET LR sequence
	std::vector<UInt_t> LR_ind = LRsequence(d);

	// Calculate the factorization sequence
	std::vector<UInt_t> FR_ind(n, 0);

	for (UInt_t i = 0; i < n; ++i) {

		std::vector<Bool_t> orig_vec = Ind2Vec(i, d);
		std::vector<Bool_t> refl_vec = Ind2Vec(LR_ind.at(i), d);

		// Calculate the target vector by bitwise OR
		std::vector<Bool_t> targ_vec(d,0);
		for (Int_t k = 0; k < d; ++k) {
			targ_vec.at(k) = (orig_vec.at(k) | refl_vec.at(k)); 
		}

		// Turn to index representation
		FR_ind.at(i) = Vec2Ind(targ_vec);
	}

	// Entropies
	//std::vector<Double_t> H(5,0);

	// Kullback-Leibler divergence values
	std::vector<Double_t> KLs(sources_.size(), 0);

	// reordered indices
	std::vector<UInt_t> sortind;

	// Loop over Data / MC sources
	for (UInt_t i = 0; i < sources_.size(); ++i) {

		std::vector<Double_t> x_count = sources_.at(i)->GetX();
		x_count.at(0) = 0; // Nullify the zeroth component
		Double_t N_sum = vsum(x_count);

		TString name = sources_.at(i)->GetName();

		std::vector<Double_t> x(n, 0);
		std::vector<Double_t> x_e(n, 0);

		std::vector<Double_t> y(n, 1); // Set to unity by default
		std::vector<Double_t> y_e(n, 0);

		// Get Binomial/Multinomial counting errors
		std::vector<Double_t> x_count_err = GetBinomError(x_count);

		for (UInt_t j = 1; j < x_count.size(); ++j) { // Loop by starting from the first

			x[j] = j;
			y[j] = TMath::Sqrt(x_count.at(j)/N_sum * x_count.at(LR_ind.at(j))/N_sum) / ((x_count.at(FR_ind.at(j)) / N_sum) + 1e-12);

			// Taylor expanded x/y error
			y_e[j] = 0;//GetRatioError(x_count.at(j), x_count_err.at(j), x_count.at(LR_ind.at(j)), x_count_err.at(LR_ind.at(j)));

			// Trivial identities (idempotent) we set to unity
			/*if (j == LR_ind.at(j)) {
				y[j] = 1;
				y_e[j] = 1e-6;
			}*/
		}

		// Data, i.e., i == 0 defines the sort order
		if (i == 0) {
			sortind = vsortind(y);
		}

		// Sort
		vsort(y, sortind);
		vsort(y_e, sortind);


		gr[i] = new TGraphErrors(x.size(), &(x[0]), &(y[0]), &(x_e[0]), &(y_e[0]));

   		gr[i]->SetMarkerColor(colors_.at(i));
		gr[i]->SetMarkerStyle(markers_.at(i));

   		if (i == 0) {
   			gr[i]->Draw("AP");
   		} else {
   			gr[i]->Draw("P");
   		}

		//gr[i]->GetYaxis()->SetRangeUser(-0.03, 0.03);
   		gr[i]->GetXaxis()->SetLimits(-1, n);
   		gr[i]->SetTitle(";Vector ID triplet (reordered); Factorization ratio #sqrt{(P(i)P(j))}/P(k)");
	}


	TLegend* leg = new TLegend(0.1408046,0.6970339,0.4204023,0.8686441);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);

	for (UInt_t i = 0; i < sources_.size(); ++i) {
		leg->AddEntry(gr[i],Form("%s", sources_.at(i)->GetName().Data() ),"LP");
	}

	//leg->Draw();
	//c1->SaveAs("./figures_xsec/%d/CrossSections/Regge_Factorization.pdf", sources_.at(0).GetRunNumber() );



	// -------------------------------------------------------------------
	// RATIO Plot
	c1->cd(2);

	std::vector<TGraphErrors*> grRatios(sources_.size(), NULL);

	std::vector<Double_t> x_ref = sources_.at(0)->GetX();
	x_ref.at(0) = 0; // Nullify the zeroth component
	std::vector<Double_t> x_ref_n = nvec(x_ref);

	std::vector<Double_t> x_ref_err   = GetBinomError(x_ref);
	std::vector<Double_t> x_ref_err_n = scalevec(x_ref_err, 1.0/vsum(x_ref));

	// Loop over different sources
	for (UInt_t i = 0; i < sources_.size(); ++i) {

		std::vector<Double_t> x_count   = sources_.at(i)->GetX();
		x_count.at(0) = 0; // Nullify the zeroth component
		std::vector<Double_t> x_count_n = nvec(x_count);

		std::vector<Float_t> x(n, 0);
		std::vector<Float_t> x_e(n, 0);

		std::vector<Float_t> y(n, 1); // Set to unity by default
		std::vector<Float_t> y_e(n, 0);

		// Get Binomial/Multinomial counting errors
		std::vector<Double_t> x_count_err   = GetBinomError(x_count);
		std::vector<Double_t> x_count_err_n = scalevec(x_count_err, 1.0/vsum(x_count));

		for (UInt_t j = 1; j < x_count.size(); ++j) {
			x[j] = j;
			y[j] = TMath::Sqrt(x_count.at(j) * x_count.at(LR_ind.at(j))) / (x_count.at(FR_ind.at(j)) + 1e-12)   /
				   (TMath::Sqrt(x_ref.at(j) * x_ref.at(LR_ind.at(j)) ) / (x_ref.at(FR_ind.at(j)) + 1e-12) + 1e-12);

			// Trivial identities (idempotent) we set to unity
			/*if (j == LR_ind.at(j)) {
				y[j] = 1;
				y_e[j] = 1e-6;
			}*/

			// Taylor expanded x/y error
			y_e[j] = 0;
		}

		grRatios[i] = new TGraphErrors(x.size(), &(x[0]), &(y[0]), &(x_e[0]), &(y_e[0]));

   		grRatios[i]->SetMarkerColor(colors_.at(i));
		grRatios[i]->SetMarkerStyle(markers_.at(i));

   		if (i == 0) {
   			grRatios[i]->Draw("ALP");
   		} else {
   			grRatios[i]->Draw("P");
   		}
			grRatios[i]->SetTitle(";Vector ID triplet (reordered); MC/Data");

   		grRatios[i]->GetYaxis()->SetRangeUser(0,2);
   		grRatios[i]->GetXaxis()->SetLimits(-1, n);
	}


	// ADD a second legend to upper plot (which has the KS values)
	c1->cd(1);

	TLegend* leg2 = new TLegend(0.735, 0.131, 0.872, 0.316);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.03);
	c1->SaveAs(Form("./figures_xsec/%d/CrossSections/Regge_Factorization.pdf", sources_.at(0)->GetRunNumber() ));


	for (UInt_t i = 0; i < sources_.size(); ++i) { 
		delete gr[i]; 
		delete grRatios[i];
	}

	delete leg2;
	delete leg;
	delete c1;
}


// Do combinatorial rate combined plots
void CombinatoricsSuper::PlotAll(TString det_or_gen, Bool_t LR_mode_on, Bool_t sort_on) {

	printf("CombinatoricsSuper::PlotAll:: \n");

	// Probability canvas
	TCanvas* cx = new TCanvas(Form("cx"), "Distribution of probabilities", 400, 300);
	cx->cd()->SetLogy();
	cx->cd()->SetLogx();

	std::vector<TH1D> hP;
	for (UInt_t i = 0; i < sources_.size(); ++i) {
		TH1D hist = TH1D(Form("test_%d",i), ";P;f(P);", 100, 0, 0.01);
		hP.push_back(hist);
	}

	// Main canvas
	TCanvas* c1 = new TCanvas(Form("c1_%d", LR_mode_on),"Partial Rates",1000,1600);
	c1->Divide(1,2);

	c1->cd(1)->SetGrid();
   	c1->cd(1)->SetLogy();
	
	const Int_t d_ = sources_.at(0)->d_;
   	const Int_t n  = pow(2,d_);

	std::vector<TGraphErrors*> gr(sources_.size(), NULL);

	// GET LR sequence
	std::vector<UInt_t> LR_ind = LRsequence(d_);

	// Entropies
	std::vector<Double_t> H(sources_.size(),0);

	// Kullback-Leibler divergence values
	std::vector<Double_t> KLs(sources_.size(), 0);


	// reordered indices
	std::vector<UInt_t> sortind;

	// Loop over Data / MC sources
	for (UInt_t i = 0; i < sources_.size(); ++i) {

		// Detector or generator level
		std::vector<Double_t> x_count;
		if (det_or_gen.Contains("Detector")) {
			x_count = sources_.at(i)->GetX();
		}
		if (det_or_gen.Contains("Generator")) {
			x_count = sources_.at(i)->GetXGen();
		}
		x_count.at(0) = 0; // Nullify the zeroth component
		Double_t N_sum = vsum(x_count);


		// Calculate Entropy
		H.at(i) = calcH(nvec(x_count));

		TString name = sources_.at(i)->GetName();

		std::vector<Double_t> x(n, 0);
		std::vector<Double_t> x_e(n, 0);

		std::vector<Double_t> y(n, 0);
		std::vector<Double_t> y_e(n, 0);


		// Get Binomial/Multinomial counting errors
		std::vector<Double_t> x_count_err = GetBinomError(x_count);

		for (UInt_t j = 0; j < x_count.size(); ++j) {
			x[j] = j;

			if (LR_mode_on) {

				y[j]   = x_count.at(j) / (x_count.at(LR_ind.at(j)) + 1e-12);

				// Taylor expanded x/y error (with multinomial Cov(x,y) = -np_ip_j)
				y_e[j] = GetRatioError(x_count.at(j), x_count_err.at(j), x_count.at(LR_ind.at(j)), x_count_err.at(LR_ind.at(j)), -N_sum*x_count.at(j)/N_sum*x_count.at(LR_ind.at(j))/N_sum );

				// Trivial identities (idempotent) we set to unity
				if (j == LR_ind.at(j)) {
					y[j] = 1;
					y_e[j] = 1e-12;
				}
				
			} else {
				y[j]   = x_count.at(j) / N_sum;
				y_e[j] = x_count_err.at(j) / N_sum;

				hP.at(i).Fill( x_count.at(j) / N_sum );
			}
		}

		// IF SORT
		if (sort_on) {

			// Data, i.e., i == 0 defines the sort order
			if (i == 0) {
				sortind = vsortind(y);
			}

			// Sort
			vsort(y, sortind);
			vsort(y_e, sortind);
		}

		// Create graph
		gr[i] = new TGraphErrors(x.size(), &(x[0]), &(y[0]), &(x_e[0]), &(y_e[0]));
   		gr[i]->SetMarkerColor(colors_.at(i));
		gr[i]->SetMarkerStyle(markers_.at(i));

   		if (i == 0) {
   			gr[i]->Draw("AP");

   			// Fit 
			//TF1* f = new TF1("f", "[2] * x * x + [1] * x + [0]");
  			//gr[i]->Fit(f);
  			//gr[i]->Draw("AL");

   		} else {
   			gr[i]->Draw("P");
   		}

   		if (LR_mode_on) {
			gr[i]->GetYaxis()->SetRangeUser(0.1, 10);
   			gr[i]->GetXaxis()->SetLimits(-1, n);
   			gr[i]->SetTitle(";Vector ID pair (reordered); C#leftrightarrowA Ratio");
   		} else {
			gr[i]->GetYaxis()->SetRangeUser(1e-6, 1);
	   		gr[i]->GetXaxis()->SetLimits(-1, n);
   			gr[i]->SetTitle(";Vector ID (reordered); Probability");
   		}
	}

	// Create legend
	TLegend* leg = new TLegend(0.168, 0.690, 0.338, 0.865);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);

	if (!LR_mode_on) {
		for (UInt_t i = 0; i < sources_.size(); ++i) {
			if (LR_mode_on) {
				leg->AddEntry(gr[i],Form("%s", sources_.at(i)->GetName().Data() ),"LP");
			} else {
				leg->AddEntry(gr[i],Form("%s, H = %0.2f", sources_.at(i)->GetName().Data(), H.at(i) ),"LP");
			}
		}
		leg->Draw();
	}

	// -------------------------------------------------------------------
	// RATIO Plot
	c1->cd(2);

	std::vector<TGraphErrors*> grRatios(sources_.size(), NULL);

	// Detector or generator level
	std::vector<Double_t> x_ref;
	if (det_or_gen.Contains("Detector")) {
		x_ref = sources_.at(0)->GetX();
	}
	if (det_or_gen.Contains("Generator")) {
		x_ref = sources_.at(0)->GetXGen();
	}
	x_ref.at(0) = 0; // Nullify the zeroth component
	std::vector<Double_t> x_ref_n = nvec(x_ref);


	std::vector<Double_t> x_ref_err   = GetBinomError(x_ref);
	std::vector<Double_t> x_ref_err_n = scalevec(x_ref_err, 1.0/vsum(x_ref));

	// Loop over all sources
	for (UInt_t i = 0; i < sources_.size(); ++i) {

		// Detector or generator level
		std::vector<Double_t> x_count;
		if (det_or_gen.Contains("Detector")) {
			x_count = sources_.at(i)->GetX();
		}
		if (det_or_gen.Contains("Generator")) {
			x_count = sources_.at(i)->GetXGen();
		}
		x_count.at(0) = 0; // Nullify the zeroth component
		std::vector<Double_t> x_count_n = nvec(x_count);

		std::vector<Double_t> x(n, 0);
		std::vector<Double_t> x_e(n, 0);

		std::vector<Double_t> y(n, 1); // Set to unity by default
		std::vector<Double_t> y_e(n, 0);

		// Get Binomial/Multinomial counting errors
		std::vector<Double_t> x_count_err   = GetBinomError(x_count);
		std::vector<Double_t> x_count_err_n = scalevec(x_count_err, 1.0/vsum(x_count));

		for (UInt_t j = 1; j < x_count.size(); ++j) {
			x[j] = j;

			if (LR_mode_on) {

				y[j]   = (x_count_n.at(j) / (x_count_n.at(LR_ind.at(j)) + 1e-12)) / (x_ref_n.at(j)/(x_ref_n.at(LR_ind.at(j)) + 1e-12) + 1e-12);

				// (Double) Taylor expanded x/y error
				y_e[j] = GetRatioError(x_count_n.at(j) / (x_count_n.at(LR_ind.at(j))+1e-12),
										GetRatioError(x_count_n.at(j), x_count_err_n.at(j), x_count_n.at(LR_ind.at(j)), x_count_err_n.at(LR_ind.at(j)), 0 ),
					                    x_ref_n.at(j) / (x_ref_n.at(LR_ind.at(j))+1e-12),
					                    GetRatioError(x_ref_n.at(j), x_ref_err_n.at(j), x_ref_n.at(LR_ind.at(j)), x_ref_err_n.at(LR_ind.at(j)), 0),
					                    0);

				// Trivial identities (idempotent) we set to unity
				if (j == LR_ind.at(j)) {
					y[j] = 1;
					y_e[j] = 1e-12;
				}
			} else {
				y[j]   = x_count_n.at(j) / (x_ref_n.at(j) + 1e-12);

				// Taylor expanded x/y error
				y_e[j] = GetRatioError(x_count_n.at(j), x_count_err_n.at(j), x_ref_n.at(j), x_ref_err_n.at(j), 0);
			}
		}

		// IF SORT
		if (sort_on) {
			// Sort
			vsort(y, sortind);
			vsort(y_e, sortind);
		}

		// Draw
		grRatios[i] = new TGraphErrors(x.size(), &(x[0]), &(y[0]), &(x_e[0]), &(y_e[0]));
   		grRatios[i]->SetMarkerColor(colors_.at(i));
		grRatios[i]->SetMarkerStyle(markers_.at(i));

   		if (i == 0) {
   			grRatios[i]->Draw("ALP");
   		} else {
   			grRatios[i]->Draw("P");
   		}
   		if (LR_mode_on) {
   			grRatios[i]->SetTitle(";Vector ID pair (reordered); MC/Data");
   		} else {
   			grRatios[i]->SetTitle(";Vector ID (reordered); MC/Data");
   		}
   		grRatios[i]->GetYaxis()->SetRangeUser(0,2);
   		grRatios[i]->GetXaxis()->SetLimits(-1, n);

   		// Calculate KL divergence (DATA | MC)
   		KLs.at(i) = calcKL(x_ref_n, x_count_n);
	}

	// ADD a second legend the to upper plot (which has the KS values)
	c1->cd(1);

	TLegend* leg2 = new TLegend(0.733, 0.692, 0.870, 0.865);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.03);

	for (UInt_t i = 1; i < sources_.size(); ++i) {
		leg2->AddEntry(grRatios[i],Form("KL = %0.3f", KLs.at(i)), "LP");
	}
	
	if (!LR_mode_on)
		leg2->Draw();

	if (LR_mode_on) {
		c1->SaveAs(Form("./figures_xsec/%d/CrossSections/Combinatorial_CA_Ratios_%s.pdf", sources_.at(0)->GetRunNumber(), det_or_gen.Data() ));
	} else {
		c1->SaveAs(Form("./figures_xsec/%d/CrossSections/Combinatorial_Rates_%s.pdf", sources_.at(0)->GetRunNumber(), det_or_gen.Data() ));		
	}

	// Delete memory allocations
	for (UInt_t i = 0; i < sources_.size(); ++i) { 
		delete gr[i]; 
		delete grRatios[i]; 
	}
	
	delete leg;
	delete leg2;
	delete c1;


	// Create test plot
	cx->cd();
	for (UInt_t i = 0; i < hP.size(); ++i) {

		hP.at(i).SetLineColor(colors_.at(i));

		Double_t norm = hP.at(i).GetEntries();
		if (norm) hP.at(i).Scale(1/norm);

		if (i == 0) {
			hP.at(i).Draw("P");
		} else {
			hP.at(i).Draw("SAME");
		}
	}

	//cx->SaveAs(Form("./figures_xsec/%d/CrossSections/ProbProb_Distribution_%s.pdf", sources_.at(0)->GetRunNumber(), det_or_gen.Data() ));
	delete cx;
}


void Plot() {

	// ------------------------------------------------------------------------------
	// Plot total process cross sections as a function of total inelastic

	TCanvas* c1 = new TCanvas(Form("%s_%s", sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()),"Lagrangian", 600, 900);
	c1->Divide(1,2);

	TGraph* g_SDL = new TGraph(inel_vals.size(), &(inel_vals[0]), &(tot_SDL[0]));
	TGraph* g_SDR = new TGraph(inel_vals.size(), &(inel_vals[0]), &(tot_SDR[0]));
	TGraph* g_DD = new TGraph(inel_vals.size(), &(inel_vals[0]), &(tot_DD[0]));
	TGraph* g_ND = new TGraph(inel_vals.size(), &(inel_vals[0]), &(tot_ND[0]));
	g_SDL->SetTitle("Unvisible #sigma_{SDL}, #sigma_{SDR}, #sigma_{DD} with dashed lines;#sigma_{inel}^{tot} [mb]; #sigma_{process}^{tot} [mb]");

	g_SDL->SetLineColor(8);
	g_SDR->SetLineColor(9);
	g_DD->SetLineColor(12);
	//g_SDR->SetMarkerStyle(markers[i]);

	Double_t y_range[2] = {0.0, 15.0};

	g_SDL->GetYaxis()->SetRangeUser(y_range[0], y_range[1]);

	c1->cd(1);
	g_SDL->Draw("AC");
	g_SDR->Draw("C");
	g_DD->Draw("C");
	g_ND->Draw("C");

	// -------------------------------------------------------------------
	// Find out the interval where SDL ~== SDR
	Double_t tolerance = 0.1; // mb
	UInt_t start_ind = 0; Bool_t start_marked = kFALSE;
	UInt_t end_ind = 0;

	for (UInt_t i = 0; i < tot_SDL.size(); ++i) {

		if ((TMath::Abs(tot_SDL.at(i) - tot_SDR.at(i)) < tolerance) && start_marked == kFALSE) {
			start_ind = i;
			start_marked = kTRUE;
		}
		if (TMath::Abs(tot_SDL.at(i) - tot_SDR.at(i)) < tolerance) {
			end_ind = i;
		}
	}

	// Draw total inelastic min and max with red dashed
	TLine* l_min = new TLine(inel_vals.at(start_ind),y_range[0],inel_vals.at(start_ind),y_range[1]); // (x1,y1,x2,y2)
	l_min->SetLineColor(2);
	l_min->SetLineStyle(3); // Points
	l_min->Draw("C");

	TLine* l_max = new TLine(inel_vals.at(end_ind),y_range[0],inel_vals.at(end_ind),y_range[1]); // (x1,y1,x2,y2)
	l_max->SetLineColor(2);
	l_max->SetLineStyle(3); // Points
	l_max->Draw("C");	

	// Draw visible inelastic with black dashed
	TLine* l3 = new TLine(inel_vals.at(0),y_range[0],inel_vals.at(0),y_range[1]); // (x1,y1,x2,y2)
	l3->SetLineStyle(3); 	// Points
	l3->Draw("C");


	// --------------------------------------------------------------------
	// Draw unvisible cross sections between feasible domain (between red lines)

	// Create efficiency vectors
	std::vector<Double_t> feasible_inel;
	std::vector<Double_t> eff_SDL;
	std::vector<Double_t> eff_SDR;
	std::vector<Double_t> eff_DD;

	for (UInt_t i = 0; i < inel_vals.size(); ++i) {
		if (i >= start_ind && i <= end_ind) {
			feasible_inel.push_back( inel_vals.at(i));
			eff_SDL.push_back(  tot_SDL.at(i) - Pvec[0][0]);
			eff_SDR.push_back(  tot_SDR.at(i) - Pvec[1][0]);
			eff_DD.push_back(   tot_DD.at(i)  - Pvec[2][0]);
		}
	}

	TGraph* g_eff_SDL = new TGraph(feasible_inel.size(), &(feasible_inel[0]), &(eff_SDL[0]));
	TGraph* g_eff_SDR = new TGraph(feasible_inel.size(), &(feasible_inel[0]), &(eff_SDR[0]));
	TGraph* g_eff_DD  = new TGraph(feasible_inel.size(), &(feasible_inel[0]), &(eff_DD[0]));


	g_eff_SDL->SetLineColor(8);
	g_eff_SDL->SetLineStyle(2);
	g_eff_SDL->Draw("C");

	g_eff_SDR->SetLineColor(9);
	g_eff_SDR->SetLineStyle(2);
	g_eff_SDR->Draw("C");

	g_eff_DD->SetLineColor(12);
	g_eff_DD->SetLineStyle(2);
	g_eff_DD->Draw("C");

	// -------------------------------------------------------------------
	// Draw Legend
	TLegend* leg = new TLegend(0.6592248,0.2228593,0.8389947,0.3948777);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);

	leg->AddEntry(l3, "Visible inelastic", "L");
	leg->AddEntry(g_SDL,Form("SDL"),"LP");
	leg->AddEntry(g_SDR,Form("SDR"),"LP");
	leg->AddEntry(g_DD, Form("DD"), "LP");

	leg->Draw();

	// REDRAW because legend might overlap
	c1->cd(1);
	g_SDL->Draw("C");
	g_SDR->Draw("C");
	g_DD->Draw("C");
	g_ND->Draw("C");

	// -------------------------------------------------------------------

	// Plot the Lagrangian cost as a function of total inelastic // Just put the last KL value
	TGraph* g = new TGraph(inel_vals.size(), &(inel_vals[0]), &(L_vals[0]));
	g->SetTitle(Form("Input: %s, Model: %s, KL(P|Q) = %0.3f;#sigma_{inel}^{tot} [mb]; Lagrangian cost", 
		sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data(), KL.at(KL.size()-1)) );

	c1->cd(2)->SetLogy();
	g->Draw("AC");

	c1->SaveAs(Form("./figures_xsec/%d/CrossSections/Total_Inelastic_Input_%s_Model_%s.pdf", sources_.at(0)->GetRunNumber(), sources_.at(data_index)->GetName().Data(),
										 									   		   sources_.at(MC_index)->GetName().Data() ));

	// SD/DD etc. ratios
	// -------------------------------------------------------------------
	TCanvas* c2 = new TCanvas(Form("%s_%s_ratios", sources_.at(data_index)->GetName().Data(), sources_.at(MC_index)->GetName().Data()), "sigmaratios", 600, 900);
	c2->Divide(1,2);

	//TGraph* g_SD_over_DD = new TGraph(feasible_inel.size(), &(feasible_inel[0]), &(eff_SDL[0]));
	//TGraph* g_SDplusDD_over_inel = new TGraph(feasible_inel.size(), &(feasible_inel[0]), &(eff_SDR[0]));

	std::vector<Double_t> SD_over_inel(inel_vals.size(), 0);
	std::vector<Double_t> DD_over_inel(inel_vals.size(), 0);

	for (UInt_t i = 0; i < inel_vals.size(); ++i) {
		SD_over_inel.at(i) = (tot_SDL.at(i) + tot_SDR.at(i)) / inel_vals.at(i);
		DD_over_inel.at(i) =  tot_DD.at(i) / inel_vals.at(i);
	}

	TGraph* g_SD_over_inel = new TGraph(inel_vals.size(), &(inel_vals[0]), &(SD_over_inel[0]));
	TGraph* g_DD_over_inel = new TGraph(inel_vals.size(), &(inel_vals[0]), &(DD_over_inel[0]));

	g_SD_over_inel->SetTitle(";#sigma_{inel}^{tot} [mb]; Ratio");

	g_SD_over_inel->SetLineColor(9);
	g_DD_over_inel->SetLineColor(12);


	// Draw ratios
	c2->cd(1);

	Double_t y_range2[2] = {0.0, 0.3};
	g_SD_over_inel->GetYaxis()->SetRangeUser(y_range2[0], y_range2[1]);

	g_SD_over_inel->Draw("AC");
	g_DD_over_inel->Draw("C");

	TLegend* leg2 = new TLegend(0.6592248,0.2228593,0.8389947,0.3948777);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.03);

	//leg2->AddEntry(l3, "Visible inelastic", "L");
	leg2->AddEntry(g_SD_over_inel,"#sigma_{SD}/#sigma_{inel}","LP");
	leg2->AddEntry(g_DD_over_inel,"#sigma_{DD}/#sigma_{inel}","LP");
	leg2->Draw();

	// Draw vertical lines
	l_min->Draw("AC");
	l_max->Draw("C");
	l3->Draw("C");

	// Draw Lagrangian cost
	c2->cd(2)->SetLogy();
	g->Draw("AC");
/*
	c2->SaveAs(Form("./figures_xsec/%d/CrossSections/SD_DD_Ratios_Input_%s_Model_%s.pdf", 
		sources_.at(0)->GetRunNumber(), sources_.at(data_index)->GetName().Data(),
										sources_.at(MC_index)->GetName().Data() ));
*/
	delete g_SDL;
	delete g_SDR;	
	delete g_DD;
	delete g_ND;

	delete l_min;
	delete l_max;
	delete l3;

	delete leg;
	delete c1;
	
}