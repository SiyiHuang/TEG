	

//  Iteration Variables	
	int i, j, itr, itrmax, iCV;
    
//  Grid Parameters	
	int Ncv;
		
    // Geometry
	double pipeDia, Ain, Aout, Abox, ARin, ARout, muarea;
	
	double x[Nx+1];
	double xinit, xend, height, width, flowLength, diaHyd, AR;
	double heightB, widthB, flowLengthB;
	double Areainit, Areaend;
	double AreaCV[Nx], AreaModule[Nx], AreaTE[Nx], AreaInsulated[Nx];
	double volume, pmin, pmax, pmin1, pmin2, pmax1, pmax2, p , dp;
	double **geometry;
	int geosize;
	
	double eta_module, eta_TE;
	
	// Temperature Arrays
	double Tavg, Tginit, Tinf;
	double T[15][Nx],R[15][Nx],Q[15][Nx];
	double Tg[Nx+1], Tc[Nx], Th[Nx];
	
	
	// Convergence 	
	double error, tol;
	    	 
	// Dummy Variables    
    double Tgavg, Tgnew, Tcinit;
	
	// Arrays for Thermo-Electrics	
	double height_TEM, length_TEM, Area_TEM;
	double rho_[Nx], lambda_[Nx], Ri_[Nx], Rl_[Nx], Kl_[Nx], S_[Nx];
	
	// Arrays for Output
	double Pout_[Nx];
	double Qh_[Nx], Qg_[Nx], Qc_[Nx], Qrad_[Nx], Qconv_[Nx], Qhtem_[Nx]; 
	double QhTot_[Nx], QgTot_[Nx], QcTot_[Nx], QradTot_[Nx], PoutTot_[Nx], QconvTot_[Nx], QhtemTot_[Nx],Vopenc_[Nx],nTETot_[Nx][2], ZT_[Nx];
	double ReDh_, delP_[Nx];
   	double PoutTotcum, QhTotcum, QgTotcum, QcTotcum, PressDropTot, QradTotcum, QconvTotcum, QhtemTotcum, nTETotcum0, nTETotcum1;
	double PDin, PDout, AllowablePD;

	double mdot;    
		
	double gasSide[Nx], TESide[Nx];
   	
	// Geometry	
	double dx;
	
	// Number of TE modules N_TE(x)
	int n_TE[Nx], nTEmax, nTEcum, nTEtot;
	
	// Arrays for Gas
	double cpg_[Nx], kg_[Nx];
	double rhoAir_[Nx], muAir_[Nx], kAir;
	double Kc, Ke, ReDh, f, delP, vel;
	
	// Variables for fin
	double kCu, kfin, kAl, kSS;
	double finSpacing, tfin, pitch;
	int Nf;
	double At, Ab, Af, Ac, mch, velch, Rtot_[Nx];
	double finSpacingB, tfinB;

	// Strings for directories
	string outputDir, currDir;

	//Radiation
	double epsInsulated, epsModule;

	// Extra Arrays for optimization
	const int size = 2500;
	int index = 0, sCV;
	double h[size], w[size], fl[size], Powermax[size], Tgdrop[size], Qgas[size], Qrad[size], Qh[size], Qc[size],Pres[size],tf[size], fs[size], nfin[size];
	double mdt[size],apd[size], pin[size], pout[size], nTEC0[size], nTEC1[size];

	double nTEM0, nTEM1;
