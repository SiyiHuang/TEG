
	// Defining TE Modules
	// Type 0 module is Skutterudite
	TEM.length[0] = 0.0508;	//m
	TEM.width[0] = 0.0508;	//m
	TEM.height[0] = 0.004;	//m TEC leg height but not the actual TEM height. Modules height = 4 mm + 1.5x2 = 7 mm
	TEM.nC[0] = 32;
	TEM.etaNP[0] = eta_TE;	
	TEM.numDensity[0]  = TEM.nC[0]/(TEM.length[0]*TEM.width[0]); 	// Number of TE couples per unit Area
	
	// Type 1 module is BiTe
	TEM.length[1] = 0.04013;	//m
	TEM.width[1] = 0.04013;	//m
	TEM.height[1] = 0.00201;	//m TEC leg height but not the actual TEM height. Modules height = 4 mm
	TEM.nC[1] = 127;
	TEM.etaNP[1] = eta_TE;	
	TEM.numDensity[1]  = TEM.nC[1]/(TEM.length[1]*TEM.width[1]); 	// Number of TE couples per unit Area

	TEM.epsilon = epsModule;

	//TEC.Area = 2.0 * 0.004 * 0.004; // TEM.Area[iCV]*TEM.etaNP[TEM.Type[iCV]];  [4] Area of one leg x 2 for SKu
	//TEC.Area = 2.0 * 0.002 * 0.002; // TEM.Area[iCV]*TEM.etaNP[TEM.Type[iCV]];  [4] Area of one leg x 2 for BiTE