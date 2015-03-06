/***********************************C++*************************************************/
	// User Inputs
	// Inlet and Outlet Diameters
	pipeDia = 0.0635; // 2.5 inch, 63.5 mm
	Ain = 3.1416*pipeDia*pipeDia/4.0;
	Aout = Ain;

	
	// Material Properties
	kCu = 401;	
	kAl = 237;
	kSS = 20;   // 16.7 0-100C, 21.5 at 500C
	kfin = kCu;
	
    // Flow Input properties
	Tginit = 550 + 273.15;
	mdot =  0.035;
	Tinf = 100.0 + 273.15;
	Tcinit = 100.0 + 273.15;
	
    // TEG Geometry Inputs
    Ncv = 10;
	volume = 0.003592148;
	//Amodule_Max = 0.002581 * 25.0;  // 50 total - 25 on either side
	
	// TEM Inputs
	eta_module = 0.8; // %area covered by Modules
	eta_TE = 0.5; 
	epsInsulated = 0.75;
	epsModule = 0.55;

	// TEG Default values
	//heightB = 0.038;  
	//widthB = 0.224;
	//flowLengthB = 0.4131;
	heightB = 0.038;
	widthB = 0.224;
	flowLengthB = 0.4131;

	tfinB = 0.003302; // 0.13" for Cu fins. 0.065" for SS
	finSpacingB = 0.00635; // 0.25" Apart fins

	// Simulation controls
	itrmax = 100;
	tol = 1e-04;

/****************************************************************************************/