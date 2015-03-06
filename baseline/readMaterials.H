
	// Read Materials
	
	// Read TEG Base
	base.t = 0.008;  // 8 mm thick wall as given in dwg
	base.k = kCu;  // Assuming to be made of Stainless Steel or Cu

	// Read Thermal Grease
	// For Omegatherm used for BiTe
	// k - 2.30765, t = 0.1 mm
	// grafoil laminate	k - 5. t - 1 mm
	tgrease.gap = 0.001;
	tgrease.k = 5;
	
	// Read Insulation Properties
	insulation.t = 0.002;
	insulation.k = 0.0334;
	insulation.epsilon = epsInsulated;

	//Read TEM
	#include "readTEM.H"
	
	//Read Air Properties
	string fname = cwd() + "\\..\\..\\..\\Props\\Air\\Cp_Air.txt";
	air.CpTable = readProp(fname.c_str(), air.size);
	
	fname = cwd() + "\\..\\..\\..\\Props\\Air\\Kt_Air.txt";
	air.kTable = readProp(fname.c_str(), air.size);
	
	fname = cwd() + "\\..\\..\\..\\Props\\Air\\density_Air.txt";
	air.densityTable = readProp(fname.c_str(), air.size);
	
	fname = cwd() + "\\..\\..\\..\\Props\\Air\\mu_Air.txt";
	air.muTable = readProp(fname.c_str(), air.size);
	
	TEMat[0].name = "ntype-SK";
	TEMat[1].name = "ptype-SK";
	TEMat[2].name = "ntype-BiTE";
	TEMat[3].name = "ptype-BiTE";
		
	// Read TE Materials
	for(int i = 0;i < 4;i++)
	{
		fname = cwd() + "\\..\\..\\..\\Props\\" + TEMat[i].name + "\\lambda.txt";
		TEMat[i].lambda = readProp(fname.c_str(), TEMat[i].lambdasize);
		
		fname = cwd() + "\\..\\..\\..\\Props\\" + TEMat[i].name + "\\alpha.txt";
		TEMat[i].alpha = readProp(fname.c_str(), TEMat[i].alphasize);
		
		fname = cwd() + "\\..\\..\\..\\Props\\" + TEMat[i].name + "\\rho.txt";
		TEMat[i].rho = readProp(fname.c_str(), TEMat[i].rhosize);
		
		fname = cwd() + "\\..\\..\\..\\Props\\" + TEMat[i].name + "\\ZT.txt";
		TEMat[i].ZT = readProp(fname.c_str(), TEMat[i].ZTsize);
	}
