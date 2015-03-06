
namespace TEG
{
	#define Nx 500
	#define SMALL 1E-16
	#define LARGE 1e16;
	 
// extern tells the compiler this variable is declared elsewhere
// Variables for TE Modules
	struct TEModules
	{
		double numDensity[2];	// Number of TE couples per unit Area
		double etaNP[2];	// Area proportion for n-p Legs
		int nC[2];	// Number of TE couples per Module
		double length[2];
		double width[2];
		double height[2];	
		double Area[Nx];
		double nTE[Nx];
		int Type[Nx];
		double epsilon;	// Radiation
	} TEM;
	
	struct TECouple
	{
		double Type[Nx];
		double height;
		double Area;
		double Alpha;
		double R;
		double Kl;
		double ZT;
	}	TEC;
	
	struct TEMaterials
	{
		string name;
		double **lambda;
		double **rho;
		double **alpha;
		double **ZT;
		int lambdasize;
		int rhosize;
		int alphasize;
		int ZTsize;
	}	TEMat[4];
		
	struct Air
	{
		double **CpTable;
		double	**kTable;
		double **densityTable;
		double **muTable;
		int size;
		
		double Cp;
		double k;
		double density;
		double mu;
		
	}	air;

	struct Insulation
	{
		double t;	// Thickness
		double k;	// Thermal Conductivity
		double Area[Nx];	// Area		
		double epsilon;
	}	insulation;

	struct TEGbase
	{
		double t;	// Thickness
		double k;	// Thermal Conductivity
		double Area[Nx];	// Area				
	} base;

	struct Thermal_grease
	{
		double gap;	// Gap
		double k;	// Thermal Conductivity
		double Area[Nx];	// Area				
	} tgrease;

	double Ri[Nx], Rl[Nx], V[Nx], I[Nx], P[Nx], ZT[Nx];
	double nTE[Nx];
	

// Thermal Resistances

// R01 - The fin resistances includes contribution from the fin-base Area and extended fin surfaces
double Req( double Tavg, double AR, double diaHyd, double height, double tfin, double kfin, int Nf, double At, double Af);

// R2- The contact resistance between fin and TEG
double Rcfin( double dummy);

// R3 - The resistance contribution by the TEG base material
double Rbasethickness(int iCV);

// R4 - The contact resistance between TEM hot junction and hot surface.
double RcTEMH( int iCV);

// R5 - The resistance by insulation material (min K)
double Rinsulation(int iCV);

// R6 - The effective internal resistance of TEC as seen from hot side
double RintTEM(double Th, double Tc, int iCV);

// R7 - Computing Radiation resistance for the TEM Hot surface
double RradTEM(double Th, double Tinf, int iCV);

// R8 - Computing Resistances for convection through Insulation
double RconvIns(double beta, double alpha, double nu, double deltaT, double L, double Area);

// R9 - Computing Radiation resistance for the Insulation material 
double RradIns(double Tins, double Tinf, int iCV);

// R10 - The contact resistance between TEM cold junction and coolant
double RcTEMC( int iCV);

// R11 - The equivalent Load resistance
double ReqLoad(double Th, double Tc, double Tcoolant, int iCV);

// Compute TE circuit
void computeTECircuit(double Th, double Tc, int iCV);

// Correlations
double fReCorr(double x);

double NuCorr(double AR);

string cwd();

double min(double a, double b);

void updateAir(double Tgavg);

void updateTECouples(double Tavg, int iCV);

double electricCurrent(double Th, double Tc, double Rl, int iCV);

void Thot
(
	double Tgavg, double Tc,
	double AR, double diaHyd, double height,double tfin, double kfin, int Nf, double At, double Af, double Ab,
	int iCV,
	double T[][Nx], double R[][Nx], double Q[][Nx]
);

double thermalConductance(double Temperature, double h_TEM, double A_TEM);

double seebeckCoefficient(double temperature, double ZT, double resistivity, double lambda);

double electricCurrent(double Th, double Tc, int iCV);

double heatTransferCoeff(double Tg, double AR, double diaHyd);


double voltageOpenCircuit(double TemperatureHotSide, double TemperatureColdSide, double seebeckCoeff);

double electricPower(double electricCurrent, double resistanceLoad);

double energyTransferHotSide(double seebeckCoeff, double temperatureHotSide, double temperatureColdSide, double current, double internalResistance, double thermalConductance);

double energyTransferGas(double heatTransferCoeff, double temperatureBulk, double temperatureHotSide, double Rtotal);

double energyTransferRadiation(double hACoeff, double temperatureHotSide, double temperatureAmbient);

double hotsideTemperature
(
	double Tg, double Th, double Tc, double Tinf,
	double AreaCV, double AreaModule, double A_TE, double h_TEM,
	double AR, double diaHyd, double height,double tfin, double kfin, int Nf, double At, double Af,
	double epsInsulated, double epsModule
);

void outputVectors(double xi, double Rii,double Rli, double Kli, double ZTi,double Si, double Ii, double Pouti, double hgi,double cpgi,double QhotSidei, double Qgi, double Qci,double Tgi,double Thi,double Tci, double PoutTot, double QoutTot);

double compPropertyTable(int TE, int property, double T);

//double min(double a, double b);

double **readProp(const char *filename, int &p);

double interpProp(double **Prop, int N, double T);

double interpolatefromFile(const char *filename, double T);

double **Geom(const char *filename, int &N);

double computePD(double density, double Ain, double Abox, double mdot, int type);

// Function Declarations
	
	// The fin resistances includes contribution from the fin-base Area and extended fin surfaces
double Req( double Tgavg, double AR, double diaHyd, double height, double tfin, double kfin, int Nf, double At, double Af)
{
	double Lc, m, nf, no, Rtot_, hbar;

	hbar = heatTransferCoeff(Tgavg, AR, diaHyd);
	Lc = height/2 + tfin/2;
	m = sqrt(2.0*hbar/(kfin*tfin));
	nf = (tanh(m*Lc))/(m*Lc);
	no = 1.0 - Nf*Af*(1.0 - nf)/At;
	Rtot_ = 1.0/(no*1.0*hbar*At);
	return Rtot_;
}

// The contact resistance between fin and TEG
double Rcfin( double dummy)
{
	return 0.0;
}

// The resistance contribution by the TEG base material
double Rbasethickness( int iCV )
{
	return base.t/(base.k * base.Area[iCV] + SMALL);
}

// The contact resistance between TEM hot junction and hot surface.
double RcTEMH( int iCV)
{
	if(TEM.Type[iCV] == 0)
	{
		tgrease.gap = 0.001;
		tgrease.k = 5;}
	else if(TEM.Type[iCV] == 1)
	{
		tgrease.gap = 0.0001;
		tgrease.k = 2.30765;
	}
	double rc = tgrease.gap/(tgrease.k * tgrease.Area[iCV] + SMALL);	
	return rc;
}

// The resistance by insulation material (min K)
double Rinsulation(int iCV)
{
	return insulation.t/(insulation.k * insulation.Area[iCV] + SMALL);
}

void computeTECircuit(double Th, double Tc, int iCV)
{
	double alpha, Kl, Tavg, nt,i,v,p,r;
	Tavg = (Th + Tc)/2.0;

	updateTECouples(Tavg, iCV);	
	
	alpha = TEC.Alpha;
	Kl = TEC.Kl;
	r = TEC.R;
	
	nTE[iCV] = TEM.nTE[iCV];	// Number of TE couples in a CV
	nt = nTE[iCV];
	
	//Ri[iCV] = nTE[iCV]*TEC.R;		// Adding the internal resistances in series for electrical circuit
	Ri[iCV] = TEC.R;
	Rl[iCV] = Ri[iCV];
	
	// The following calculations are done by adding all the TEMs in series electrically and parallel thermally in a CV
	//V[iCV] = nTE[iCV]*alpha*(Th - Tc);	// Adding all the voltages generated due to deltaT
	V[iCV] = alpha*(Th - Tc);
	I[iCV] = V[iCV]/(Ri[iCV] + Rl[iCV]);	// Matching Internal Impedance
	//I[iCV] = SMALL;
	P[iCV] = I[iCV]*I[iCV]*Rl[iCV]; // Power by one TEC
	//P[iCV] = SMALL;
	i = I[iCV];
	v = V[iCV];
	p = P[iCV];

	//cout<<"\n NTE "<<nTE[iCV];
	//cout<<"\n V "<<V[iCV];
}

// The effective internal resistance of TEC as seen from hot side
double RintTEM(double Th, double Tc, int iCV)
{
	// Modeling Effective TEM Reistance
	// Hot side HT Qh = f(Th,Tc)
	double Qh,Pout,t1,t2,t3,i,alpha;	

	//cout<<"\n RintTEM "<<iCV<<"\t"<<Th<<"\t"<<Tc;
	
	// Updating TE Circuit
	computeTECircuit(Th,Tc,iCV);	
	
	i = I[iCV];
	alpha = TEC.Alpha;

	t1 = TEC.Alpha*I[iCV]*Th;
	t2 = -0.5*I[iCV]*I[iCV]*TEC.R;
	t3 = TEC.Kl*(Th - Tc);
	Qh = t1 + t2 + t3;	// Heat Transfer through one TE Couple
	Qh = Qh*nTE[iCV];	// Adding All ThermoElectric couples in one CV
	Pout = P[iCV]*nTE[iCV]; // Adding all power by ALl TECs in one CV
	return (Th - Tc)/(Qh - Pout + SMALL);	
}

// Computing Radiation resistance for the TEM Hot surface
double RradTEM(double Th, double Tinf, int iCV)
{
	double sigma, hArad, area;	
	sigma = 5.6703E-08;
	area = TEM.Area[iCV] - nTE[iCV]*TEC.Area + SMALL; // Area which is uncovered by TECs
	hArad = TEM.epsilon * area * sigma * ( pow(Th,3.0) + pow(Tinf,3.0) + (Th + Tinf)*Th*Tinf);
	return 1.0/hArad;
}

// Computing Resistances for convection through Insulation
double RconvIns(double beta, double alpha, double nu, double deltaT, double L, double Area)
{
	double Ra;
	double g = -9.81, hconv = 1e-10;
	double k;
		
	Ra = (g*beta)/(alpha*nu)*deltaT*pow(L,3);
	if(Ra >= 1e5 && Ra <= 2e7 )
		hconv = 0.54 * k * pow(Ra,0.25)/L;
	else if( Ra >= 2e7 && Ra <= 3e10 )
		hconv = 0.14 * k * pow(Ra,1/3)/L;

	return 1.0/(hconv*Area);

	//For hot plates facing down
	//	hconv11 = 0.27 * k * power(Ra,1/3)/L
}

// Computing Radiation resistance for the Insulation material 
double RradIns(double Tins, double Tinf, int iCV)
{
	double sigma, hArad;	
	sigma = 5.6703E-08;
	hArad = insulation.epsilon * insulation.Area[iCV] * sigma * ( pow(Tins,3.0) + pow(Tinf,3.0) + (Tins + Tinf)*Tins*Tinf);
	return 1.0/hArad;
}

// The contact resistance between TEM cold junction and coolant
double RcTEMC( int iCV)
{
	if(TEM.Type[iCV] == 0) // Grafoil
	{
		tgrease.gap = 0.001;
		tgrease.k = 5;
	}
	else if(TEM.Type[iCV] == 1) //Omegatherm
	{
		tgrease.gap = 0.0001;
		tgrease.k = 2.30765;
	}
	return tgrease.gap/(tgrease.k * tgrease.Area[iCV]);
}


double ReqLoad(double Th, double Tc, double Tcoolant, int iCV)
{
	double t1,t2,t3,pp;
	t1 = Th;
	t2 = Tc; 
	t3 = Tcoolant;
	computeTECircuit(Th, Tc, iCV);
	pp = nTE[iCV]*P[iCV];	// Total TEC power in one CV
	return (Th - Tcoolant)/(pp + SMALL);	
}


double fReCorr(double x)
{
    // interpolation on 8.1 (Incropera)for constant temperature
    if(x<=1.0)
        return 57;
    else if(x>1.0 && x <= 1.43)
         return 57.0 + (x - 1.0)*(59.0 - 57.0)/(1.43 - 1.0);                    
    else if(x>1.43 && x <= 2.0)
         return 59.0 + (x - 1.43)*(62.0 - 59.0)/(2.0 - 1.43);
    else if(x>2.0 && x <= 3.0)
         return 62.0 + (x - 2.0)*(69.0 - 62.0)/(3.0 - 2.0);
    else if(x>3.00 && x <= 4.0)
         return 69.0 + (x - 3.0)*(73.0 - 69.0)/(4.0 - 3.0);
    else if(x>4.0 && x <= 8.0)
         return 73.0 + (x - 4.0)*(82.0 - 73.0)/(8.0 - 4.0);
    else 
         return min(82.0 + (x - 8.0)*(96.0 - 82.0)/(1000), 96.0);
                             
}

// Need to be reviewed
// Look for correlations from Incropera book or do simulations in fluent to account for the values
double NuCorr(double AR)
{
    // cubic interpolation on 8.1 (Incropera) for constant temperature
    double p1 = 0.0009736, p2 = -0.040208, p3 = 0.67904, p4 = 2.2533;
	double x = AR;
	double val = p1*x*x*x + p2*x*x + p3*x + p4;
	if(val > 8.23)
	    return 8.23; 
    else 
    return val;
//    return min(p1*x*x*x + p2*x*x + p3*x + p4,8.23) ;
}

double heatTransferCoeff(double Tg, double AR, double diaHyd)
{
	double htf = air.k*NuCorr(AR)/diaHyd;
    return htf;
}

string cwd()
{
	string cwDir = _getcwd(NULL,0);
	return cwDir;
}

double min(double a, double b)
{
       if( a < b)
           return a;
       else
           return b;
}

// Function to Update Air Properties based on Average Bulk Temperature
void updateAir(double Tgavg)
{
	air.Cp = interpProp(air.CpTable,air.size, Tgavg);	
	air.k = interpProp(air.kTable,air.size, Tgavg);	
	air.density = interpProp(air.densityTable,air.size, Tgavg);	
	air.mu = interpProp(air.muTable,air.size, Tgavg);	
}

	//	Update TEcouple Parameters
//Function to update TE couples
void updateTECouples(double Tavg, int iCV)
{
	int i,j;
	double lambdan, lambdap;
	double rhon, rhop;
	double alphan, alphap;
	double ZTn,ZTp;
	double height, area, R, Kl, alpha;

	if(TEM.Type[iCV] == 0)	// Skutterudites - Use Material 0 an n-type and 1 as p-type
	{
		i = 0;
		j = 1;
		
		// Interpolating properties of n-type TE Materials for Tavg		
		lambdan = interpProp(TEMat[i].lambda,TEMat[i].lambdasize, Tavg);		// Thermal conductivity
		rhon = interpProp(TEMat[i].rho,TEMat[i].rhosize, Tavg);				    // Electrical resistivity
		alphan = interpProp(TEMat[i].alpha,TEMat[i].alphasize,Tavg);  			// Seebeck Coefficient
		ZTn = interpProp(TEMat[i].ZT,TEMat[i].ZTsize,Tavg);  					// Figure of Merit
		
		// Interpolating properties of p-type TE Materials for Tavg
		lambdap = interpProp(TEMat[j].lambda,TEMat[j].lambdasize, Tavg);		// Thermal conductivity
		rhop = interpProp(TEMat[j].rho,TEMat[j].rhosize, Tavg);				    // Electrical resistivity
		alphap = interpProp(TEMat[j].alpha,TEMat[j].alphasize,Tavg);  			// Seebeck Coefficient
		ZTp = interpProp(TEMat[j].ZT,TEMat[j].ZTsize,Tavg);  					// Figure of Merit
	}
	else if(TEM.Type[iCV] == 1)	// Tellurides - Use Material 2 an n-type and 3 as p-type
	{
		i = 2;
		j = 3;
		
		// Interpolating properties of n-type TE Materials for Tavg
		lambdan = interpProp(TEMat[i].lambda,TEMat[i].lambdasize, Tavg);		// Thermal conductivity
		rhon = interpProp(TEMat[i].rho,TEMat[i].rhosize, Tavg);				    // Electrical resistivity
		alphan = interpProp(TEMat[i].alpha,TEMat[i].alphasize,Tavg);  			// Seebeck Coefficient
		ZTn = interpProp(TEMat[i].ZT,TEMat[i].ZTsize,Tavg);  					// Figure of Merit
		
		// Interpolating properties of p-type TE Materials for Tavg
		lambdap = interpProp(TEMat[j].lambda,TEMat[j].lambdasize, Tavg);		// Thermal conductivity
		rhop = interpProp(TEMat[j].rho,TEMat[j].rhosize, Tavg);				    // Electrical resistivity
		alphap = interpProp(TEMat[j].alpha,TEMat[j].alphasize,Tavg);  			// Seebeck Coefficient
		ZTp = interpProp(TEMat[j].ZT,TEMat[j].ZTsize,Tavg);  					// Figure of Merit
	}
	
	TEC.height = TEM.height[TEM.Type[iCV]];// TEC Leg Height
	if(TEM.Type[iCV] == 0)
		TEC.Area = 2.0 * 0.004 * 0.004; // TEM.Area[iCV]*TEM.etaNP[TEM.Type[iCV]];  [4] Area of one leg x 2 for SKu
	else	
		TEC.Area = 2.0 * 0.002 * 0.002; // TEM.Area[iCV]*TEM.etaNP[TEM.Type[iCV]];  [4] Area of one leg x 2 for BiTE
	
	TEC.Alpha = abs(alphan) + abs(alphap);	// 
	TEC.R = rhon*TEC.height/(TEC.Area/2.0) + rhop*TEC.height/(TEC.Area/2.0);	// Internal Electrical Resistance
	// Adding TEC COntact Resistance = 0.53 ohm[4]
	TEC.R += 0.0;  //[0.108 for Sk, 0.025 for BiTE]


	TEC.Kl = lambdan*(TEC.Area/2.0)/TEC.height + lambdap*(TEC.Area/2.0)/TEC.height;	// Internal Electrical Resistance 
	TEC.ZT = pow((alphap - alphan),2.0)/pow((sqrt(lambdap*rhop) + sqrt(lambdan*rhon)),2)*Tavg; // Figure of Merit for TEC [Refer to Rowe, 3,3]
	ZT[iCV] = TEC.ZT;

	height = TEC.height;
	area = TEC.Area;
	R = TEC.R;
	Kl = TEC.Kl;
	alpha = TEC.Alpha;
	return;
}



// //Compute electric current & power produced by TE
double electricCurrent(double Th, double Tc, int iCV)
{	
	double alpha, Ri, Rl, I, nTE, V, Tavg;
	Tavg = (Th + Tc)/2.0;
	updateTECouples(Tavg, iCV);	
	alpha = TEC.Alpha;
	Ri = TEC.R;	
	nTE = TEM.nTE[iCV];	// Number of TE couples in a CV

	// Assumption is Rl is equal to Ri
	Rl = Ri;
	
	// The following calculations are done by adding all the TEMs in series electrically and parallel thermally in a CV
	V = nTE*alpha*(Th - Tc);	// Adding all the voltages generated due to deltaT
	I = V/(Ri + Rl);		
	return I;
}

// Open Circuit Voltage
double voltageOpenCircuit(double temperatureHotSide, double temperatureColdSide, double seebeckCoeff)
{
	return seebeckCoeff*(temperatureHotSide - temperatureColdSide);
}

double electricPower(double electricCurrent, double resistanceLoad)
{
	return electricCurrent*electricCurrent*resistanceLoad;
}

double energyTransferHotSide(double seebeckCoeff, double temperatureHotSide, double temperatureColdSide, double current, double internalResistance, double thermalConductance)
{
	return seebeckCoeff*temperatureHotSide*current - 0.5*current*current*internalResistance + thermalConductance*(temperatureHotSide - temperatureColdSide);
}

// Need to review equation
double energyTransferGas(double heatTransferCoeff, double temperatureBulk, double temperatureHotSide, double Rtotal)
{
	//return heatTransferCoeff*( temperatureBulk - temperatureHotSide)*Area;
	
	return ( temperatureBulk - temperatureHotSide)/(Rtotal);
}


double energyTransferRadiation(double hACoeff, double temperatureHotSide, double temperatureAmbient)
{
	return hACoeff*( temperatureHotSide - temperatureAmbient);
}

void outputVectors
(
	double xi, double Rii,double Rli, double Kli, 
	double ZTi,double Si, double Ii, double Pouti, 
	double hgi,double cpgi,double QhotSidei, double Qgi, double Qci,
	double Tgavgi,double Thi,double Tci,double QhTot, double PoutTot
)
{
        cout<<"\n Displaying Vectors for distance : "<<xi;
		cout<<"\n Resistance --> Internal : "<<Rii<<" Load "<<Rli;
		cout<<"\n Conductance "<<Kli<<" ZT "<<ZTi<<" S "<<Si;
		cout<<"\n I "<<Ii<<" Pout "<<Pouti;
		cout<<"\n hg "<<hgi<<" cpg "<<cpgi;
		cout<<"\n QhotSide "<<QhotSidei<<" Qg "<<Qgi<<" Qc "<<Qci;
		cout<<"\n PoutTot "<<PoutTot<<" QhTot "<<QhTot;
		cout<<"\n Tgavg "<<Tgavgi<<" Th "<<Thi<<" Tc "<<Tci;
}

	
// Solution of equation for hotsideTemperature. 
void Thot
(
	double Tgavg, double Tavg, double Tc,
	double AR, double diaHyd, double height,double tfin, double kfin, int Nf, double At, double Af,
	int iCV,
	double T[][Nx], double R[][Nx], double Q[][Nx]
)
{
	int ctr, i;
	double Kl, Ri, Rl, S, hg, hArad, Thold, Rtotal;
	double tol, error;
	double root1, root2, a, b, c;
	double Rhg, Rhc, Rhctt, Rhct, Rhcbm, Rhcbr, Rhcbl, Rhcb;
	double R01;
	double Th, Therr, terr, Rerr;
	double t[10], r[15];

	
	tol = 1e-6;
	ctr = 0;
	
	do
	{
		if(ctr == 0)  // For first iteration, Guess the Temperatures using Tavg, Tc, Tgavg
		{
			// Update Temperature on Guess Values			
			updateTECouples(Tavg, iCV);	//	Update TECouple		
			updateAir(Tgavg);	// Update Air
			T[0][iCV] =	Tgavg;	// Tg
			T[1][iCV] = (Tavg + 2.0*Tgavg)/3.0;
			T[2][iCV] = (Tavg + 2.0*Tgavg)/3.0;
			T[3][iCV] =	Tavg;
			T[4][iCV] =	(Tavg + Tgavg)/2.0;
			T[5][iCV] = (Tavg + Tc)/2.0;
			T[6][iCV] = (Tavg + Tc)/2.0;
			T[7][iCV] = Tc;

			R01 =	Req(Tgavg, AR, diaHyd, height, tfin, kfin, Nf, At, Af);		
		}
		else	// Update Temperatures on basis of Resitances using Th as Reference
		{
			T[0][iCV] =	Tgavg;	// Tg

			T[1][iCV] =	T[0][iCV]/R01 + T[3][iCV]/(R[2][iCV] + R[3][iCV] + SMALL);
			T[1][iCV] /= ( 1.0/R01 + 1.0/(R[2][iCV] + R[3][iCV] + SMALL) );

			T[2][iCV] = T[0][iCV]/(R01 + R[2][iCV] + SMALL) + T[3][iCV]/(R[3][iCV] + SMALL);
			T[2][iCV] /= ( 1.0/(R01 + R[2][iCV] + SMALL) + 1.0/(R[3][iCV] + SMALL) );

			//T[4][iCV] = T[3][iCV]/(R[4][iCV] + R[12][iCV] + SMALL) + T[7][iCV]/(Rhcb + SMALL);
			//T[4][iCV] /= ( 1.0/(R[4][iCV] + R[12][iCV] + SMALL) + 1.0/(Rhcb + SMALL));
			T[4][iCV] = T[3][iCV]/(Rhcbl + SMALL) + T[7][iCV]/(Rhcbr + SMALL);
			T[4][iCV] /= ( 1.0/(Rhcbl + SMALL) + 1.0/(Rhcbr + SMALL));
			
			//cout<<" \n"<<iCV<<" Resistance Value "<<R[1][iCV]<<"\t"<<R[2][iCV]<<"\t"<<R[3][iCV]<<"\t"<<R[4][iCV]<<"\t"<<R[5][iCV]<<"\t"<<R[6][iCV]<<"\t"<<R[7][iCV]<<"\t"<<R[8][iCV]<<"\t"<<R[9][iCV];
			//cout<<" \n"<<iCV<<" t "<<R[10][iCV]<<"\t"<<R[11][iCV]<<"\t"<<R[12][iCV]<<"\t"<<R[13][iCV];

			T[5][iCV] = T[3][iCV]/(R[5][iCV] + SMALL) + T[7][iCV]/(R[9][iCV] + SMALL);
			T[5][iCV] /= (1.0/(R[5][iCV] + SMALL) + 1.0/(R[9][iCV] + SMALL));

			// Writing in terms of 4 and 7
			T[6][iCV] = T[4][iCV]/(R[6][iCV] + SMALL) + T[7][iCV]/(R[10][iCV] + R[12][iCV] + SMALL);
			T[6][iCV] /= (1.0/( R[6][iCV] + SMALL) + 1.0/(R[10][iCV] + R[12][iCV] + SMALL));

			T[7][iCV] = Tc;

			T[8][iCV] = T[3][iCV]/(R[4][iCV] + SMALL) + T[4][iCV]/(R[12][iCV] + SMALL);
			T[8][iCV] /= ( 1.0/(R[4][iCV] + SMALL) + 1.0/(R[12][iCV] + SMALL));

			T[9][iCV] = T[4][iCV]/(R[12][iCV] + R[6][iCV] + SMALL) + T[7][iCV]/(R[10][iCV] + SMALL);
			T[9][iCV] /= ( 1.0/(R[12][iCV] + R[6][iCV] + SMALL) + 1.0/(R[10][iCV] + SMALL));
			//cout<<" \n"<<iCV<<" Thot "<<T[1][iCV]<<"\t"<<T[2][iCV]<<"\t"<<T[3][iCV]<<"\t"<<T[4][iCV]<<"\t"<<T[5][iCV]<<"\t"<<T[6][iCV]<<"\t"<<T[7][iCV]<<"\t"<<T[8][iCV]<<"\t"<<T[9][iCV];
			
			updateAir(Tgavg);	// Update Air
		}

/*		
		// Update All Temperatures and Current
		R[0][iCV] = 0.0;		// fin Resistance
		R[1][iCV] = 0.0;	// Fin - Base Resistance
		// R01 = R[0][iCV]*R[1][iCV]/(R[0][iCV] + R[1][iCV]);
		R01 =	Req( height, tfin, kfin,Nf, At, Af);
		R[2][iCV] =	Rcfin(0.0); 		// Fin-Contact Resistance
		R[3][iCV] =	Rbasethickness(iCV); 		// Base TEG Resistance
		R[4][iCV] = RcTEMH(iCV);	// Contact Resistance TEM - Th
		R[5][iCV] =	Rinsulation(iCV);	// Insulation Resistance
		R[6][iCV] =	RintTEM(T[4][iCV],T[6][iCV], iCV);	// Effective TEM Resistance = (T4-T6)/Qh(T4,T6)
		R[7][iCV] = RradTEM(T[4][iCV], T[7][iCV], iCV);	// TEM Radiation Resistance
//		R[8][iCV] =	RconvIns(beta, alpha, nu, deltaT, L, Area);// INS Convection Resistance
		R[8][iCV] = 0.0; // Assuming there is no convection
		R[9][iCV] =	RradIns(T[5][iCV], T[7][iCV], iCV);	// INS Radiation Resistance
		R[10][iCV] = RcTEMC(iCV); // Contact Resistance TEM - Th
		R[11][iCV] = ReqLoad( T[4][iCV],T[7][iCV], iCV);	// Equivalent Load Resistance
*/
	
		// Update All Resistances			
		R01 =	Req(Tgavg, AR, diaHyd, height, tfin, kfin, Nf, At, Af);
		R[1][iCV] = R01;
		R[2][iCV] =	Rcfin(0.0); 		// Fin-Contact Resistance
		R[3][iCV] =	Rbasethickness(iCV); 	// Base TEG Resistance
		R[4][iCV] = RcTEMH(iCV);		// Contact Resistance TEM - Th
		R[5][iCV] =	Rinsulation(iCV);		// Insulation Resistance		
		R[6][iCV] =	RintTEM(T[4][iCV],T[6][iCV], iCV);	// Effective TEM Resistance = (T4-T6)/Qh(T4,T6)
		R[7][iCV] = RradTEM(T[4][iCV], T[7][iCV], iCV);		// TEM Radiation Resistance		
		R[8][iCV] = LARGE;		// Assuming there is no convection
		R[9][iCV] =	RradIns(T[5][iCV], T[7][iCV], iCV);		// INS Radiation Resistance
		R[10][iCV] = RcTEMC(iCV); // Contact Resistance TEM - Th
		R[11][iCV] = ReqLoad( T[4][iCV], T[6][iCV], T[7][iCV], iCV);	// Equivalent Load Resistance
		
		// Computing Resistance Brances
		// Adding Ceramic plate thermal resistance
		// Alumina Plate -  1.5 mm, 20 W/m-K
		//double rc = tgrease.gap/(tgrease.k * tgrease.Area[iCV]);
		R[12][iCV] = 0.0015/(20.0*tgrease.Area[iCV] + SMALL);
		
	
		// Resistances between 0 and 3
		Rhg = R[3][iCV] + R[2][iCV] + R01; //R[0][iCV]*R[1][iCV]/((R[0][iCV] + R[1][iCV]); // Resistance between hot surface and Gas 

		// Resistances between top Branches
		// Insulation Part
		//Rhctt = R[9][iCV]; //1.0/(1.0/(R[8][iCV] + SMALL) + 1.0/(R[9][iCV] + SMALL));
		Rhct = R[9][iCV] +  R[5][iCV] ;	// Top branch resistance between coolant and hot surface

		// Resistances between bottom branches		
		Rhcbm = R[6][iCV] + R[10][iCV] + R[12][iCV]; // Medium Branch
		Rhcbr = 1.0/(1.0/(R[7][iCV] + SMALL) + 1.0/(Rhcbm + SMALL) + 1.0/(R[11][iCV] + SMALL));  // Bottom branch including TEM + Radiation
		Rhcbl = R[4][iCV] + R[12][iCV]; // Bottom excluding TEM
		Rhcb = Rhcbr + Rhcbl;

		Rhc = 1.0/(1.0/Rhct + 1.0/Rhcb + SMALL);

		// Update All the Temperatures
		// Update T3
		Th =  (T[0][iCV]/Rhg + T[7][iCV]/Rhc )/(1.0/(Rhg + SMALL) + 1.0/(Rhc + SMALL));
		Therr = abs( (Th - T[3][iCV])/(T[3][iCV] + SMALL));
		T[3][iCV] = Th;	
		//cout<<"\n Error " << Therr;
		//cout<<"\n "<<ctr<<"\t"<<iCV<<" Th "<<T[3][iCV];
		
		terr = 0;
		Rerr = 0;
		for(i = 0; i < 10; i++)
		{
			terr += abs(t[i] - T[i][iCV])/abs(T[i][iCV]);
			t[i] = T[i][iCV];			
		}
		for(i = 0; i < 13; i++)
		{
			if(i != 8 && i != 2)
			{
				Rerr += abs(r[i] - R[i][iCV])/abs(R[i][iCV]);
			}
			r[i] = R[i][iCV];
		}


		ctr++;	
	}
	//while(Therr > tol);	
	while(Rerr + terr > tol);
	//while(ctr < 100);
	/*
	{
		// Updating Temps
		T[0][iCV] =	Tgavg;	// Tg

		T[3][iCV] = (T[0][iCV]/Rhg + T[7][iCV]/Rhc )/(1.0/(Rhg + SMALL) + 1.0/(Rhc + SMALL));

		T[1][iCV] =	T[0][iCV]/R01 + T[3][iCV]/(R[2][iCV] + R[3][iCV] + SMALL);
		T[1][iCV] /= ( 1.0/R01 + 1.0/(R[2][iCV] + R[3][iCV] + SMALL) );

		T[2][iCV] = T[0][iCV]/(R01 + R[2][iCV] + SMALL) + T[3][iCV]/(R[3][iCV] + SMALL);
		T[2][iCV] /= ( 1.0/(R01 + R[2][iCV] + SMALL) + 1.0/(R[3][iCV] + SMALL) );

		T[4][iCV] = T[3][iCV]/(R[4][iCV] + R[12][iCV] + SMALL) + T[7][iCV]/(Rhcb + SMALL);
		T[4][iCV] /= ( 1.0/(R[4][iCV] + R[12][iCV] + SMALL) + 1.0/(Rhcb + SMALL));

		T[5][iCV] = T[3][iCV]/(R[5][iCV] + SMALL) + T[7][iCV]/(Rhct + SMALL);
		T[5][iCV] /= (1.0/(R[5][iCV] + SMALL) + 1.0/(Rhct + SMALL));

		// Writing in terms of 4 and 7
		T[6][iCV] = T[4][iCV]/(R[6][iCV] + SMALL) + T[7][iCV]/(R[10][iCV] + R[12][iCV] + SMALL);
		T[6][iCV] /= (1.0/( R[6][iCV] + SMALL) + 1.0/(R[10][iCV] + R[12][iCV] + SMALL));

		T[7][iCV] = Tc;

		T[8][iCV] = T[3][iCV]/(R[4][iCV] + SMALL) + T[4][iCV]/(R[12][iCV] + SMALL);
		T[8][iCV] /= ( 1.0/(R[4][iCV] + SMALL) + 1.0/(R[12][iCV] + SMALL));

		T[9][iCV] = T[6][iCV]/(R[12][iCV] + SMALL) + T[7][iCV]/(R[10][iCV] + SMALL);
		T[9][iCV] /= ( 1.0/(R[12][iCV] + SMALL) + 1.0/(R[10][iCV] + SMALL));

		// Update All Resistances			
		R01 =	Req(Tgavg, AR, diaHyd, height, tfin, kfin, Nf, At, Af);
		R[1][iCV] = R01;
		R[2][iCV] =	Rcfin(0.0); 		// Fin-Contact Resistance
		R[3][iCV] =	Rbasethickness(iCV); 	// Base TEG Resistance
		R[4][iCV] = RcTEMH(iCV);		// Contact Resistance TEM - Th
		R[5][iCV] =	Rinsulation(iCV);		// Insulation Resistance
		R[6][iCV] =	RintTEM(T[4][iCV],T[6][iCV], iCV);	// Effective TEM Resistance = (T4-T6)/Qh(T4,T6)
		R[7][iCV] = RradTEM(T[4][iCV], T[7][iCV], iCV);		// TEM Radiation Resistance		
		R[8][iCV] = LARGE;		// Assuming there is no convection
		R[9][iCV] =	RradIns(T[5][iCV], T[7][iCV], iCV);		// INS Radiation Resistance
		R[10][iCV] = RcTEMC(iCV); // Contact Resistance TEM - Th
		R[11][iCV] = ReqLoad( T[4][iCV], T[6][iCV], T[7][iCV], iCV);	// Equivalent Load Resistance
		
		// Computing Resistance Brances
		// Adding Ceramic plate thermal resistance
		// Alumina Plate -  0.5 mm, 20 W/m-K
		//double rc = tgrease.gap/(tgrease.k * tgrease.Area[iCV]);
		R[12][iCV] = 0.0005/(20.0*tgrease.Area[iCV]);
	
	}
	*/
	// After Temperatures & Resistances have converged, compute the Heat Fluxes
	// These are total fluxes in a CV
	//Q[0][iCV] =	( T[0][iCV] - T[3][iCV] )/(R01 + R[3][iCV] + R[2][iCV] + SMALL);
	Q[0][iCV] =	( T[0][iCV] - T[3][iCV] )/(Rhg + SMALL);
	//Q[1][iCV] =	( T[3][iCV] - T[5][iCV] )/(R[5][iCV] + SMALL);		
	Q[1][iCV] =	( T[3][iCV] - T[7][iCV] )/(Rhct + SMALL);		
	//Q[2][iCV] =	( T[3][iCV] - T[4][iCV] )/(R[4][iCV] + R[12][iCV] + SMALL);
	Q[2][iCV] =	( T[3][iCV] - T[7][iCV] )/(Rhcb + SMALL);
	//Q[1][iCV] = Q[0][iCV] - Q[2][iCV];
	//Q[3][iCV] =	Q[1][iCV]; //
	Q[3][iCV] = ( T[5][iCV] - T[7][iCV] )/(R[9][iCV] + SMALL);
	Q[4][iCV] =	0.0; //( T[5][iCV] - T[7][iCV] )/(R[8][iCV] + SMALL);
	Q[5][iCV] =	( T[4][iCV] - T[7][iCV] )/(R[7][iCV] + SMALL);
	//Q[6][iCV] =	( T[4][iCV] - T[6][iCV] )/(R[6][iCV] + SMALL);   // Qc as seen from hot side. Effective Resistance = (T4-T6)/Qh(T4,T6)
	Q[6][iCV] =	( T[4][iCV] - T[7][iCV] )/(Rhcbm + SMALL);   // Qc as seen from hot side. Effective Resistance = (T4-T6)/Qh(T4,T6)
	Q[7][iCV] = ( T[4][iCV] - T[7][iCV] )/(R[11][iCV] + SMALL);
	//P[iCV] = Q[7][iCV];

	return;
}

double interpolatefromFile(const char *filename, double T)
{
	ifstream ifs;
	int i = 0, N;
	double x[100],y[100];
  
	//cout<<"\n FIle name " << filename;

	ifs.open(filename);
	if ( ! ifs.fail() )
	{
		do 
		{
			ifs >> x[i] >> y[i];
			i++;
		}
		while (!ifs.eof());
		ifs.close();
	} 
	else if (ifs.fail())
	{
		cout << "Opening of txt file "<<filename<<" failed ! \n";
		cin.get();
		exit(1);
	}
	
	//Computing size of array
	N  = i;
	
	// Temperature value is in absolute scale
	
		
	if ( T > x[0] && T < x[N-1])
	{
		for( i = 0; i < N-1; i++)
		{
			if(T > x[i] && T <= x[i+1] )
				return y[i] + (y[i+1] - y[i])/(x[i+1] - x[i])*(T - x[i]);
		}
	}
	else if (T >= x[N-1])
		return y[N-1];
	else if (T <= x[0])
		return y[0];
	else 
	{
		cout<<"\n Interpolation error in text file "<<filename<<endl;
		cin.get();
		exit(1);
	}
}

double interpProp(double **Prop, int N, double T)
{
	double x[100],y[100];
	for(int i = 0;i < N;i++)
	{
		x[i] = Prop[i][0];
		y[i] = Prop[i][1];
	}
	
	// Temperature value is in absolute scale
	if ( T > x[0] && T < x[N-1])
	{
		for(int i = 0; i < N-1; i++)
		{
			if(T > x[i] && T <= x[i+1] )
				return y[i] + (y[i+1] - y[i])/(x[i+1] - x[i])*(T - x[i]);
		}
	}
	else if (T >= x[N-1])
		return y[N-1];
	else if (T <= x[0])
		return y[0];
	else 
	{
		cout<<"\n Interpolation error with T value "<<T<<endl;
		cin.get();
		exit(1);
	}
}

double **readProp(const char *filename, int &p)
{
	ifstream ifs;
	int i = 0;
	double **var;
	

	var = new double*[500];
    for( int j = 0 ; j < 500 ; j++ )
		var[j] = new double [2];

	ifs.open(filename);
	if ( ! ifs.fail() )
	{
		do 
		{			
			ifs >> var[i][0] >> var[i][1];
			i++;
		}
		while (!ifs.eof());
		ifs.close();
	} 
	else if (ifs.fail())
	{
		cout << "Opening of txt file "<<filename<<" failed ! \n";
		cin.get();
		exit(1);
	}
	
	//Computing size of array
	p = i;
	
	return var;
}

double **Geom(const char *filename, int &p)
{
	ifstream ifs;
	int i = 0;
	double **geom;
	

	geom = new double*[500];
    for( int j = 0 ; j < 500 ; j++ )
		geom[j] = new double [3];

	ifs.open(filename);
	if ( ! ifs.fail() )
	{
		do 
		{			
			ifs >> geom[i][0] >> geom[i][1] >> geom[i][2];
			i++;
		}
		while (!ifs.eof());
		ifs.close();
	} 
	else if (ifs.fail())
	{
		cout << "Opening of txt file "<<filename<<" failed ! \n";
		cin.get();
		exit(1);
	}
	
	//Computing size of array
	p = i;
	
	return geom;
}

// Cpmoutes Pressure Drops. Drops are +ve values
double computePD(double density, double Ain, double Abox, double mdot, int type)
{
	double AR, PD, vpipe, vbox, v1, mu;

	vpipe = mdot/(density*Ain);
	vbox =  mdot/(density*Abox);

	if(type == 0 ) // Inlet
	{
		AR = Ain/Abox;	
		v1 = vpipe;
	}
	else
	{
		v1 = vbox;
		AR = Abox/Ain;		
	}
	mu = 0.63 + 0.37*pow((1.0/AR),3.0);
	
	// Check  if it is expansion or contraction
	// Assuming sudden Expansion for area ratios > 4.0 
	// Using Boarda - Carnot Equation 
	if( AR < 1.0 ) // Diffuser .. Expansion   Matches well with tht of conical diffuser thing
	{
		PD = density*AR*(1.0 - AR)*pow(v1,2.0); // Sudden expansion in Pipe 
	}
	else if (AR > 1.0 ) // Nozzle ... Contraction
	{
		PD = 0.5*density*pow((1.0/mu - 1.0),2)*pow(AR,2)*pow(v1,2.0); // This is always a head loss
	}
	else
		PD = 0.0; // If AR = 0

	return PD;
}
				







/**************************************************************************************************************/

}