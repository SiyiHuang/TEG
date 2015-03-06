	
	// Parameters per CV
	Pout_[iCV] = Q[7][iCV];	// Power output and heat transfer from TE module in one CV
	Qg_[iCV] = Q[0][iCV];		// Total Heat absorbed by the setup from the hot gas
	Qrad_[iCV] = Q[5][iCV] + Q[3][iCV];	// Total Radiation loss - from Insulation and exposed Module surface
	Qconv_[iCV] = Q[4][iCV];	// Convection loss through Insulation surface
	Qc_[iCV] = Q[6][iCV];	// Total heat lost to Coolant
	Qhtem_[iCV] = Q[6][iCV] + Q[7][iCV]; // Heat Transfer through TEM from hot side
	
	// Twice because it appears on both the sides
	PoutTot_[iCV] = 2.0*Pout_[iCV]; 
	QgTot_[iCV] = 2.0 * Qg_[iCV];
	QradTot_[iCV] = 2.0 * Qrad_[iCV];
	QconvTot_[iCV] = 2.0 * Qconv_[iCV];	
	QcTot_[iCV] = 2.0 * Qc_[iCV];
	QhtemTot_[iCV] = 2.0*Qhtem_[iCV];
	if(TEM.Type[iCV] == 0 )  // SKU
	{
		nTETot_[iCV][0] = 2.0*nTE[iCV];
		nTETot_[iCV][1] = 0.0;
	}
	else // BiTE
	{
		nTETot_[iCV][0] = 0.0;
		nTETot_[iCV][1] = 2.0*nTE[iCV];
	}
	
	// Define a strategy to define the Open Voltage Calculation
	//Vopenc_[iCV] = voltageOpenCircuit(Th[iCV], Tc[iCV], S_[iCV]);
	Vopenc_[iCV] = V[iCV]; // Voltage of One TEC
	ZT_[iCV] = ZT[iCV];
	
	// Pressure Drop Calculations
	Tgavg = T[0][iCV];
	updateAir(Tgavg); // Update Air Properties

	Kc = 0.45;	// Contraction Coefficient
	Ke = 0.0;	//Expansion Coefficient
	ReDh_= mch*diaHyd/(Ac*air.mu);
	f = fReCorr(AR)/ReDh_;
	delP_[iCV] = (f*dx/diaHyd)*0.5*air.density*(mch/(Ac*air.density))*(mch/(Ac*air.density));
	//cout<<"\n Re "<<mdot<<"\t"<<ReDh_<<"\t"<<mch/(Ac*air.density)<<"\t"<<finSpacing<<"\t"<<tfin;
						
	// Computing cumulative parameters
	if(iCV == 0)
	{
		PoutTotcum = PoutTot_[iCV];
		QgTotcum = QgTot_[iCV];
		QradTotcum = QradTot_[iCV];
		QconvTotcum = QconvTot_[iCV];
		QcTotcum = QcTot_[iCV];
		QhtemTotcum = QhtemTot_[iCV];
		nTETotcum0 = nTETot_[iCV][0];
		nTETotcum1 = nTETot_[iCV][1];
				
		// Pressure drops due to expansion and contraction to be added only once for channel		
		PressDropTot = delP_[iCV] + ( Kc + Ke )*0.5*(mch/Ac)*(mch/Ac)/air.density;
		
		// Also add expansion Loss from Pipe to channel and contraction loss
		// Using borda carnot equation for pipes. assuming rectangular secction of be equivalent pipe dia
		//PDin = air.density*ARin*(1.0 - ARin)*(mdot/(Ain*air.density))*(mdot/(Ain*air.density));
		PDin = computePD(air.density,Ain,Abox,mdot,0);		
		//PressDropTot =
	}
	else
	{
		PoutTotcum += PoutTot_[iCV];
		QgTotcum += QgTot_[iCV];
		QradTotcum += QradTot_[iCV];
		QconvTotcum += QconvTot_[iCV];
		QcTotcum += QcTot_[iCV];
		QhtemTotcum += QhtemTot_[iCV];
		PressDropTot += delP_[iCV];		
		nTETotcum0 += nTETot_[iCV][0];
		nTETotcum1 += nTETot_[iCV][1];
	}
	

	//cout<<"\n Heat Trf "<<heatTransferCoeff(Tgavg, AR, diaHyd);