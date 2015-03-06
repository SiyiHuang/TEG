	
	// Initializinf 'x'		
	
		
	dx = 0.01; // Size set to include one TEC Leg
	Ncv = int(flowLength/dx);
	dx = (flowLength/Ncv);
	x[0] = 0.0; // Inlet	

	Abox = width*height;
	ARin = Ain/Abox;
	ARout = Abox/Aout;
	
	//Set up Area Progression	
	for(i = 0; i < Ncv; i++)
	{		
		x[i+1] = x[i] + dx;	
		AreaCV[i] = width*dx;   // Area covered by one CV
		AreaModule[i] = eta_module*AreaCV[i]; // Area covered by TE Modules
		//AreaTE[i] = eta_TE*AreaModule[i];	// Area covered by all TE couples or unicouples
		AreaInsulated[i] = AreaCV[i] - AreaModule[i]; // Area covered by Insulated Material		
		base.Area[i] = AreaCV[i];
		tgrease.Area[i] = AreaModule[i]; // Grease is only applied on Modules surfaces
		insulation.Area[i] = AreaInsulated[i];
	}
		
	//Check if available area for TEM exceeds the Maximum Area allowed 
	/*	
	// Choose minimum of 0.8 or ratio
	if( Amodule_Max < 0.8*(width*flowLength))
		{
			eta_module = Amodule_Max/(width*flowLength);
		}
		else
		{
			eta_module = 0.8;
		}
	*/
	
	
	
	
		