	
	// Setup Bulk and coolant temperatures
	Tc[iCV] = Tcinit;	// Constant coolant temp 	
	Th[iCV] = Tg[iCV];	// Guess value for Th
	Tg[iCV+1] = Tg[iCV];  //Guess value of Tg -- Tg represents Temperatures at the boundaries of each CVs				
	
	itr = 0;
	error = 100;

	// Set TEM Type SKU on basis of average TE Temperature as guess
	TEM.Type[iCV] = 0;
	TEM.nTE[iCV] = (TEM.numDensity[TEM.Type[iCV]]*TEM.Area[iCV]);
	Tavg = Tg[iCV];	// Guess value of Tavreage

	while((error > tol) && (++itr < itrmax))
	{
		//cout << "\n" << itr <<" " <<i;
		
		Tgavg = (Tg[iCV] + Tg[iCV+1])/2.0;	// Mean bulk Temperature of gas (Air)		
		
		if(itr == 1)
			Tavg = (Tg[iCV] + Tc[iCV])/2.0;
		else
			Tavg = (T[4][iCV] + T[6][iCV])/2.0; // Average TE Temperature

		// Set TEM Type depending on the Average TEM
		
		// Write a logic for change in Module Type after checking the Tavg and ZT curve
		if(Tavg <= 550.0 ) // Use BiTE Modules
		{
			TEM.Type[iCV] = 0;
			sCV = iCV;
		}
		else	// Use SK Modules
		{
			TEM.Type[iCV] = 0;
		}

		TEM.Area[iCV] = AreaModule[iCV];
		TEM.nTE[iCV] = (TEM.numDensity[TEM.Type[iCV]]*TEM.Area[iCV]);	// Number of TEM couples in the CV
		nTE[iCV] = TEM.nTE[iCV];
		
		// Write a logic to control total number of TEMs
		// Check if total number of thermocouples till last CV is more than 50
		// Technique 1
		nTEM0 = (nTETotcum0/32.0);
		nTEM1 = (nTETotcum1/127.0);
		
		/*if(nTEM0 + nTEM1 > 50)
		{
			TEM.Area[iCV] = 0.0;
			AreaModule[iCV] = 0.0;
			AreaInsulated[iCV] = AreaCV[iCV] - AreaModule[iCV]; // Area covered by Insulated Material		
			base.Area[iCV] = AreaCV[iCV];
			tgrease.Area[iCV] = AreaModule[iCV]; // Grease is only applied on Modules surfaces
			insulation.Area[iCV] = AreaInsulated[iCV];
			TEM.nTE[iCV] = 0.0;
			nTE[iCV] = 0.0;
		}
		*/

		/*
		// Technique2
		if(TEM.Type[iCV] == 0)
		{
			if(50 - nTEM0  > 0) // if there is possibility to add modules
			{
				nTE[iCV] = int(25.0*TEM.nC[0]/Ncv + 0.26);
				TEM.nTE[iCV] = nTE[iCV];
				TEM.Area[iCV] = nTE[iCV]/TEM.numDensity[0];
				AreaModule[iCV] = TEM.Area[iCV];
				AreaInsulated[iCV] = AreaCV[iCV] - AreaModule[iCV]; // Area covered by Insulated Material		
				base.Area[iCV] = AreaCV[iCV];
				tgrease.Area[iCV] = AreaModule[iCV]; // Grease is only applied on Modules surfaces
				insulation.Area[iCV] = AreaInsulated[iCV];
			}
			else
			{
				TEM.Area[iCV] = 0.0;
				AreaModule[iCV] = 0.0;
				AreaInsulated[iCV] = AreaCV[iCV] - AreaModule[iCV]; // Area covered by Insulated Material		
				base.Area[iCV] = AreaCV[iCV];
				tgrease.Area[iCV] = AreaModule[iCV]; // Grease is only applied on Modules surfaces
				insulation.Area[iCV] = AreaInsulated[iCV];
				TEM.nTE[iCV] = 0;
				nTE[iCV] = 0;
			}

		}
		else if(TEM.Type[iCV] == 1)
		{
			if(50 - nTEM0 > 0) // if there is possibility to add modules
			{
				int j = 0;
				while(!TEM.Type[j++]);
				nTE[iCV] = (((50 - nTEM0)/2)*TEM.nC[1]/(Ncv-j) + 0.26);
				//nTE[iCV] = (((50 - nTEM0)/2)*TEM.nC[1]/(Ncv-iCV) + 0.26);
				TEM.nTE[iCV] = nTE[iCV];
				TEM.Area[iCV] = nTE[iCV]/TEM.numDensity[1];
				AreaModule[iCV] = TEM.Area[iCV];
				AreaInsulated[iCV] = AreaCV[iCV] - AreaModule[iCV]; // Area covered by Insulated Material		
				base.Area[iCV] = AreaCV[iCV];
				tgrease.Area[iCV] = AreaModule[iCV]; // Grease is only applied on Modules surfaces
				insulation.Area[iCV] = AreaInsulated[iCV];
			}
			else
			{
				TEM.Area[iCV] = 0.0;
				AreaModule[iCV] = 0.0;
				AreaInsulated[iCV] = AreaCV[iCV] - AreaModule[iCV]; // Area covered by Insulated Material		
				base.Area[iCV] = AreaCV[iCV];
				tgrease.Area[iCV] = AreaModule[iCV]; // Grease is only applied on Modules surfaces
				insulation.Area[iCV] = AreaInsulated[iCV];
				TEM.nTE[iCV] = 0;
				nTE[iCV] = 0;
			}
		}
		*/
				
		//cout<<"\n NCV, Tavg, NTE "<<iCV<<"\t"<<Tavg<<"\t"<<TEM.nTE[iCV];
		//cin.get();
			
		// Solve for Th or T4
		Thot(Tgavg, Tavg, Tc[iCV], AR, diaHyd, height, tfin, kfin, Nf, At, Af, iCV, T, R, Q);
		Th[iCV] = T[3][iCV];
		
		// Update Air Properties
		updateAir(Tgavg);
				
		Tgnew = Tg[iCV] -  Q[0][iCV]/((mdot/2.0)*air.Cp);	// Only half side is being modeled
		error = fabs(Tgnew - Tg[iCV+1]);				
		//cout<<"\n "<<itr<<" Tg : old = " << Tg[iCV+1] <<" new = " << Tgnew;
		Tg[iCV+1] = Tgnew;
	}
		
	// After convergence of Tg. Recalculate Parameters
	Tgavg = (Tg[iCV] + Tg[iCV+1])/2.0;