/**************************************************************************************************************************************************************************/
/* Code deveoped by Sumeet on 17th June, 2011 to compute variation in cross-section for the TEG's
/
/ Solution methodology (As explained by Prof Heister by email )
/ 1)  Gather inputs on TE properties
/ 2) Setup the flow path Ag = Ag(x) and n = n(x)
/ 3)  Get the hg (convective heat transfer coefficienct ) from appropriate correlations / (option 
      to run fluent to get the value) / Make sure correlations incorporate effects of the developing boundary layer
/ 4)  Set the gas inflow conditions. Tg = Tg(i) at x = 0. It should include the Cp variation as the deltaT is pretty large.
/ 5) Step along the flowpath solving equations 1 - 10 sequentially and storing the results in vectors for display . Should 
     also keep the track of the open circuit voltage Voc =  S(Th - Tc)
 /****************************************************************************************************************************************************************************/

// Generating functions for the equations 1 - 10

// Function prototype declarations

// Assumption:-
// 1) rho, Kl, Ri are constant values
// 2) ZT versus T(hot) is available so S = S(T) is available
// 3) Heat transfer coefficient (hg)   
// Assuming constant properties except seebeck coefficient S = S(T)....

// %6th December%
// Revision to multiple functionalities 


#include <conio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <direct.h>

using namespace std;

#include "TEG.H"

int main()
{
	using namespace TEG;
	
	// Define Variables

	# include "defineVars.H"

	//Reset Variables including Global variables
	# include "resetVars.H"

	// User Inputs
	# include "userInputs.H"
	
	// Read Directory Paths
	# include "readPaths.H"

	// Read Geometries configurations from file
	string geomFile = currDir + "\\..\\..\\..\\Source\\TEMs.txt";
	geometry = Geom(geomFile.c_str(), geosize);

	// Read allowable Pressure Drop vs mass flow rate [GM]
	string mpdFile = currDir + "\\..\\..\\..\\Props\\Mass-PD.txt";
	
	// Read Material Property Data and TEM
	#include "readMaterials.H"
	
	string test = outputDir + "\\test.txt";
	ofstream outtest(test.c_str(),ios_base::out);
	//for(mdot = 0.02; mdot <= 0.1; mdot = mdot + 0.02)
	//for(Tginit = 400 + 273.15; Tginit <= 700 + 273.16; Tginit += 50)
	{
		
	ostringstream ms;
	ms << mdot;				
	string opt = outputDir + "\\Optimum_" + ms.str() + ".txt";
	ofstream outopt(opt.c_str(),ios_base::out);

	// Space for Output
	outopt <<"Index\th\tw\tl\tmdot\tPower\tTgd\tQg\tQr\tQh\tQc\tPDTotal\tPDin\tPDout\tAllwPD\tNTEC\ttf\tfs\tNf";
	
	// Looping over all geometry configurations
	int gi = 374;
	//for( gi = 0,index = 0; gi < geosize;  gi++)
	{
		// Computing (w,l,h) for each geometry configuration

		// Inputs for Basic Configuration
		width = widthB;
		height = heightB;
		flowLength = flowLengthB;
	
		finSpacing = finSpacingB;
		tfin = tfinB;
		


		//width = geometry[gi][0];
		//flowLength = geometry[gi][1];
		//height = geometry[gi][2];
		
		
		Powermax[index] = -1000;
		h[index] = height;
		w[index] = width;
		fl[index] = flowLength;						
		mdt[index] = mdot;
		

		// For Output file for Varying flowLength, width and height
		
		ostringstream ss;
		ss << gi << "_" << mdot << "_" <<Tginit;
	
		string finaloutput = outputDir +"\\Output_" + ss.str() + ".txt";
		
		ofstream outff(finaloutput.c_str(),ios_base::out);
		outff	<<"No\theight\twidth\tlength\tmdot\tTgDrop\tQgas\tQhTEM\tQrad\tQc\tPout\tPTot\tPDin\tPDout\tnTESK\tnTEBiTE\tNf\ttfin\tfinSpacing"<<endl;
			
		cout<<"\n gi "<<gi<<" mdot= "<<mdot<<" Height = " << height << " width = " << width << " Length " << flowLength;
		
    			
		#include "computeArea.H"
		
		// Looping over fin-spacing and fin-width for one-configuration
		//for( tfin = 0.001; tfin <= 0.008; tfin+=0.001)
		{
			//for ( finSpacing = 0.001; finSpacing <= 0.008; finSpacing += 0.001)
			//for ( Nf = 10; Nf < 151; Nf= Nf+10)
			{			
				
				#include "computeFin.H"	// Computing Parameters based on fin
				
				Ac = finSpacing*height;	// Cross-Sectional Area in the Channel
				diaHyd = 4.0 * Ac / (2.0 * (height + finSpacing));	// Hydraulic Diameter in channel
				AR = finSpacing/height; //Aspect Ratio for flow channel
									
				mch = mdot/Nf;	// Mass Flow Rate per channel
				
				// Setup inflow gas conditions at x = 0
				Tg[0] = Tginit;
				
				// cout<<"\n ***********************************************************************";
				// Traverse through all the CVs

				for(iCV = 0; iCV < Ncv; iCV++)
				{															
					//Compute Tg at the end of the each CV
					#include "solveTG.H"

					//cout<<"\n iCV "<<iCV<<"\t"<<Tg[iCV]<<"\t"<<Tg[iCV+1];
					//cin.get();
				
					// After All Temperatures are converged
					
					// Compute and Save Parameters of Interests per CV
					#include "computePOI.H"				

				}	// Ncv Loop ends here	

				// Adding pressure drop due to contraction near the exit
				muarea = 0.63 + 0.37*pow(1.0/ARout,3);
				PDout = computePD(air.density, Ain, Abox, mdot, 1);
				//PDout = -0.5*air.density*pow((1.0/muarea - 1.0),2)*pow(ARout,2)*pow(mdot/(Abox*air.density),2); // -ve sign as this will gain pressure
				PressDropTot = PressDropTot + PDin + PDout;

				/*
				// Test Output
				outtest << "\n T " ;
				for(iCV = 0; iCV < Ncv; iCV++){	for(j = 0; j <= 7; j++)	outtest << "\t" << T[j][iCV]; outtest <<"\n";}
				outtest << "\n R " ;
				for(iCV = 0; iCV < Ncv; iCV++){	for(j= 0; j <= 11; j++)	outtest << "\t" << R[j][iCV]; outtest <<"\n";}
				outtest << "\n Q " ;
				for(iCV = 0; iCV < Ncv; iCV++){	for(j = 0; j <= 8; j++)	outtest << "\t" << Q[j][iCV]; outtest << "\n";}
				outtest << "\n I\tV\tP\tRi\tRl\tnTE" ;
				for(iCV = 0; iCV < Ncv; iCV++)	outtest << "\n" << I[iCV]<<"\t"<<V[iCV]<<"\t"<<P[iCV]<<"\t"<<Ri[iCV]<<"\t"<<Rl[iCV]<<"\t"<<nTE[iCV];		
				outtest <<"\n Pressure Drops: Total\tIN\nOUT";
				outtest << "\n"<<PressDropTot<<"\t"<<PDin<<"\t"<<PDout<<endl;
				*/
				
				
				/*
				// For output per a fixed spacing and fin thickness
				ostringstream s;
				s << tfin << "_" << finSpacing;


				string outputfile = outputDir + "\\output_fin_" + s.str() + ".txt";
				ofstream outf(outputfile.c_str(),ios_base::out);

				// Output TE Parameters
				string matTEfile = outputDir + "\\matTE_" + s.str() + ".txt";
				ofstream matTEf(matTEfile.c_str(),ios_base::out);
				
				outf<<"I\tx\tTg\tTh\tQg\tQrad\tQh\tQc\tPout\tVopen\tDelP\tNf"<<endl;
				matTEf<<"I\tx\tRi\tRl\tKl\tZT\tS\tRtotfin"<<endl;
					
				for( i = 0 ; i < Ncv; i++)
				{
					outf	<<i<<"\t"<<(x[i]+x[i+1])/2.0<<"\t"
							<<T[0][i]<<"\t"<<Th[i]<<"\t"<<QgTot_[i]<<"\t"<<QradTot_[i]<<"\t"<<QhTot_[i]<<"\t"<<QcTot_[i]<<"\t"<<PoutTot_[i]
							<<"\t"<<Vopenc_[i]<<"\t"<<delP_[i]<<"\t"<<Nf<<"\t"<<gasSide[i]<<"\t"<<TESide[i]<<endl;
							
					matTEf	<<i<<"\t"<<(x[i]+x[i+1])/2.0<<"\t"
							<<Ri_[i]<<"\t"<<Rl_[i]<<"\t"<<Kl_[i]<<"\t"<<ZT_[i]<<"\t"<<S_[i]<<"\t"<<Rtot_[i]<<endl;
				}
				outf.close();
				matTEf.close();
				
					
				//cout<<"\n Pressure Drop Total "<< PressDropTot;
				
				// For final output files
				*/

				// Test Output
				
				for(iCV = 0; iCV < Ncv; iCV++){	outtest << "\t" << (x[iCV] + x[iCV+1])/2.0; outtest <<"\n";}
				outtest << "\n T " ;
				for(iCV = 0; iCV < Ncv; iCV++){	for(j = 0; j < 10; j++)	outtest << "\t" << T[j][iCV]; outtest <<"\n";}
				outtest << "\n R " ;
				for(iCV = 0; iCV < Ncv; iCV++){	for(j= 0; j < 13; j++)	outtest << "\t" << R[j][iCV]; outtest <<"\n";}
				outtest << "\n Q " ;
				for(iCV = 0; iCV < Ncv; iCV++){	for(j = 0; j < 8; j++)	outtest << "\t" << Q[j][iCV]; outtest << "\n";}
				outtest << "\n I\tV\tP\tRi\tRl\tnTE\tZT" ;
				for(iCV = 0; iCV < Ncv; iCV++)	outtest << "\n" << I[iCV]<<"\t"<<V[iCV]<<"\t"<<P[iCV]<<"\t"<<Ri[iCV]<<"\t"<<Rl[iCV]<<"\t"<<nTE[iCV]<<"\t"<<ZT[iCV]<<"\t"<<TEM.Type[iCV];		
				outtest <<"\n Pressure Drops: Total\tIN\nOUT";
				outtest << "\n"<<PressDropTot<<"\t"<<PDin<<"\t"<<PDout<<endl;
				
				
				outff	<<gi<<"\t"<<height<<"\t"<<width<<"\t"<<flowLength<<"\t"<<mdot
						<<"\t"<<(Tg[Ncv]-Tg[0])<<"\t"<<QgTotcum<<"\t"<<QhtemTotcum<<"\t"<<QradTotcum<<"\t"<<QcTotcum<<"\t"<<PoutTotcum
						<<"\t"<<PressDropTot<<"\t"<<PDin<<"\t"<<PDout
						<<"\t"<<nTETotcum0<<"\t"<<nTETotcum1<<"\t"<<Nf<<"\t"<<tfin<<"\t"<<finSpacing<<endl;
					 
				

				// Updating arrays for each geometry point...				
				AllowablePD = interpolatefromFile(mpdFile.c_str(), mdot);	
				//cout<<" \n pd " << AllowablePD;
				if(PressDropTot <= AllowablePD)
				{
					if(Powermax[index] < PoutTotcum)
					{
						Powermax[index]  = PoutTotcum;
						h[index] = height;
						w[index] = width;
						fl[index] = flowLength;						
						mdt[index] = mdot;
						apd[index] = AllowablePD;
						Pres[index] = PressDropTot;
						pin[index] = PDin;
						pout[index] = PDout;
						Tgdrop[index] = Tg[Ncv]-Tg[0];
						Qgas[index] = QgTotcum;
						Qrad[index] = QradTotcum;
						Qh[index] = QhtemTotcum;
						Qc[index] = QcTotcum;						
						nTEC0[index] = (nTETotcum0/32);	// SKU TEMs
						nTEC1[index] = (nTETotcum1/127);	// BiTE TEMs
						tf[index] = tfin;
						fs[index] = finSpacing;
						nfin[index] = Nf;						
	  					//PoutTotcum = PoutTot_[iCV];
					}
				}

			} // Fin Spacing Loop
		} // Fin Thickness Loop
		
		
		//outff.close();
				
		// Output Optimum Values for a geometry Triplet
		outopt	<<"\n" << gi
				<<"\t" << h[index]
				<<"\t" << w[index]
				<<"\t" << fl[index]
				<<"\t" << mdt[index]
				<<"\t" << Powermax[index]
				<<"\t" << Tgdrop[index]
				<<"\t" << Qgas[index]
				<<"\t" << Qrad[index]
				<<"\t" << Qh[index]
				<<"\t" << Qc[index]
				<<"\t" << Pres[index]
				<<"\t" << pin[index]
				<<"\t" << pout[index]
				<<"\t" << apd[index]
				<<"\t" << nTEC0[index]
				<<"\t" << nTEC1[index]
				<<"\t" << tf[index]
				<<"\t" << fs[index]
				<<"\t" << nfin[index];	

		 
		 //cout<<"\n "<<h[index]<<"\t"<<w[index]<<"\t"<<fl[index];
		 index++;
		
	}	// geom - Loop ends here

	
	
	outopt.close();	
	
	}
	// End of mdot loop
	
	cout <<"\n Program Ended " ;
	
	cin.get();
	return 0;	
}
