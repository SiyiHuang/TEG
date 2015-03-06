

	// For Global Variables
	TEM.numDensity[0] = 0;
	TEM.numDensity[1] = 0;
	TEM.etaNP[0] = 0;
	TEM.etaNP[1] = 0;
	TEM.nC[0] = 0;
	TEM.nC[1] = 0;
	TEM.length[0] = 0;
	TEM.length[1] = 0;
	TEM.width[0] = 0;
	TEM.width[1] = 0;
	TEM.height[0] = 0;
	TEM.height[1] = 0;
	TEM.epsilon = 0;

	TEC.Kl = 0;
	TEC.R = 0;
	TEC.Alpha = 0;
	TEC.Area = 0;
	TEC.height = 0;

	for(iCV = 0; iCV < Nx; iCV++)
	{
		for(j= 0; j <= 15; j++)	R[j][iCV] = 0.001;
		I[iCV] = SMALL;
		Ri[iCV] = SMALL;
		Rl[iCV] = SMALL;
		V[iCV] = SMALL;
		P[iCV] = SMALL;
		nTE[iCV] = 1;
		for(j = 0; j < 10; j++)	Q[j][iCV] = 0.001;
		for(j = 0; j < 10; j++)	T[j][iCV] = 0.001;
		nTETot_[iCV][0] = 0;
		nTETot_[iCV][1] = 0;
	}