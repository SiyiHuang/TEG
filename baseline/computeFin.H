
	pitch = finSpacing + tfin; // Fin Array Pitch
	Nf = int(width/pitch);	// Total number of Fins
	pitch = width/Nf;
	finSpacing = pitch - tfin;
	//if(finSpacing <= 0.001)
	//	break;

	//cout<<"\n "<<Nf<<"\t"<<tfin<<"\t"<<finSpacing<<endl;
	
	Af = Nf*(2.0*(height/2.0))*dx;	// Total Fin area in the CV
	Ab = (width - Nf*tfin)*dx;	// Remaining Unfinnned Base area in the CV
	At = Nf*Af + Ab; // Total area = Base Area + fin Area