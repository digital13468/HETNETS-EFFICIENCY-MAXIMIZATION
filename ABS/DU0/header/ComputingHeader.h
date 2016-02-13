#ifndef _COMPUTING_HEADER_H
#define _COMPUTING_HEADER_H 

double * SCoverageFinder(int *dBM, double *D, int Tbs, int Ntp, int Modulation){
	int i,j,k;
	double *coverage=malloc(Tbs*sizeof(double));
	SArrayInitialization(coverage,Tbs);
	for(i=0;i<Tbs;i++)
		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
				if(*(dBM+i*Ntp*Modulation+j*Modulation+k)==1)
					if(*(D+i*Ntp+j)/20>coverage[i])
						coverage[i]=*(D+i*Ntp+j)/20;
					//printf("%g\n",coverage[i]);
	return coverage;
}

Stemp_dBM_Initiator(int *ptr, int Ntp){	//initiate the array of temp_dBM for the new heuristic algorithm
	int i,j;
	for(i=0;i<2;i++)
		for(j=0;j<Ntp;j++)
			*(ptr+i*Ntp+j)=-1;
}

SArrayInitialization(double *ptr, int Tbs){
       	int i;         
       	for(i=0;i<Tbs;i++)
       		*(ptr+i)=0.0;
}

SArrayInitialization1(int *ptr, int Tbs){
       	int i;         
       	for(i=0;i<Tbs;i++)
       		*(ptr+i)=0;
}     

SArrayInitialization2(int *ptr, int Tbs, int Ntp, int Modulation){
       	int i,j,k;         
       	for(i=0;i<Tbs;i++)
       		for(j=0;j<Ntp;j++)
			for(k=0;k<Modulation;k++)
       				*(ptr+i*Ntp*Modulation+j*Modulation+k)=0;
}           

#endif
