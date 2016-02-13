//#include "ComparisonModel.h"
SpecParti(int Tbs, int* DS, int* minM, int nUser, int Modulation, int Ntp, int* ptrAdjacentBSs, double mSBP, double mABP, double* P, double ABP, double SBP, int mNbs, double BP, double MP, double mP, int DSt, double* BR, int* ConnectionAssignment){
	FILE *PrintToFile;
	if((PrintToFile=fopen("SpectrumPartitionModel.txt","w"))==NULL)
		printf("\nerror!Fail to open SpectrumPartitionModel.txt!");
	else
		printf("\nOpen SpectrumPartitionModel.txt successfully!\n");
	fprintf(PrintToFile,"This is the input to CPLEX for Spectrum Partition Model.\n");
	
	cobjective(Tbs,PrintToFile,Ntp,Modulation,BR,DS);
	
	int i,j,k,l,m,n;
	int macroDSt = 34; 
	int smallDSt = 9;
	for(i=0; i<Tbs-mNbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(PrintToFile,"%d u%d_%d_%d", *(DS+j*Modulation+k),i, j,k);
					/*for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									fprintf(PrintToFile," + %d u%d_%d_%d",*(DS+m*Modulation+n),l,m,n);*/
				}		
				else 
					fprintf(PrintToFile,"%d u%d_%d_%d + ",*(DS+j*Modulation+k), i, j,k);
		}
	fprintf(PrintToFile," - %d t <= 0\n",macroDSt);	
	}
	for(i=Tbs-mNbs; i<Tbs; i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++)
				if(j==Ntp-1&&k==Modulation-1){
					fprintf(PrintToFile,"%d u%d_%d_%d", *(DS+j*Modulation+k),i, j,k);
					/*for(l=Tbs-mNbs;l<Tbs;l++)
						if(*(minM+i*Tbs+l)==1)
							for(m=0;m<Ntp;m++)
								for(n=0;n<Modulation;n++)
									fprintf(PrintToFile," + %d u%d_%d_%d",*(DS+m*Modulation+n),l,m,n);*/
				}		
				else 
					fprintf(PrintToFile,"%d u%d_%d_%d + ",*(DS+j*Modulation+k), i, j,k);
		}
	fprintf(PrintToFile," - %d t <= 0\n",smallDSt);	
	}
	cconstraint2(Tbs,PrintToFile, mNbs, Ntp, Modulation, P, MP, mP);
	cconstraint3(Tbs,PrintToFile,nUser, Modulation);
	cconstraint4(Tbs,PrintToFile, Ntp, Modulation, BP);
	cconstraint5(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint10(Tbs,PrintToFile, Ntp, Modulation, ptrAdjacentBSs);
	cconstraint6(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint7(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint8(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint9(Tbs,PrintToFile, mNbs, ABP, SBP, Ntp, Modulation, P, mABP, mSBP);
	cconstraint11(Tbs,PrintToFile,nUser, Ntp, Modulation);
// 	int i, j, k;
// 	for (i=0; i<Tbs; i++)
// 		for (j=0; j<Ntp; j++)
// 			if (*(ConnectionAssignment+j) == i)
// 				for (k=0; k<Modulation; k++)
// 					if (k == Modulation-1)
// 						fprintf (PrintToFile, "u%d_%d_%d - t = 0\n", i, j, k);
// 					else
// 						fprintf (PrintToFile, "u%d_%d_%d + ", i, j, k);
// 	fprintf(NoExistingUser,"bounds\n");
// fprintf(PowerSavingCPLEX,"bounds\n");
//  printf("bounds\n");
	cbounds(Tbs,PrintToFile, Ntp, Modulation);
	cspecifyTypes(Tbs,PrintToFile, Ntp, Modulation);
	fclose(PrintToFile);
}  
 
 
