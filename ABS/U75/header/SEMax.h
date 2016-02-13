 
//#include "ComparisonModel.h"
SEMax(int Tbs, int* DS, int* minM, int nUser, int Modulation, int Ntp, int* ptrAdjacentBSs, double mSBP, double mABP, double* P, double ABP, double SBP, int mNbs, double BP, double MP, double mP, int DSt, double* BR, int* ConnectionAssignment){
	FILE *PrintToFile;
	if((PrintToFile=fopen("SEMax.txt","w"))==NULL)
		printf("\nerror!Fail to open SEMax.txt!");
	else
		printf("\nOpen SEMax.txt successfully!\n");
	fprintf(PrintToFile,"This is the input to CPLEX for Spectrum Efficiency Optimization Model.\n");
	
	int count,i,k,j;
 
//	printf("maximize\n");
	//fprintf(PowerSavingCPLEX,"minimize\n");
	fprintf(PrintToFile,"maximize\n");
// 	for(count=0;count<Ntp;count++)
// 		for(i=0;i<Tbs;i++)
// 			for(k=0;k<Modulation;k++)
// 				if(count==Ntp-1&&i==Tbs-1&&k==Modulation-1)
// 					printf("%lf u%d_%d_%d",BR[count]/DS[count][k],i,count,k);
// 				else 
// 					printf("%lf u%d_%d_%d + ",BR[count]/DS[count][k],i,count,k);
	for(count=0;count<Ntp;count++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(count==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(PrintToFile,"%lf u%d_%d_%d",*(BR+count),i,count,k);
				else 
					fprintf(PrintToFile,"%lf u%d_%d_%d + ",*(BR+count),i,count,k);
	fprintf(PrintToFile,"\nst\n");
 
	cconstraint1(Tbs,DS,PrintToFile,minM,mNbs,Modulation,DS,DSt,Ntp);
	cconstraint2(Tbs,PrintToFile, mNbs, Ntp, Modulation, P, MP, mP);
	cconstraint3(Tbs,PrintToFile,nUser, Modulation);
	cconstraint4(Tbs,PrintToFile, Ntp, Modulation, BP);
	cconstraint5(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint10(Tbs,PrintToFile, Ntp, Modulation, ptrAdjacentBSs);
	cconstraint6(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint7(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint8(Tbs,PrintToFile, Ntp, Modulation);
	for(count=0;count<Ntp;count++)
		for(i=0;i<Tbs;i++)
			for(k=0;k<Modulation;k++)
				if(count==Ntp-1&&i==Tbs-1&&k==Modulation-1)
					fprintf(PrintToFile,"%d u%d_%d_%d",(*(DS+count*Modulation+k)),i,count,k);
				else 
					fprintf(PrintToFile,"%d u%d_%d_%d + ",(*(DS+count*Modulation+k)),i,count,k);
	fprintf(PrintToFile," = 1\n");
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
 
 
