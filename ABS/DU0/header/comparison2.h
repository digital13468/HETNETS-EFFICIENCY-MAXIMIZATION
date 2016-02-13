//#include "ComparisonModel.h"
comparison2(int Tbs, int* DS, int* minM, int nUser, int Modulation, int Ntp, int* ptrAdjacentBSs, double mSBP, double mABP, double* P, double ABP, double SBP, int mNbs, double BP, double MP, double mP, int DSt, double* BR){
	FILE *PrintToFile;
	if((PrintToFile=fopen("comparison2.txt","w"))==NULL)
		printf("\nerror!Fail to open comparison4.txt!");
	else
		printf("\nOpen comparison2.txt successfully!\n");
	fprintf(PrintToFile,"This is the input to CPLEX for comparison2.\n");
	
	cobjective(Tbs,PrintToFile,Ntp,Modulation,BR,DS);
	
//  	fprintf(NoExistingUser,"\nst\n");
// fprintf(PowerSavingCPLEX,"\nst\n");
//  printf("\nst\n");
 
	cconstraint1(Tbs,DS,PrintToFile,minM,mNbs,Modulation,DS,DSt,Ntp);
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
	int i, j, k;
	for (i=0; i<Tbs; i++)
		for (j=0; j<Ntp; j++)
			for (k=1; k<Modulation; k++)
				fprintf(PrintToFile, "u%d_%d_%d = 0\n", i, j, k);
// 	fprintf(NoExistingUser,"bounds\n");
// fprintf(PowerSavingCPLEX,"bounds\n");
//  printf("bounds\n");
	cbounds(Tbs,PrintToFile, Ntp, Modulation);
	cspecifyTypes(Tbs,PrintToFile, Ntp, Modulation);
	fclose(PrintToFile);
}  
