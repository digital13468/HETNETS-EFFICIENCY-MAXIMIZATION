//#include "ComparisonModel.h"
void OptFixedABSpattern(int Tbs, int* DS, int* minM, int nUser, int Modulation, int Ntp, int* ptrAdjacentBSs, double mSBP, double mABP, double* P, double ABP, double SBP, int mNbs, double BP, double MP, double mP, int DSt, double* BR, int* ConnectionAssignment, double ABS, double CSB){
  
	char name[80];
	char fileName[80];

	char str1[10];
	char str2[10];
	sprintf(str1,"%g", ABS);
	sprintf(str2,"%g", CSB);
	//printf("%lf", ABS);
	strcpy(name, "Optimal Fixed(");
	
	strcat(name, str1);
	strcat(name, ",");
	
	strcat(name, str2);
	strcat(name, ")");
	//char fileName = "Fixed(" + char(ABS) + "," + char(CSB) + ").txt";
	strcpy(fileName, "OptFixed(");
	strcat(fileName, str1);
	strcat(fileName, ",");
	strcat(fileName, str2);
	strcat(fileName, ").txt");
	
	FILE *PrintToFile;
	if((PrintToFile=fopen( fileName,"w"))==NULL)
		printf("\nerror!Fail to open %s!", fileName);
	else
		printf("\nOpen %s successfully!\n", fileName);
	fprintf(PrintToFile,"This is the input to CPLEX for %s.\n", name);
	
	cobjective(Tbs,PrintToFile,Ntp,Modulation,BR,DS);
	
	int i, j, k;
	for (i = 0; i < Tbs - mNbs; i++)
	  fprintf(PrintToFile, "x%d - %g t = 0\n", i, 1 - ABS);
	for (i = Tbs - mNbs; i < Tbs; i++)
	  fprintf(PrintToFile, "x%d - %g t = 0\n", i, ABS);
//  	fprintf(NoExistingUser,"\nst\n");
// fprintf(PowerSavingCPLEX,"\nst\n");
//  printf("\nst\n");
	ABSconstraint1(Tbs,DS,PrintToFile,minM,mNbs,Modulation,DS,DSt,Ntp);
	//cconstraint1(Tbs,DS,PrintToFile,minM,mNbs,Modulation,DS,DSt,Ntp);
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
 	for (i=0; i<Tbs; i++)
 		for (j=0; j<Ntp; j++)
 			if (*(ConnectionAssignment+j) == i)
 				for (k=0; k<Modulation; k++)
 					if (k == Modulation-1)
 						fprintf (PrintToFile, "u%d_%d_%d - t = 0\n", i, j, k);
 					else
 						fprintf (PrintToFile, "u%d_%d_%d + ", i, j, k);
// 	fprintf(NoExistingUser,"bounds\n");
// fprintf(PowerSavingCPLEX,"bounds\n");
//  printf("bounds\n");
	ABSbounds(Tbs,PrintToFile, Ntp, Modulation, mNbs);
	//cbounds(Tbs,PrintToFile, Ntp, Modulation);
	cspecifyTypes(Tbs,PrintToFile, Ntp, Modulation);
	fclose(PrintToFile);
}  
 
