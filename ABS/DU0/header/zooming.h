//#include "ComparisonModel.h"
zooming(int Tbs, int* DS, int* minM, int nUser, int Modulation, int Ntp, double mSBP, double mABP, double* P, double ABP, double SBP, int mNbs, double BP, double MP, double mP, int DSt, double* BR, int space, int* Xbs, int* Ybs, int* Xtp, int* Ytp, int Rbs, int rbs){
	FILE *PrintToFile;
	if((PrintToFile=fopen("zooming.txt","w"))==NULL)
		printf("\nerror!Fail to open zooming.txt!");
	else
		printf("\nOpen zooming.txt successfully!\n");
	fprintf(PrintToFile,"This is the input to CPLEX for no zooming out mechanism..\n");
	
	cobjective(Tbs,PrintToFile,Ntp,Modulation,BR,DS);
	
  int i,j;
  int* ptrAdjBSs;
  int AdjBSs[Tbs][Ntp];
  ptrAdjBSs=&AdjBSs[0][0];
  for(i=0;i<Tbs;i++)
	for(j=0;j<Ntp;j++)
		*(ptrAdjBSs+i*Ntp+j)=0;
  for(i=0;i<Tbs-mNbs;i++)
	for(j=0;j<Ntp;j++){
		
		//if((space*(distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j])))<=Rbs){
		if(space*sqrt( pow((Xbs[i]-Xtp[j]),2) + pow((Ybs[i]-Ytp[j]),2) )<=Rbs){
			//printf("Distance between user[%d] and BS[%d] is %lf\n",j,i,space*sqrt( pow((Xbs[i]-Xtp[j]),2) + pow((Ybs[i]-Ytp[j]),2) ));
			*(ptrAdjBSs+i*Ntp+j)=1;
		}
	}
  for(i=Tbs-mNbs;i<Tbs;i++)
	for(j=0;j<Ntp;j++){
		if(space*sqrt( pow((Xbs[i]-Xtp[j]),2) + pow((Ybs[i]-Ytp[j]),2) )<=rbs){
		//if((space*(distance(Xbs[i],Ybs[i],Xtp[j],Ytp[j])))<=rbs){
			//printf("Distance between user[%d] and BS[%d] is %lf\n",j,i,space*sqrt( pow((Xbs[i]-Xtp[j]),2) + pow((Ybs[i]-Ytp[j]),2) ));
			*(ptrAdjBSs+i*Ntp+j)=1;
		}
	}
//   for(i=0;i<Tbs;i++)
// 	for(j=0;j<Ntp;j++)
// 		if(*(ptrAdjBSs+i*Ntp+j)==1)
// 			printf("User[%d] is covered by BS[%d]\n",j,i);
  cconstraint10(Tbs,PrintToFile, Ntp, Modulation, ptrAdjBSs);
	cconstraint1(Tbs,DS,PrintToFile,minM,mNbs,Modulation,DS,DSt,Ntp);
	cconstraint2(Tbs,PrintToFile, mNbs, Ntp, Modulation, P, MP, mP);
	cconstraint3(Tbs,PrintToFile,nUser, Modulation);
	cconstraint4(Tbs,PrintToFile, Ntp, Modulation, BP);
	cconstraint5(Tbs,PrintToFile, Ntp, Modulation);
	
	cconstraint6(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint7(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint8(Tbs,PrintToFile, Ntp, Modulation);
	cconstraint9(Tbs,PrintToFile, mNbs, ABP, SBP, Ntp, Modulation, P, mABP, mSBP);
	cconstraint11(Tbs,PrintToFile,nUser, Ntp, Modulation);
// 	int i, j, k;
// 	for (i=0; i<Tbs; i++)
// 		for (j=0; j<Ntp; j++)
// 			for (k=1; k<Modulation; k++)
// 				fprintf(PrintToFile, "u%d_%d_%d = 0\n", i, j, k);
// 	fprintf(NoExistingUser,"bounds\n");
// fprintf(PowerSavingCPLEX,"bounds\n");
//  printf("bounds\n");
	cbounds(Tbs,PrintToFile, Ntp, Modulation);
	cspecifyTypes(Tbs,PrintToFile, Ntp, Modulation);
	fclose(PrintToFile);
}  
 
