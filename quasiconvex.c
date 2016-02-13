/* --------------------------------------------------------------------------

 * --------------------------------------------------------------------------
 */

/*  - Reading in and optimizing a problem.

 */

/* To run this example, no command line arguments are required.
   This program reads a problem from a file name "feasibility.lp" */

/* Bring in the CPLEX function declarations and the C library 
   header file stdio.h with the following single include. */

#include <ilcplex/cplex.h>

/* Bring in the declarations for the string, math and character functions 
   and malloc */
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "header/bisec.h"
#include <time.h>
//#define limit 0.429

void main (void)
{
	int Tbs, Ntp, DSt, mNbs, Modulation, nUser, tempI;
	double BP, MP, mP, SBP, mSBP, mABP, ABP, tempD;
	FILE *CplexIn, *parameterIn;
	int solstat;
	int Xbs, Ybs, Xms, Yms, ds, temp, FC;
	int BS, MS, Mod;
	double p,br;
	char id1[1000];
	int i, j, k, l, m, n;
	if((parameterIn=fopen("parameter.txt","r"))==NULL){
		printf("Fail to open parameter.txt!");
		exit(1);
	}
	else
		printf("Open parameter.txt successfully!\n");
	
	fscanf(parameterIn,"Tbs=%d\n",&tempI);
	//printf("%d",tempI);
	Tbs=tempI;
	fscanf(parameterIn,"Ntp=%d\n",&tempI);
	//printf("%d",tempI);
	Ntp=tempI;
	fscanf(parameterIn,"BP=%lf\n",&tempD);
	BP=tempD;
	fscanf(parameterIn,"DSt=%d\n",&tempI);
	DSt=tempI;
	fscanf(parameterIn,"mNbs=%d\n",&tempI);
	mNbs=tempI;
	fscanf(parameterIn,"Modulation=%d\n",&tempI);
	Modulation=tempI;
	fscanf(parameterIn,"mP=%lf\n",&tempD);
	mP=tempD;
	fscanf(parameterIn,"MP=%lf\n",&tempD);
	MP=tempD;
	fscanf(parameterIn,"nUser=%d\n",&tempI);
	nUser=tempI;
	fscanf(parameterIn,"SBP=%lf\n",&tempD);
	SBP=tempD;
	fscanf(parameterIn,"mSBP=%lf\n",&tempD);
	mSBP=tempD;
	fscanf(parameterIn,"mABP=%lf\n",&tempD);
	mABP=tempD;
	fscanf(parameterIn,"ABP=%lf\n",&tempD);
	ABP=tempD;
	double epsilon=0.1;
	double upper=1.2;
	double lower=0.2;
	double limit=0.0;
	int lastRound=0;

	int *DS=malloc(Ntp*Modulation*sizeof(int));
	double *P=malloc(Tbs*Ntp*Modulation*sizeof(double));
	double *BR=malloc(Ntp*sizeof(double));
	int *minM=malloc((Tbs-mNbs)*Tbs*sizeof(int));
	int *AdjacentBSs=malloc(Tbs*Ntp*sizeof(int));
	int iteration=0;
	i=0;
	while(i<(Tbs-mNbs)*Tbs){
		//printf("loop");
		fscanf(parameterIn,"minM[%d][%d]=%d\n",&BS,&FC,&temp);
		*(minM+BS*Tbs+FC)=temp;
		//printf("%d. minM[%d][%d]=%d\n",i,BS,FC,*(minM+BS*Tbs+FC));
		i++;
	}
	i=0;
	while(i<Ntp*Tbs){
		fscanf(parameterIn,"AdjacentBSs[%d][%d]=%d\n",&BS,&MS,&temp);
		*(AdjacentBSs+BS*Ntp+MS)=temp;
		//printf("%d. AdjacentBSs[%d][%d]=%d\n",i,BS,MS,*(AdjacentBSs+BS*Ntp+MS));
		i++;
	}

	fclose(parameterIn);
	
//	int dBM[Tbs][Ntp][Modulation];
//	int *ptrdBM;
//	int *dBA=malloc(Tbs*sizeof(int));

	
// 	double w[Tbs][Ntp][Modulation];
// 	double CumulativeP[Tbs]={0};
// 	double CumulativeBR[Tbs]={0};
// 	double v[Tbs]={0};
// 	int CumulativeDS[Tbs]={0};
// 	int CumSerMS[Tbs]={0};
// 	double TotalP=0.0;
// 	double TotalBR=0.0;
// 	double t=0.0;
// 	double NOV=0.0;
// 	int TotalDS=0;
// 	int TotSerMS=0;
	
	if ((CplexIn=fopen("outfile1.txt", "r")) == NULL){
		printf("Fail to open file!");
		exit(1);
	}
	else
		printf("Open outfile1.txt successfully!\n");
	//if ((CplexOut=fopen("SpectrumEfficiencyCPLEX_Output.txt", "w")) == NULL)
	//printf("Fail to open file!");
	//else
	//printf("Open file successfully!\n");
	fscanf(CplexIn, "%s", &id1);
	while (fscanf(CplexIn, "%s", &id1)!=EOF){
		//fprintf(CplexOut, "%s\n",id1);
		fscanf(CplexIn, " BS[%d]=%d,%d", &BS, &Xbs, &Ybs);
		//Xb[BS]=Xbs;
		//Yb[BS]=Ybs;
		//fprintf(CplexOut,"BS[%d]=(%d,%d) ", BS, Xbs, Ybs);
		fscanf(CplexIn, " MS[%d]=%d,%d", &MS, &Xms, &Yms);
		//Xtp[MS]=Xms;
		//Ytp[MS]=Yms;
		//fprintf(CplexOut,"MS[%d]=(%d,%d) ", MS, Xms, Yms);
		fscanf(CplexIn, " with modulation %d", &Mod);
		fscanf(CplexIn, " DS=%d", &ds);
		*(DS+MS*Modulation+Mod)=ds;
		//printf("DS:%d ", ds);
		fscanf(CplexIn, " power=%lf mW", &p);
		*(P+BS*Ntp*Modulation+MS*Modulation+Mod)=p;
		fscanf(CplexIn," BW=%lf",&br);
		*(BR+MS)=br;
		//printf("power:%g\n", p);
		fscanf(CplexIn," %s",&id1);
	}
	clock_t start_time, end_time;
	double timeTotal=0.0;
	start_time=clock();
	do{
		limit=(upper+lower)/2;
		printf("upper=%lf lower=%lf limit=%lf\n",upper,lower,limit);
		FILE *feasibility;
		if((feasibility=fopen("feasibility.lp","w"))==NULL){
			printf("\nerror! Fail to oopen the file of feasibility.lp!");
			exit(1);
		}
		else{
			printf("\nOPen feasibility.lp successfully!\n");
			printf("Iteration %d:\n",++iteration);
		}
		fprintf(feasibility,"'\'This is the input to CPLEX for feasibility check.(iteration %d)\nupper=%lf lower=%lf limit=%lf\n",iteration,upper,lower,limit);
		//objective function
		fprintf(feasibility,"maximize\n0\nst\n");
		//constraint 1
		for(i=0; i<Tbs-mNbs; i++){
			for(j=0;j<Ntp;j++){
				for(k=0;k<Modulation;k++)
					if(j==Ntp-1&&k==Modulation-1){
						fprintf(feasibility,"%d dBM%d_%d_%d", *(DS+j*Modulation+k),i, j,k);
						for(l=Tbs-mNbs;l<Tbs;l++)
							if(*(minM+i*Tbs+l)==1)
								for(m=0;m<Ntp;m++)
									for(n=0;n<Modulation;n++)
										fprintf(feasibility," + %d dBM%d_%d_%d",*(DS+m*Modulation+n),l,m,n);
					}		
					else 
						fprintf(feasibility,"%d dBM%d_%d_%d + ",*(DS+j*Modulation+k), i, j,k);
			}
		fprintf(feasibility," <= %d\n",DSt);	
		}
		//constraint 2
		for(i=0; i<Tbs-mNbs; i++){
			for(j=0;j<Ntp;j++){
				for(k=0;k<Modulation;k++)
					if(j==Ntp-1&&k==Modulation-1)
						fprintf(feasibility,"%lf dBM%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i, j,k);
					else 
						fprintf(feasibility,"%lf dBM%d_%d_%d + ",*(P+i*Ntp*Modulation+j*Modulation+k), i, j,k);
			}
			fprintf(feasibility," <= %lf\n",MP);
		}
		for(i=Tbs-mNbs; i<Tbs; i++){
			for(j=0;j<Ntp;j++){
				for(k=0;k<Modulation;k++)
					if(j==Ntp-1&&k==Modulation-1)
						fprintf(feasibility,"%lf dBM%d_%d_%d", *(P+i*Ntp*Modulation+j*Modulation+k),i, j,k);
					else 
						fprintf(feasibility,"%lf dBM%d_%d_%d + ",*(P+i*Ntp*Modulation+j*Modulation+k), i, j,k);
			}
			fprintf(feasibility," <= %lf\n",mP);
		}		
		//constraint 3
		for(i=0; i<nUser; i++)
			for(j=0; j<Tbs; j++) 
				for(k=0;k<Modulation;k++)
					if(j==Tbs-1&&k==Modulation-1)
						fprintf(feasibility,"dBM%d_%d_%d <= 1\n",j,i,k);
					else
						fprintf(feasibility,"dBM%d_%d_%d +",j,i,k);
		
		//constraint 4
		for(i=0; i<Tbs; i++) 
			for(j=0; j<Ntp; j++)
				for(k=0; k<Modulation; k++)
					if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
						fprintf(feasibility,"dBM%d_%d_%d",i,j,k);
					else 
						fprintf(feasibility,"dBM%d_%d_%d +",i,j,k);
		fprintf(feasibility," >= %.2lf\n",Ntp*(1-BP));		
		//constraint 5
		for(i=0;i<Tbs;i++)
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					if(j==Ntp-1&&k==Modulation-1)
						fprintf(feasibility,"dBM%d_%d_%d - dBA%d >= 0\n",i,j,k,i);
					else
						fprintf(feasibility,"dBM%d_%d_%d + ",i,j,k);
		//constraint 8
		for(i=0;i<Tbs;i++)
			for(j=0;j<Ntp;j++){
				for(k=0;k<Modulation;k++){
					if(k==Modulation-1)
						fprintf(feasibility,"dBM%d_%d_%d ",i,j,k);
					else
						fprintf(feasibility,"dBM%d_%d_%d + ",i,j,k);
				}
				fprintf(feasibility,"<= %d\n",*(AdjacentBSs+i*Ntp+j));	
			}
		//constraint 6
		for(i=0;i<Tbs;i++)
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					fprintf(feasibility,"dBM%d_%d_%d - dBA%d <= 0\n",i,j,k,i);	
		//constraint 7
		for(j=0;j<Ntp;j++)
			for(i=0;i<Tbs;i++)
				for(k=0;k<Modulation;k++)
					if(j==Ntp-1&&i==Tbs-1&&k==Modulation-1)
						fprintf(feasibility,"%lf dBM%d_%d_%d",*(BR+j)/(*(DS+j*Modulation+k)),i,j,k);
					else 
						fprintf(feasibility,"%lf dBM%d_%d_%d + ",*(BR+j)/(*(DS+j*Modulation+k)),i,j,k);
		fprintf(feasibility," - ");
		for(i=0;i<Tbs-mNbs;i++){
			fprintf(feasibility,"%lf dBA%d - ",(ABP-SBP)*limit,i);
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					if(i==Tbs-mNbs-1&&j==Ntp-1&&k==Modulation-1)
						fprintf(feasibility,"%lf dBM%d_%d_%d",*(P+i*Ntp*Modulation+j*Modulation+k)*limit,i,j,k);
					else
						fprintf(feasibility,"%lf dBM%d_%d_%d - ",*(P+i*Ntp*Modulation+j*Modulation+k)*limit,i,j,k);
		}
		for(i=Tbs-mNbs;i<Tbs;i++){
			fprintf(feasibility," - %lf dBA%d",(mABP-mSBP)*limit,i);
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					if(i==Tbs-1&&j==Ntp-1&&k==Modulation-1)
						fprintf(feasibility," - %lf dBM%d_%d_%d",*(P+i*Ntp*Modulation+j*Modulation+k)*limit,i,j,k);
					else
						fprintf(feasibility," - %lf dBM%d_%d_%d",*(P+i*Ntp*Modulation+j*Modulation+k)*limit,i,j,k);
				
		}
		fprintf(feasibility," - %lf dBS >= 0\n",(((Tbs-mNbs)*SBP)+(mNbs*mSBP))*limit);		

		//constraint 9
		for(j=nUser;j<Ntp;j++)
			for(i=0;i<Tbs;i++)
				for(k=0;k<Modulation;k++)
					if(i==Tbs-1&&k==Modulation-1)
						fprintf(feasibility,"dBM%d_%d_%d = 1\n",i,j,k);
					else
						fprintf(feasibility,"dBM%d_%d_%d + ",i,j,k);		
		//bounds
		fprintf(feasibility,"bounds\ndBS=1\n");
		//binary
		fprintf(feasibility, "binaries\n");
		for(i=0;i<Tbs;i++)
			fprintf(feasibility,"dBA%d\n",i);
		for(i=0;i<Tbs;i++)
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					fprintf(feasibility,"dBM%d_%d_%d\n",i,j,k);	
		
		//end
		fprintf(feasibility,"end\n");
		
	//}while();

	//return;
		fclose(feasibility);
		solstat=Bisection();
		printf("%d\n",solstat);
		if(solstat==103)
			upper=limit;
		else
			lower=limit;
		if((100*((upper-lower)/upper))<=epsilon&&lastRound==0)
			lastRound=1;
		else if(lastRound==1&&solstat==103)
			lastRound=1;
		else
			lastRound=0;
	}while(((100*((upper-lower)/upper))>epsilon)||lastRound==1);
	end_time=clock();
	if(solstat!=101){
		printf("No solutions found. Program Stops!\n");
		exit(1);
	}
	FILE *biOptR;
	if ((biOptR=fopen("bisectionOptimality.txt", "r")) == NULL){
		printf("Fail to open file!");
		exit(1);
	}
	else{
		//printf("text line 7");
		printf("Open bisectionOptimality.txt successfully!\n");
		//fflush(stdout);
		//printf("text line 6\n");
	}
	char tempC;
	//int value;
	int *dBM=malloc(Tbs*Ntp*Modulation*sizeof(int));
	int reachEnd=0;
	int *dBA=malloc(Tbs*sizeof(int));
	//fprintf(stderr, "error string");
	//printf("text line 5\n");
	while(!reachEnd){
		//printf("text line 0\n");
		fscanf(biOptR,"dB%c",&tempC);
		//printf("%c",tempC);
		if(tempC=='M'){
			fscanf(biOptR,"%d_%d_%d:",&BS,&MS,&Mod);
			//printf("%d %d %d ",BS,MS,Mod);
			fscanf(biOptR,"%s",&id1);
			fscanf(biOptR,"%s",&id1);
			//printf("%s",id1);
			fscanf(biOptR,"%d\n",&tempI);
			//printf("%d\n",tempI);
			//fprintf(stderr, "error string");
			*(dBM+BS*Ntp*Modulation+MS*Modulation+Mod)=tempI;
		}
		else if(tempC=='A'){
			fscanf(biOptR,"%d:",&BS);
			//printf("%d",BS);
			fscanf(biOptR,"%s",&id1);
			fscanf(biOptR,"%s",&id1);
			//printf("%s",id1);
			fscanf(biOptR,"%d\n",&tempI);
			//printf("%d\n",tempI);
			*(dBA+BS)=tempI;
		}
		else if(tempC=='S'){
			tempI=0;
			while(!feof(biOptR)){
				fscanf(biOptR,"%s",&id1);
				//printf("\n%s",id1);
				tempI++;
			}
			//printf("%d\n",tempI);
			if(tempI>5){
				printf("\ndBS is not at the end! Program Stops!\n");
// 				printf("%s\n",id1);
				exit(1);
			}
// 			fscanf(biOptR,"%s",&id1);
// 			fscanf(biOptR,"%s",&id1);
// 			fscanf(biOptR,"%s",&id1);
// 			printf("%s",id1);
// 			fscanf(biOptR,"%d",&tempI);
// 			printf("%d\n",tempI);
// 			//if(fscanf(biOptR,"%s",&id1)!=EOF){
// 				fscanf(biOptR,"%s",&id1);
// 				printf("dBS is not at the end! Program Stops!");
// 				printf("%s\n",id1);
// 			//	exit(1);
// 			//}
 			else
				reachEnd=1;
				//printf("dBS is not at the end! Program Stops!");
// 				printf("%s\n",id1);
				//exit(1);
			//}
		}
	}
	fclose(biOptR);
	
	double *CumulativeP=malloc(Tbs*sizeof(double));
	int *CumulativeDS=malloc(Tbs*sizeof(int));
	double *CumulativeBR=malloc(Tbs*sizeof(double));
	int *CumSerMS=malloc(Tbs*sizeof(int));
	double TotalP=0.0;
	int TotalDS=0;
	double TotalBR=0.0;
	int TotSerMS=0;
	double NOV=0.0;
	for(i=0;i<Tbs;i++){
		*(CumulativeP+i)=0.0;
		*(CumulativeDS+i)=0;
		*(CumulativeBR+i)=0.0;
		*(CumSerMS+i)=0;
	}
	for(i=0;i<Tbs;i++){
		for(j=0;j<Ntp;j++){
			for(k=0;k<Modulation;k++){
				if(*(dBM+i*Ntp*Modulation+j*Modulation+k)==1){
					//printf("i=%d j=%d k=%d\n",i,j,k);
					//printf("1=%lf\n",*(CumulativeP+i));
					*(CumulativeP+i)=*(CumulativeP+i)+*(P+i*Ntp*Modulation+j*Modulation+k);
					//printf("1=%lf\n",*(CumulativeP+i));
					TotalP=TotalP+*(P+i*Ntp*Modulation+j*Modulation+k);
					//printf("2=%lf\n",TotalP);
					*(CumulativeDS+i)=*(CumulativeDS+i)+*(DS+j*Modulation+k);
					//printf("3=%d\n",*(CumulativeDS+i));
					*(CumulativeBR+i)=*(CumulativeBR+i)+*(BR+j);
					//printf("4=%lf\n",*(CumulativeBR+i));
					TotalDS=TotalDS+*(DS+j*Modulation+k);
					//printf("5=%d\n",TotalDS);
					TotalBR=TotalBR+*(BR+j);
					//printf("6=%lf\n",TotalBR);
					*(CumSerMS+i)=*(CumSerMS+i)+1;
					//printf("7=%f\n",*(CumSerMS+i));
					TotSerMS=TotSerMS+1;
					//printf("8=%d\n",TotSerMS);
					NOV=NOV+(*(BR+j)/(*(DS+j*Modulation+k)));
					//printf("9=%lf\n",NOV);
					if(*(dBA+i)!=1){
						printf("BS%d should operate in active mode! Program stops.",i);
						exit(1);
					}
				}
			}
		}
	}
	for(i=0;i<Tbs;i++){
		if(*(CumulativeP+i)>0&&i<Tbs-mNbs){
			*(CumulativeP+i)=*(CumulativeP+i)+ABP;
			TotalP=TotalP+ABP;
		}
		else if(*(CumulativeP+i)>0&&i>=Tbs-mNbs){
			*(CumulativeP+i)=*(CumulativeP+i)+mABP;
			TotalP=TotalP+mABP;
		}
		else if(*(CumulativeP+i)==0&&i<Tbs-mNbs){
			*(CumulativeP+i)=*(CumulativeP+i)+SBP;
			TotalP=TotalP+SBP;
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					if(*(dBM+i*Ntp*Modulation+j*Modulation+k)==1){
						printf("BS%d should run in sleep mode! Program stops.",i);
						exit(1);
					}
		}
		else if(*(CumulativeP+i)==0&&i>=Tbs-mNbs){
			*(CumulativeP+i)=*(CumulativeP+i)+mSBP;
			TotalP=TotalP+mSBP;
			for(j=0;j<Ntp;j++)
				for(k=0;k<Modulation;k++)
					if(*(dBM+i*Ntp*Modulation+j*Modulation+k)==1){
						printf("BS%d should run in sleep mode! Program stops.",i);
						exit(1);
					}
		}
	}
	FILE *biOptW;
	if ((biOptW=fopen("finalBi.txt", "w")) == NULL){
		printf("Fail to open file!");
		exit(1);
	}
	else
		printf("Open finalBi,txt successfully!\n");

	for(i=0;i<Tbs;i++){
		if(i==0)
			fprintf(biOptW,"Macro BS Statistics\n");
		if(i==Tbs-mNbs)
			fprintf(biOptW,"Micro BS Statistics\n");
	fprintf(biOptW,"Power from BS %d: %lf\n",i,*(CumulativeP+i));
	fprintf(biOptW,"DS from BS %d: %d\n",i,*(CumulativeDS+i));
		fprintf(biOptW,"Bandwidth from BS %d: %lf\n",i,*(CumulativeBR+i));
	fprintf(biOptW,"MSs served by BS %d: %d\n\n",i,*(CumSerMS+i)); 
	}
	fprintf(biOptW,"Total Power %lf\n",TotalP);
	fprintf(biOptW,"Total DSs: %d\n",TotalDS);
	fprintf(biOptW,"Total served MSs: %d\n",TotSerMS);
	fprintf(biOptW,"Total Bandwidth: %lf\n",TotalBR);
	fprintf(biOptW,"Objective Value: %lf\n",TotalBR/TotalDS/TotalP);
	fprintf(biOptW,"New Objective Value: %lf\n",NOV/TotalP);
	timeTotal=(double)(end_time-start_time)/CLOCKS_PER_SEC;
	fprintf(biOptW,"Time spent in solving the problem: %lf\n",timeTotal);
	free(AdjacentBSs);
	free(minM);
	free(DS);
	free(P);
	free(BR);
	free(CumulativeP);
	free(CumulativeDS);
	free(CumulativeBR);
	free(CumSerMS);
	free(dBM);
	free(dBA);
	fclose(biOptW);
}  /* END main */







 
