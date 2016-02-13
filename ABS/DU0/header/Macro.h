//#include "NoMod.h"
int* Macro(int Tbs, int* DS, double* P, double* BR, int* Xbs, int* Ybs, int* Xtp, int* Ytp, int* Connection, int Ntp, double BP, int Modulation, int MP, int DSt, double mABP, double mSBP, int mNbs, double mp, int* minM, int Rbs, int rbs, int space, double* ptrD, double ABP, double SBP, int* MacroUserMod){
	int i,j,k;
	//double ShortestDistance=1000000000;
	//int ConnectionAssignment[Ntp]={-1};
	//int* ptrConnectionAssignment=&ConnectionAssignment[0];
// 	int mnbs=mNbs;
// 	double msbp=mSBP;
// 	double mabp=mABP;
// 	double sbp=SBP;
// 	double abp=ABP;
	
// 	for(j=0;j<Ntp;j++){
// 		for(i=0;i<Tbs;i++){
// 			if(D[i][j]<ShortestDistance){
// 				
// 				if(i<Tbs-mNbs&&D[i][j]<=Rbs){
// 					ShortestDistance=D[i][j];
// 					ConnectionAssignment[j]=i;
// 				}
// 				else if(i>=Tbs-mNbs&&D[i][j]<=rbs){
// 					ShortestDistance=D[i][j];
// 					ConnectionAssignment[j]=i;
// 				}
// 			}
// 		}
// 		if(ShortestDistance>Rbs){
// 			printf("%lf MS%d [%d][%d] locates outside the cells!\n Porgram Stops!!\n",ShortestDistance,j,Xtp[j],Ytp[j]);
// 			exit(1);
// 		}
// 		
// 		ShortestDistance=1000000000;
// 	}
// 	
// 	int MConnectionAssignment[Ntp]={-1};
// 	int* ptrMConnectionAssignment=&MConnectionAssignment[0];
// 	for(j=0;j<Ntp;j++){
// 		for(i=0;i<Tbs-mNbs;i++){
// 			if(D[i][j]<ShortestDistance){
// 				if(D[i][j]<=Rbs){
// 					ShortestDistance=D[i][j];
// 					MConnectionAssignment[j]=i;
// 				}
// 			}
// 		}
// 		if(ShortestDistance>Rbs){
// 			printf("%lf MS%d [%d][%d] locates outside the cells!\n Porgram Stops!!\n",ShortestDistance,j,Xtp[j],Ytp[j]);
// 			exit(1);
// 		}
// 		ShortestDistance=1000000000;
// 	}
// 	Macro(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space);
// 	free(ptrMConnectionAssignment);
// 	BPSK(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mabp, msbp, mnbs, mP, minM, Rbs, rbs, space);
// 	QAM64(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mNbs, mP, minM, Rbs, rbs, space);
// 	QAM16(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mNbs, mP, minM, Rbs, rbs, space);
// 	RANDOM(Tbs,DS,P,ptrBR,Xbs,Ybs,Xtp,Ytp,ptrConnectionAssignment, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mnbs, mP, minM, Rbs, rbs, space);
// 	free(ptrConnectionAssignment);
	double TotalBandwidth=0.0, TotalPower=0.0, NOV=0.0;
	int TotalDS=0;
	double *Power=malloc(Tbs*sizeof(double));
	int *DataSub=malloc(Tbs*sizeof(int));
	for(i=0;i<Tbs;i++){
		*(Power+i)=0;
		*(DataSub+i)=0;
	}
	FILE *Macro;
 	if ((Macro=fopen("Macro.txt", "w")) == NULL){
		printf("\n\nerror!Fail to open file!");
		exit(1);
	}
	else
		printf("\n\nOpen Macro.txt successfully!\n");

	fprintf(Macro,"\n#This file is meant to display the Macro system.\n");
	fprintf(Macro,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
	for(i=0;i<Tbs-mNbs;i++)
		fprintf(Macro," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],Rbs/space);
	fprintf(Macro,"\n\n");
	for(i=Tbs-mNbs;i<Tbs;i++)
		fprintf(Macro," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],rbs/space);
	fprintf(Macro,"\n\n");
	int* UserModulation=malloc(sizeof(int)*Ntp);
        for(i=0;i<Tbs-mNbs;i++){
		fprintf(Macro,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
        	for(j=0;j<Ntp;j++)
		//	for(k=0;k<Modulation;k++)
				if(*(Connection+j)<0){
					printf("MS has no coonections! Program Stops!!\n");
					exit(1);
				}
				else if(*(Connection+j)==i){
					if(*(ptrD+i*Ntp+j)<1250)
						k=3;//k=6;
					else if(*(ptrD+i*Ntp+j)<2500)
						k=2;//k=5;
					else if(*(ptrD+i*Ntp+j)<3750)
						k=1;//k=4;
					//else if(*(ptrD+i*Ntp+j)<2800)
					//	k=2;//k=3;
					//else if(*(ptrD+i*Ntp+j)<3500)
					//	k=1;//k=2;
					//else if(*(ptrD+i*Ntp+j)<4200)
					//	k=1;
					else
						k=0;
//for(i=0;i<Ntp;i++){
					
		//for(j=0;j<Tbs;j++)
			//if(ConnectionAssignment[i]==j){
					Power[i]+=*(P+i*Ntp*Modulation+j*Modulation+k);
					if(i<Tbs-mNbs&&Power[i]>MP){
						printf("BS%d has too many power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}
					else if(i>=Tbs-mNbs&&Power[i]>mp){
						printf("BS%d has too many power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}
					DataSub[i]+=*(DS+j*Modulation+k);
						if(DataSub[i]>DSt){
							printf("BS%d has too many subcarrier demands. Program Stops!\n",i);
							exit(1);
						}
					
					TotalBandwidth+=BR[j];
					TotalDS+=*(DS+j*Modulation+k);
					//DataSub[j]+=*(DS+i);
					
					TotalPower+=*(P+i*Ntp*Modulation+j*Modulation+k);
			//}		
					NOV+=(*(BR+j)/(*(DS+j*Modulation+k)));
					
					fprintf(Macro,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
					fprintf(Macro,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW BR %lf with modulation %d\n",i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],*(DS+j*Modulation+k),*(P+i*Ntp*Modulation+j*Modulation+k),BR[j],k);
					*(UserModulation+j)=k;
	//}
				}
// 		if(Power[i]>0){
// 			if(i<Tbs-mNbs){
				Power[i]+=ABP;
				TotalPower+=ABP;
// 			}
// 			else{
// 				Power[i]+=mABP;
// 				TotalPower+=mABP;
// 			}
// 		}
// 		else{	
// 			if(i<Tbs-mNbs){	
// 				Power[i]+=SBP;
// 				TotalPower+=SBP;
// 			}
// 			else{
// 				Power[i]+=mSBP;
// 				TotalPower+=mSBP;
// 			}
// 		}
	}
	j=0;
	while(j<Modulation){
		fprintf(Macro,"\n\n");
		for(i=0;i<Ntp;i++)
			if(*(UserModulation+i)==j)
				fprintf(Macro,"MS[%d] %d %d DS %d BR %lf with modulation %d\n",i,Xtp[i],Ytp[i],*(DS+i*Modulation+j),BR[i],j);
		j++;
	}
	fprintf(Macro,"\n\n");
	int temp;
	for(i=0;i<Tbs-mNbs;i++){
		temp=DataSub[i];
		for(j=Tbs-mNbs;j<Tbs;j++)
			if(*(minM+i*Tbs+j)==1)
				temp+=DataSub[j];
		if(temp>DSt){
			printf("The group of BS%d has too many subcarrier demands. Program Stops!\n",i);
			exit(1);
		}				
	}
	fprintf(Macro,"Total Power %lf\n",TotalPower);
        fprintf(Macro,"Total DSs: %d\n",TotalDS);
	fprintf(Macro,"Total BR: %lf\n",TotalBandwidth);
        //fprintf("Total served MSs: %d\n",TotalServedMSs);
	//fprintf(FixedCellSize,"Objective Value: %lf\n",TotalBandwidth/TotalDS/TotalPower);
	fprintf(Macro,"Spectrum-energy Efficiency: %lf\n",NOV/TotalPower);
	fprintf(Macro,"\n\n");
	//int MacroUserMod[Ntp];
	//int* ptrMacroUserMod = &MacroUserMod[0];
	for (i=0; i<Ntp; i++)
	  MacroUserMod[i] = *(UserModulation+i);
	free(Power);
	free(DataSub);
	free(UserModulation);	
	for(i=0;i<Ntp;i++)
	  printf("%d\t",*(MacroUserMod+i));
	printf("\n");
	fprintf(Macro," plot \"coordinates.txt\" index 0:0 using 2:3 title \"users\" with points, \"FixedCellSize.txt\" index 0:0 using 2:3:1 notitle with labels, \"FixedCellSize.txt\" index 0:0 using 2:3:4 title \"Macro Cell\" with circles, \"FixedCellSize.txt\" index 1:1 using 2:3:1 notitle with labels, \"FixedCellSize.txt\" index 1:1 using 2:3:4 title \"micro cell\" with circles");
	for(i=0;i<Tbs;i++)
		fprintf(Macro,", \"FixedCellSize.txt\" index %d:%d using 5:6 title \"BS[%d]\" with lines",i+2,i+2,i);
	for(i=0;i<Modulation;i++)
		fprintf(Macro,", \"FixedCellSize.txt\" index %d:%d using 2:3 title \"Modulation%d\" with points",Tbs+2+i,Tbs+2+i,i);
	fclose(Macro);
	//NoMod(Tbs, DS, P, BR, Xbs, Ybs, Xtp, Ytp, Connection, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mNbs, mp, minM, Rbs, rbs, space, ptrD, ABP, SBP);
	//return ptrMacroUserMod;
}
