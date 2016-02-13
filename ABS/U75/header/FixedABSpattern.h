//#include "NoMod.h"
void FixedABSpattern(double ABS, double CSB, int Tbs, int* DS, double* P, double* BR, int* Xbs, int* Ybs, int* Xtp, int* Ytp, int* preConnection, int Ntp, double BP, int Modulation, int MP, int DSt, double mABP, double mSBP, int mNbs, double mp, int* minM, int Rbs, int rbs, int space, double* ptrD, double ABP, double SBP, double maxCellSize, int* preUserMod, int nUser, int* ptrAdjacentBSs){
	int i,j, l;
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
  //int* UserMod = malloc(Ntp*sizeof(int));
  int *UserMod=malloc(Ntp*sizeof(int));
  int *Connection=malloc(Ntp*sizeof(int));
  
	double TotalBandwidth=0.0, TotalPower=0.0, NOV=0.0, TotalTransmitPower = 0.0;
	int TotalDS=0;
	double *Power=malloc(Tbs*sizeof(double));
	int *DataSub=malloc(Tbs*sizeof(int));
	int *userNum=malloc(Tbs*sizeof(int));
	double *transmitPower=malloc(Tbs*sizeof(double));
	for(i=0;i<Tbs;i++){
		*(Power+i)=0;
		*(DataSub+i)=0;
		*(userNum+i) = 0;
		*(transmitPower+i) = 0;
		
	}
	
	for (i = 0; i < Ntp; i++){
	  *(UserMod + i) = *(preUserMod + i);
	  *(Connection + i) = *(preConnection + i);
	}
	//char name[80] = "Fixed(" + char(ABS) + "," + char(CSB) + ")";
	char name[80];
	char fileName[80];

	char str1[10];
	char str2[10];
	sprintf(str1,"%g", ABS);
	sprintf(str2,"%g", CSB);
	//printf("%lf", ABS);
	strcpy(name, "Fixed(");
	
	strcat(name, str1);
	strcat(name, ",");
	
	strcat(name, str2);
	strcat(name, ")");
	//char fileName = "Fixed(" + char(ABS) + "," + char(CSB) + ").txt";
	strcpy(fileName, "Fixed(");
	strcat(fileName, str1);
	strcat(fileName, ",");
	strcat(fileName, str2);
	strcat(fileName, ").txt");

	
	FILE *Macro;
 	if ((Macro=fopen(fileName, "w")) == NULL){
		printf("\n\nerror!Fail to open file!");
		exit(1);
	}
	else
		printf("\n\nOpen %s successfully!\n", fileName);
  //printf("%s\n", name);
	fprintf(Macro,"\n#This file is meant to display the %s system.\n", name);
	fprintf(Macro,"#BS%dMS%dBP%g\n",Tbs,Ntp,BP);
	for(i=0;i<Tbs-mNbs;i++)
		fprintf(Macro," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],Rbs/space);
	fprintf(Macro,"\n\n");
	for(i=Tbs-mNbs;i<Tbs;i++)
		fprintf(Macro," \"BS[%d]\" %d %d %d\n",i,Xbs[i],Ybs[i],rbs/space);
	fprintf(Macro,"\n\n");
	//int* UserModulation=malloc(sizeof(int)*Ntp);
        for(i=0; i<Tbs; i++){
		//fprintf(Macro,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
        	for(j=0;j<Ntp;j++)
		//	for(k=0;k<Modulation;k++)
				if(*(Connection+j)<0){
					printf("MS has no coonections! Program Stops!!\n");
					exit(1);
				}
				else if(*(Connection+j)==i){/*
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
						k=0;*/
//for(i=0;i<Ntp;i++){
					
		//for(j=0;j<Tbs;j++)
			//if(ConnectionAssignment[i]==j){
					Power[i]+=*(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j));/*
					if(i<Tbs-mNbs&&Power[i]>MP){
					  //for (l = 0; l < Tbs-mNbs; l++)
					    //if (l != i && *(ptrD+l*Ntp+j) <= Rbs)
						printf("BS%d has too many power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}
					else if(i>=Tbs-mNbs&&Power[i]>mp){
						printf("BS%d has too many power demands (%lf)! Program Stops!\n",i,Power[i]);
						exit(1);
					}*/
					DataSub[i]+=*(DS+j*Modulation+*(UserMod+j));/*
						if(DataSub[i]>DSt){
							printf("BS%d has too many subcarrier demands. Program Stops!\n",i);
							exit(1);
						}*/
					
					//TotalBandwidth+=BR[j];
					//TotalDS+=*(DS+j*Modulation+k);
					//DataSub[j]+=*(DS+i);
					
					//TotalPower+=*(P+i*Ntp*Modulation+j*Modulation+k);
			//}		
					//NOV+=(*(BR+j)/(*(DS+j*Modulation+k)));
					
					//fprintf(Macro,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
					//fprintf(Macro,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW BR %lf with modulation %d\n",i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],*(DS+j*Modulation+k),*(P+i*Ntp*Modulation+j*Modulation+k),BR[j],k);
					//*(UserModulation+j)=k;
	//}
				}
// 		if(Power[i]>0){
// 			if(i<Tbs-mNbs){
				//Power[i]+=ABP;
				//TotalPower+=ABP;
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
	
	for(i=Tbs-1;i>-1;i--){
		//fprintf(Macro,"\n\nBS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
        	for(j=0;j<Ntp;j++)
		//	for(k=0;k<Modulation;k++)
				if(*(Connection+j)<0){
					printf("MS has no coonections! Program Stops!!\n");
					exit(1);
				}
				else if(*(Connection+j)==i){/*
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
						k=0;*/
//for(i=0;i<Ntp;i++){
					
		//for(j=0;j<Tbs;j++)
			//if(ConnectionAssignment[i]==j){
					//Power[i]+=*(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j));
					if(i<Tbs-mNbs){
					  for (l = Tbs-mNbs; l < Tbs; l++){
					    //if (*(ptrD+l*Ntp+j) <= rbs){
					      int n = -1;
					      if(*(ptrD+l*Ntp+j)< (rbs/Modulation))
						n=3;//k=6;
					      else if(*(ptrD+l*Ntp+j) < (2 * (rbs / Modulation)))
						n=2;//k=5;
					      else if(*(ptrD+l*Ntp+j) < (3 * (rbs / Modulation)))
						n=1;//k=4;
					      else
						n=0;
					      
					      if (*(ptrD+l*Ntp+j) <= maxCellSize && *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j)) >= *(P+l*Ntp*Modulation+j*Modulation+n) - CSB && Power[l] + *(P+l*Ntp*Modulation+j*Modulation+n) <= mp && DataSub[l] + *(DS+j*Modulation+n) <= DSt * ABS){
						fprintf(Macro, "#UE%d is associated to new eNodeB%d from %d(%lf %lf %lf).\n", j, l, i, *(P+l*Ntp*Modulation+j*Modulation+n), *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j)), CSB);
						Power[i] -= *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j));
						DataSub[i] -= *(DS+j*Modulation+*(UserMod+j));
						*(Connection+j) = l;
						*(UserMod+j) = n;
						Power[l] += *(P+l*Ntp*Modulation+j*Modulation+n);
						DataSub[l] += *(DS+j*Modulation+n);
						if (*(ptrD+l*Ntp+j) > maxCellSize){
						  fprintf(Macro, "something wrong (%lf %lf)\n", *(ptrD+l*Ntp+j), maxCellSize);
						  exit(1);
						}
					      }
					    //}
					    
					  }
					  if (Power[i] > MP || DataSub[i] > DSt * (1 - ABS)){
					  //s&&Power[i]>MP)
					    for (l = 0; l < Tbs-mNbs; l++){
					      if (l != i && *(ptrD+l*Ntp+j) <= Rbs){
						int n = -1;
						if(*(ptrD+l*Ntp+j) < (Rbs/Modulation))
						  n=3;//k=6;
						else if(*(ptrD+l*Ntp+j) < (2 * (Rbs / Modulation)))
						  n=2;//k=5;
						else if(*(ptrD+l*Ntp+j) < (3 * (Rbs / Modulation)))
						  n=1;//k=4;
						else
						  n=0;
					      
						if (Power[l] + *(P+l*Ntp*Modulation+j*Modulation+n) <= MP && DataSub[l] + *(DS+j*Modulation+n) <= DSt * (1 - ABS)){
						  fprintf(Macro, "#UE%d is associated to new eNodeB%d from %d(%lf %lf).\n", j, l, i, *(P+l*Ntp*Modulation+j*Modulation+n), *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j)));
						  Power[i] -= *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j));
						  DataSub[i] -= *(DS+j*Modulation+*(UserMod+j));
						  *(Connection+j) = l;
						  *(UserMod+j) = n;
						  Power[l] += *(P+l*Ntp*Modulation+j*Modulation+n);
						  DataSub[l] += *(DS+j*Modulation+n);
						  if (*(ptrD+l*Ntp+j) > Rbs){
						    printf("something wrong (%lf/%lf)\n", *(ptrD+l*Ntp+j), Rbs);
						    exit(1);
						  }
						}
					      }
					      
					      if (l == Tbs-mNbs-1 && Power[i] > MP || DataSub[i] > DSt * (1 - ABS)){
						printf("BS%d has too much power (%lf, %lf) or many RB demand (%d, %d)! Program Stops!\n",i,Power[i], MP, DataSub[i], DSt * (1 - ABS));
						exit(1);
					      }
					      else if (Power[i] <= MP && DataSub[i] <= DSt * (1 - ABS))
						l = Tbs;
					    }
					    
					  }
					}
					else {
					  if(Power[i]>mp || DataSub[i]>DSt * ABS){
					    for (l = 0; l < Tbs-mNbs; l++){
						if (*(ptrD+l*Ntp+j) <= Rbs){
						  int n = -1;
						  if(*(ptrD+l*Ntp+j) < (Rbs/Modulation))
						    n=3;//k=6;
						  else if(*(ptrD+l*Ntp+j) < (2 * (Rbs / Modulation)))
						    n=2;//k=5;
						  else if(*(ptrD+l*Ntp+j) < (3 * (Rbs / Modulation)))
						    n=1;//k=4;
						  else
						    n=0;
						
						  if (Power[l] + *(P+l*Ntp*Modulation+j*Modulation+n) <= MP && DataSub[l] + *(DS+j*Modulation+n) <= DSt * (1 - ABS)){
						    fprintf(Macro, "#UE%d is associated to new eNodeB%d from %d(%lf %lf).\n", j, l, i, *(P+l*Ntp*Modulation+j*Modulation+n), *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j)));
						    Power[i] -= *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j));
						    DataSub[i] -= *(DS+j*Modulation+*(UserMod+j));
						    *(Connection+j) = l;
						    *(UserMod+j) = n;
						    Power[l] += *(P+l*Ntp*Modulation+j*Modulation+n);
						    DataSub[l] += *(DS+j*Modulation+n);
						    if (*(ptrD+l*Ntp+j) > Rbs){
						      printf("something wrong (%lf/%lf)\n", *(ptrD+l*Ntp+j), Rbs);
						      exit(1);
						    }
						  }
						}
						if (l == Tbs-mNbs-1 && Power[i] > mp || DataSub[i] > DSt * ABS){
						  printf("BS%d has too much power (%lf, %lf) or many RB demand (%d, %d)! Program Stops!\n",i,Power[i], mp, DataSub[i], DSt * ABS);
						  exit(1);
						}
						else if (Power[i] <= mp && DataSub[i] <= DSt * ABS)
						  l = Tbs;
					      }
					  
						//printf("BS%d has too many power demands (%lf)! Program Stops!\n",i,Power[i]);
						//exit(1);
					  }
					//DataSub[i]+=*(DS+j*Modulation+k);
						/*if(DataSub[i]>DSt){
							printf("BS%d has too many subcarrier demands. Program Stops!\n",i);
							exit(1);
						}*/
					}
					//TotalBandwidth+=BR[j];
					//TotalDS+=*(DS+j*Modulation+k);
					//DataSub[j]+=*(DS+i);
					
					//TotalPower+=*(P+i*Ntp*Modulation+j*Modulation+k);
			//}		
					//NOV+=(*(BR+j)/(*(DS+j*Modulation+k)));
					
					//fprintf(Macro,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
					//fprintf(Macro,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW BR %lf with modulation %d\n",i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],*(DS+j*Modulation+k),*(P+i*Ntp*Modulation+j*Modulation+k),BR[j],k);
					//*(UserModulation+j)=k;
	//}
				}
// 		if(Power[i]>0){
// 			if(i<Tbs-mNbs){
				//Power[i]+=ABP;
				//TotalPower+=ABP;
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
	
	//int userNum[Tbs]={0};
	for (i = 0; i < Tbs; i++){
	  fprintf(Macro,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
	  for (j = 0; j < Ntp; j++){
	    if(*(Connection+j)<0){
		printf("MS has no coonections! Program Stops!!\n");
		exit(1);
	    }
	    else if(*(Connection+j)==i){
	      *(userNum+i) = *(userNum) + 1;
	      TotalBandwidth += *(BR + j);
	      TotalDS += *(DS+j*Modulation+*(UserMod+j));
	      TotalPower += *(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j));
	      NOV += *(BR+j)/(*(DS+j*Modulation+*(UserMod+j)));
	      fprintf(Macro,"BS[%d] %d %d BS[%d] %d %d\n",i,Xbs[i],Ybs[i],i,Xbs[i],Ybs[i]);
	      fprintf(Macro,"BS[%d] %d %d MS[%d] %d %d DS %d power %g mW BR %lf with modulation %d\n",i,Xbs[i],Ybs[i],j,Xtp[j],Ytp[j],*(DS+j*Modulation+*(UserMod+j)),*(P+i*Ntp*Modulation+j*Modulation+*(UserMod+j)),BR[j],*(UserMod+j));
	    }
	    
	  }
	  fprintf(Macro, "\n\n");
	  transmitPower[i] = Power[i];
	  TotalTransmitPower += Power[i];
	  if (*(userNum+i) > 0){
	    if (i < Tbs - mNbs){
	     // transmitPower[i] = Power[i];
	      Power[i] += ABP;
	      TotalPower += ABP;
	    }
	    else{
	      //transmitPower[i] = Power[i];
	      Power[i] += mABP;
	      TotalPower += mABP;
	    }
	  }
	  else{
	    if (i < Tbs - mNbs){
	     // transmitPower[i] = Power[i];
	      Power[i] += SBP;
	      TotalPower += SBP;
	    }
	    else{
	      //transmitPower[i] = Power[i];
	      Power[i] += mSBP;
	      TotalPower += mSBP;
	    }
	  }
	}
	
	j=0;
	while(j<Modulation){
		fprintf(Macro,"\n\n");
		for(i=0;i<Ntp;i++)
			if(*(UserMod+i)==j)
				fprintf(Macro,"MS[%d] %d %d DS %d BR %lf with modulation %d\n",i,Xtp[i],Ytp[i],*(DS+i*Modulation+j),BR[i],j);
		j++;
	}
	fprintf(Macro,"\n\n");/*
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
	}*/
	fprintf(Macro,"Total Power %lf\n",TotalPower);
        fprintf(Macro,"Total DSs: %d\n",TotalDS);
	fprintf(Macro,"Total BR: %lf\n",TotalBandwidth);
        //fprintf("Total served MSs: %d\n",TotalServedMSs);
	//fprintf(FixedCellSize,"Objective Value: %lf\n",TotalBandwidth/TotalDS/TotalPower);
	fprintf(Macro,"Spectrum-energy Efficiency: %lf\n",NOV/TotalPower);
	fprintf(Macro,"Transmit Spectrum-energy Efficiency: %lf\n",NOV/TotalTransmitPower);
	fprintf(Macro,"\n\n");
	//int MacroUserMod[Ntp];
	//int* ptrMacroUserMod = &MacroUserMod[0];
//	for (i=0; i<Ntp; i++)
//	  MacroUserMod[i] = *(UserModulation+i);
	free(Power);
	free(DataSub);
	//free(UserModulation);	
	free(transmitPower);
	free(userNum);
//	for(i=0;i<Ntp;i++)
//	  printf("%d\t",*(MacroUserMod+i));
	printf("\n");
	fprintf(Macro," plot \"coordinates.txt\" index 0:0 using 2:3 title \"users\" with points, \" %s \" index 0:0 using 2:3:1 notitle with labels, \" %s \" index 0:0 using 2:3:4 title \"Macro Cell\" with circles, \" %s \" index 1:1 using 2:3:1 notitle with labels, \" %s \" index 1:1 using 2:3:4 title \"micro cell\" with circles", fileName, fileName, fileName, fileName);
	for(i=0;i<Tbs;i++)
		fprintf(Macro,", \" %s \" index %d:%d using 5:6 title \"BS[%d]\" with lines", fileName, i+2,i+2,i);
	for(i=0;i<Modulation;i++)
		fprintf(Macro,", \" %s \" index %d:%d using 2:3 title \"Modulation%d\" with points",fileName, Tbs+2+i,Tbs+2+i,i);
	fclose(Macro);
	//NoMod(Tbs, DS, P, BR, Xbs, Ybs, Xtp, Ytp, Connection, Ntp, BP, Modulation, MP, DSt, mABP, mSBP, mNbs, mp, minM, Rbs, rbs, space, ptrD, ABP, SBP);
	//return ptrMacroUserMod;
	//free(UserMod);
	OptFixedABSpattern(Tbs, DS, minM, nUser, Modulation, Ntp, ptrAdjacentBSs, mSBP, mABP, P, ABP, SBP, mNbs, BP, MP, mp, DSt, BR, Connection, ABS, CSB);
	free (UserMod);
	free (Connection);
}
