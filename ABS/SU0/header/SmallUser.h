 
void SmallUser(int* Xtp, int* Ytp, int Px, int Py, int Tbs, int* Xbs, int* Ybs, int nUser, int mNbs, int space, int maxCellSize, int Rbs, int Ntp){
int i,j,k;

	double ShortestDistance;
	double tempDistance[Tbs];
	double smallUserRatio = 0.912;
	for(i=0;i<Tbs;i++)
		tempDistance[i]=0.0;
/*initialize random seed:*/
srand(time(NULL));
  
/*generate number:*/
for(i=0;i<nUser*smallUserRatio;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=Tbs-mNbs;j<Tbs;j++){
		tempDistance[j]=space*sqrt( pow((Xbs[j]-Xtp[i]),2) + pow((Ybs[j]-Ytp[i]),2) );
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>maxCellSize||ShortestDistance==0){
		
			i=i-1;
		
		printf("%d renew\n",i);
		//continue;
	}
	else{
	  ShortestDistance=1000000000;
	  for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*sqrt( pow((Xbs[j]-Xtp[i]),2) + pow((Ybs[j]-Ytp[i]),2) );
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	  }
	  if(ShortestDistance>Rbs||ShortestDistance==0){
		
			i=i-1;
		
		printf("%d renew\n",i);
		//continue;
	  }
	}
	  
	for(k=i-1;k>-1;k--)
		//if(i!=j)
		if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
			i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
                
}
for(i=nUser*smallUserRatio;i<nUser;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*sqrt( pow((Xbs[j]-Xtp[i]),2) + pow((Ybs[j]-Ytp[i]),2) );
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs||ShortestDistance==0){
		
			i=i-1;
		
		printf("%d renew\n",i);
		//continue;
	}
	//else if((Xtp[i]>X1urban&&Xtp[i]<X2urban)&&(Ytp[i]>Y1urban&&Ytp[i]<Y2urban)){
	//	i=i-1;
	//	printf("%d renew\n",i);
	//}
		
	for(k=i-1;k>-1;k--)
		//if(i!=j)
		if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
			i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
                
}/*
for(i=nUser*(UrbanUser+SubUser);i<nUser;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		
			i=i-1;
		
		printf("%d renew\n",i);
		//continue;
	}
	else if((Xtp[i]>X1suburban&&Xtp[i]<X2suburban)&&(Ytp[i]>Y1suburban&&Ytp[i]<Y2suburban)){
		i=i-1;
		printf("%d renew\n",i);
	}
	for(k=i-1;k>-1;k--)
		//if(i!=j)
		if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
			i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
                
}*/
for(i=nUser;i<nUser+ceil((Ntp-nUser)*smallUserRatio);i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	
	for(j=Tbs-mNbs;j<Tbs;j++){
		tempDistance[j]=space*sqrt( pow((Xbs[j]-Xtp[i]),2) + pow((Ybs[j]-Ytp[i]),2) );
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>maxCellSize||ShortestDistance==0){
		
			i=i-1;
		
		printf("%d renew\n",i);
		//continue;
	}
	else{
	  ShortestDistance=1000000000;
	  for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*sqrt( pow((Xbs[j]-Xtp[i]),2) + pow((Ybs[j]-Ytp[i]),2) );
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	  }
	  if(ShortestDistance>Rbs||ShortestDistance==0){
		
			i=i-1;
		
		printf("%d renew\n",i);
		//continue;
	  }
	}
	for(k=i-1;k>-1;k--)
		//if(i!=j)
		if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
			i=i-1;
				//printf("error!!MS[%d] and [%d] are at the same point.\n",i,j);
                		//exit(1);
                
}

for(i=nUser+ceil((Ntp-nUser)*smallUserRatio);i<Ntp;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*sqrt( pow((Xbs[j]-Xtp[i]),2) + pow((Ybs[j]-Ytp[i]),2) );
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs||ShortestDistance==0){
		i=i-1;
		printf("%d renew\n",i);
	}
	//else if((Xtp[i]>X1urban&&Xtp[i]<X2urban)&&(Ytp[i]>Y1urban&&Ytp[i]<Y2urban)){
	//	i=i-1;
	//	printf("%d renew\n",i);
	//}
	else{
		for(k=i-1;k>-1;k--)
		//if(i!=j)
			if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
				i=i-1;
	}
}/*
for(i=nUser+ceil((Ntp-nUser)*(SubUser+UrbanUser));i<Ntp;i++){	//randomly generate
	ShortestDistance=1000000000;
	Xtp[i]=rand()%(Px+1);
	Ytp[i]=rand()%(Py+1);
	for(j=0;j<Tbs-mNbs;j++){
		tempDistance[j]=space*distance(Xbs[j],Ybs[j],Xtp[i],Ytp[i]);
		if(tempDistance[j]<ShortestDistance)
			ShortestDistance=tempDistance[j];
		//printf("%lf\n",ShortestDistance);
	}
	if(ShortestDistance>Rbs){
		i=i-1;
		printf("%d renew\n",i);
	}
	else if((Xtp[i]>X1suburban&&Xtp[i]<X2suburban)&&(Ytp[i]>Y1suburban&&Ytp[i]<Y2suburban)){
		i=i-1;
		printf("%d renew\n",i);
	}
	else{
		for(k=i-1;k>-1;k--)
		//if(i!=j)
			if(Xtp[i]==Xtp[k]&&Ytp[i]==Ytp[k])
				i=i-1;
	}
}*//*
for(i=41;i<=64;i++){
	Xtp[i]=(rand()%(376))+250;
	Ytp[i]=(rand()%(376))+250;
}*/
return;
}

