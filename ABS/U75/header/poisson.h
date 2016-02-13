const int PoissonRandomNumber(const double lambda, int max_k){
	int k=0;	//counter
	
	double u;	//uniform random number
	double L = exp(-lambda);	//probability
	double p = 1;
	
	while (1){
		k++;
		u = (double) rand()/RAND_MAX;
		p = p * u;
		if (p < L)
			break;
	}
	if(k-1>max_k){
		printf("Poisson Random Number is too high!");
		exit(1);
	}
return (k-1);
}
