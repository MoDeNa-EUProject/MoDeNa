void momentsConverter(const state_type &y , const double t);

void momentsConverter(const state_type &y , const double t)
{
// momentsConverter - converts moments based on the unit volume of foam
// @param - const state_type &y -  vector of all the variables
// @param - const double t - time
	double kappa 		= 1.0 - (y[8]/(1.0+y[8]));
	double M[4] 		= {};
	M[0]				= kappa*y[7];
	M[1]				= kappa*y[8];
	M[2]				= kappa*y[9];
	M[3]				= kappa*y[10];

	ofstream MM[4];
	MM[0].open("../results/M0.txt", std::ios::app);
	MM[1].open("../results/M1.txt", std::ios::app);
	MM[2].open("../results/M2.txt", std::ios::app);
	MM[3].open("../results/M3.txt", std::ios::app);

	for (int i = 0; i < 4; i++)
	{
		MM[i] << t << '\t' << M[i] << '\n';
		MM[i].close();
	}
}
