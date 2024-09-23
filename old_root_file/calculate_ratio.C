void calculate_ratio(){
	double x1 = 2.666;
	double x2 = 2.635;
	
	double err1 = 0.138;
	double err2 = 0.149;

	double ratio;
	double errratio;

	ratio = (x1 - 2.4955)/(x2 - 2.4955);

	errratio = ratio * sqrt(pow((err1/x1),2) + pow((err2/x2),2));

	cout << "ratio is " << ratio << "errratio is " << errratio << endl;

}