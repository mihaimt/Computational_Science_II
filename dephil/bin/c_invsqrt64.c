double c_invsqrt64_( double number )
{
	long i;
	double x2, y;
	const double threehalfs = 1.5F;
 
	x2 = number * 0.5F;
	y  = number;
	i  = * ( double * ) &y;                       // evil floating point bit level hacking
	i  = 0x5fe6eb50c7b537a9 - ( i >> 1 );		// wtf?!
	y  = * ( double * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
//      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
 
	return y;
}
