Include "Parameter.pro";
Function{
  b = {0.0, 0.14285714285714285, 0.2857142857142857, 0.42857142857142855, 0.5714285714285714, 0.7142857142857142, 0.8571428571428571, 1.0} ;
  mu_real = {3000.0, 2978.4390619047617, 2967.756242857143, 1.0, 1.0, 1.0, 1.0, 1.0} ;
  mu_imag = {1.0, 343.4662714285714, 477.2211857142857, 486.9375, 486.9375, 486.9375, 486.9375, 486.9375} ;
  mu_imag_couples = ListAlt[b(), mu_imag()] ;
  mu_real_couples = ListAlt[b(), mu_real()] ;
  f_mu_imag_d[] = InterpolationLinear[Norm[$1]]{List[mu_imag_couples]};
  f_mu_real_d[] = InterpolationLinear[Norm[$1]]{List[mu_real_couples]};
  f_mu_imag[] = f_mu_imag_d[$1];
  f_mu_real[] = f_mu_real_d[$1];
 }  