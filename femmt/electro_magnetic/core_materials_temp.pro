Include "Parameter.pro";
Function{
  b = {0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4} ;
  mu_real = {3000.0, 2993.0, 2986.0, 2978.0, 2972.0, 2966.0, 2969.5, 2973.0, 1487.000000000001, 1.0} ;
  mu_imag = {1.0, 202.0, 288.0, 360.0, 405.0, 450.0, 485.0, 520.0, 260.00000000000017, 0.0} ;
  mu_imag_couples = ListAlt[b(), mu_imag()] ;
  mu_real_couples = ListAlt[b(), mu_real()] ;
  f_mu_imag_d[] = InterpolationLinear[Norm[$1]]{List[mu_imag_couples]};
  f_mu_real_d[] = InterpolationLinear[Norm[$1]]{List[mu_real_couples]};
  f_mu_imag[] = f_mu_imag_d[$1];
  f_mu_real[] = f_mu_real_d[$1];
 }  