Function{
  // Liste von Lukas hinterlegen
  // Mu_imag muss noch aus mur und


  // --------------------------
  // TDK N95 100 Celsius

  N95_b = {
  0, 0.1, 0.2, 0.3, 0.4
  } ;

  N95_mu_imag = {
  200.0, 400.0, 600.0, 800.0, 1000.0
  } ;

  N95_mu_imag_couples = ListAlt[N95_b(), N95_mu_imag()] ;

  f_N95_mu_imag[] = InterpolationLinear[Norm[$1]]{List[N95_mu_imag_couples]};
}
