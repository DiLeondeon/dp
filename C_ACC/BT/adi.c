#include "header.h"
//#include "solve_subs.h"
//#include "x_solve.h"

void adi()
{
  
  compute_rhs();

  x_solve();

  y_solve();

  z_solve();

  add();
}
