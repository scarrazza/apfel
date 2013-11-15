#include <cmath>
#include "common.h"


double ComputeAVG(int n, int ix, double **x)
{
  double sum = 0.0;
  for (int i = 0; i < n; i++)
    sum += x[i][ix];
  return sum / n;
}


double ComputeStdDev(int n, int ix, double **x)
{
  double sum = 0.0;
  double avg = ComputeAVG(n, ix, x);
  for (int i = 0; i < n; i++)
    sum += (x[i][ix]-avg)*(x[i][ix]-avg);

  sum /= n-1;

  return sqrt(sum);
}
