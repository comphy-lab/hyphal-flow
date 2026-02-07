/**
# getData-elastic-nonCoalescence.c

Sample scalar diagnostics from a Basilisk snapshot and print tabulated
pointwise values for downstream plotting.

## Computed Fields
- `D2c`: log-scaled viscous dissipation proxy
- `vel`: velocity magnitude
- `trA`: trace-like conformation/stress metric
*/

#include "utils.h"
#include "output.h"

scalar f1[], f2[];
vector u[];
symmetric tensor tau_p[];
scalar tau_qq[];

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay, Oh1, Oh2, Oh3;

scalar D2c[], vel[], trA[];
scalar * list = NULL;

/**
## main()

Usage:
`./getData-elastic-nonCoalescence snapshot xmin ymin xmax ymax ny Oh1 Oh2 Oh3`
*/
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);

  Oh1 = atof(arguments[7]);
  Oh2 = atof(arguments[8]);
  Oh3 = atof(arguments[9]);

  list = list_add (list, D2c);
  list = list_add (list, vel);
  list = list_add (list, trA);

  /**
  Restore the snapshot and evaluate derived fields on cell centers.
  */
  restore (file = filename);

  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/y);
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = 2*( Oh1*clamp(f1[], 0., 1.) + Oh2*clamp(f2[], 0., 1.) + Oh3*clamp(1-f1[]-f2[], 0., 1.) )*D2;
    
    if (D2c[] > 0.){
      D2c[] = log(D2c[])/log(10);
    } else {
      D2c[] = -10;
    }

    vel[] = sqrt(sq(u.x[])+sq(u.y[]));

    trA[] = (tau_p.x.x[] + tau_p.y.y[] + tau_qq[])/2.0;

    if (trA[] > 0.){
      trA[] = log(trA[])/log(10);
    } else {
      trA[] = -10;
    }

  }

  FILE * fp = ferr;
  Deltay = (double)((ymax-ymin)/(ny));
  nx = (int)((xmax - xmin)/Deltay);
  Deltax = (double)((xmax-xmin)/(nx));
  len = list_len(list);
  /**
  Interpolate all diagnostics onto a uniform `nx x ny` sampling grid.
  */
  double ** field = (double **) matrix_new (nx, ny+1, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  fclose (fp);
  matrix_free (field);
}
