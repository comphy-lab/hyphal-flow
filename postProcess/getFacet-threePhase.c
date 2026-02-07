/**
# getFacet-threePhase.c

Extract interface facets from a saved snapshot for either the drop
(`f1`) or film (`f2`) phase.
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[], f1[], f2[];
char filename[80];
bool includeCoat;

/**
## main()

Usage:
`./getFacet-threePhase snapshot true|false`

When `true`, facets are extracted from `f1`; otherwise `f2` is used.
*/
int main(int a, char const *arguments[]){
  sprintf(filename, "%s", arguments[1]);

  /**
  Interpret the second argument as a boolean selector.
  */
  if (strcmp(arguments[2], "true") == 0) {
    includeCoat = true;
  } else {
    includeCoat = false;
  }

  restore (file = filename);

  if (includeCoat == true){
    foreach(){
      f[] = f1[];
    }
  } else {
    foreach(){
      f[] = f2[];
    }
  }

  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  fclose (fp);
}
