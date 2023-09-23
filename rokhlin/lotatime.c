/*
 *
 */

#include <stdio.h>
#include <time.h>

/*
*
*
*
*
*/

double clotatim_(void)
{
        long clock(void);
	double t;

        /* 
        *      return to the user the execution time for this process
        */

        t=clock();
        t=t/CLOCKS_PER_SEC;

        return t;
}
