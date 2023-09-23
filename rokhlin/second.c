/*
        Returns the time in seconds used by the process.
        The system call he uses is getrusage(2), which returns a
        structure containing various pieces of resource usage information.
 */
/*********************************************************************/
#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>


double second_()
{
        static struct rusage temp;
        double foo, foo1;

        getrusage(RUSAGE_SELF,&temp)    ;
        foo     = temp.ru_utime.tv_sec  ;       /* seconds */
        foo1    = temp.ru_utime.tv_usec ;       /* uSecs */
    return  foo  + (foo1/1000000.0)         ;       /* milliseconds */
}
