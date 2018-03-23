#include <sys/types.h>
#include <unistd.h>

void newdir_ (int* myid)
{
  char path[100];
  if (*myid > 0)
    {
      sprintf (path, "./%d", *myid);
/*      printf ("myid = %d, path = %s\n", *myid, path); */
      chdir (path);
      getcwd (path, 100);
/*      printf ("myid = %d, cwd = %s\n", *myid, path); */
    }

}

void xgetpid_ (int* pid)
{
  *pid = getpid ();
}
