/* error.C
   error handling

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: error.C 91 2007-05-06 17:47:45Z bakker $
*/

#define SVN_DATE "$Date: 2007-05-06 19:47:45 +0200 (Sun, 06 May 2007) $"

#ifdef USE_MPI
  #include"mympi.h"
#endif

// #include <sys/types.h>
// #include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>

#include "error.h"
#include "libprf/prf.h"

static int
  firstTime   = 1,
  reportLevel = 0,   // verbosity level
  nodeID      = 0;   // required to be MPI compatible

static char
  *fName = NULL;

void initErrorMsg(int level,const char *fbase) {
  stdprf ();    /* set standard functions */
  extprf ();    /* set extended standard functions */
  fltprf ();    /* set floating standard functions */

  reportLevel = level;

#ifdef USE_MPI
  nodeID = mpi_rank;
#endif

  if (fbase) {
    initErrorFilename(fbase);
  }
} /* initErrorMsg() */

void initErrorFilename(const char *fbase) {
  if (fbase && (strlen(fbase) > 0)) {
    char *cmd = (char *) malloc(sizeof(char)*(strlen(fbase) + 100));

    if (fName) {
      free(fName);
    }

    fName = (char *) malloc(sizeof(char)*(strlen(fbase)+10));
    if (fName) {
#ifdef USE_MPI
      nodeID = mpi_rank;
      sprintf(fName,"%s.%03d.msg",fbase,mpi_rank);
#else
      sprintf(fName,"%s.msg",fbase);
#endif

      sprintf(cmd,"rm -f %s",fName);
      system(cmd);
      free(cmd);
      firstTime = 1;
    } else {
      fprintf(stderr,"ERROR in initErrorFilename: %s (%s)\n",
              "Insufficient memory to allocate filename",fbase);
    }
  }
} /* initErrorFileName() */

void setReportLevel(int level) {
  reportLevel = level;
}

void writeError(ErrorMode m,ErrorType t,const char* fmt, ...) {
  if ((m == errModeAll) ||
      ((m == errModeMaster) && (nodeID == 0)) ||
      ((m == errModeSlave)  && (nodeID  > 0))) {

    va_list ap;
    char    str[4096],mpiStr[20];

    va_start(ap,fmt);

    switch (t) {
    case errMessage: sprf(str,"MSG:.. "); break;
    case errWarning: sprf(str,"WAR:.. "); break;
    case errGeneral: sprf(str,"ERR:.. "); break;
    case errFatal:   sprf(str,"FATAL: "); break;
    }

#ifdef USE_MPI
    sprintf(mpiStr,"%2d ",mpi_rank);
    sprf(str,mpiStr);
#endif

    sprfv(str, fmt, &ap);
    va_end(ap);

    fprintf(stderr,"%s\n",str);
    fflush(stderr);

    if (fName) {
        FILE *ofp = fopen(fName,(firstTime==1?"w":"a"));

        if (ofp) {
            fprintf(ofp,"%s\n",str);
            fclose(ofp);
        } else if (firstTime == 1) {
            int myerr = errno;
            fprintf(stderr,
                    "writeError cannot open %s (%d)\n%s\n",
                    fName,myerr,strerror(myerr));
        }
        firstTime = 0;
    }
  }
  if (t == errFatal) {
#ifdef USE_MPI
    fprintf(stderr,"FATAL ERROR reported on node %d\n",mpi_rank);
    MPI::COMM_WORLD.Abort(1);
#endif
    exit(1);
  }
}

void writeError(int level,ErrorMode m,ErrorType t,const char* fmt, ...) {
  if ((m == errModeAll) ||
      ((m == errModeMaster) && (nodeID == 0)) ||
      ((m == errModeSlave)  && (nodeID  > 0))) {

    if ((level <= reportLevel) || (t == errFatal)) {
      va_list ap;
      char    str[4096],mpiStr[20];

      va_start(ap,fmt);

      switch (t) {
      case errMessage: sprf(str,"MSG:.. "); break;
      case errWarning: sprf(str,"WAR:.. "); break;
      case errGeneral: sprf(str,"ERR:.. "); break;
      case errFatal:   sprf(str,"FATAL: "); break;
      }

#ifdef USE_MPI
      sprintf(mpiStr,"%2d ",mpi_rank);
      sprf(str,mpiStr);
#endif

      sprfv(str, fmt, &ap);
      va_end(ap);

      fprintf(stderr,"%s\n",str);
      fflush(stderr);

      if (fName) {
          FILE *ofp = fopen(fName,(firstTime==1?"w":"a"));

          if (ofp) {
              fprintf(ofp,"%s\n",str);
              fclose(ofp);
          } else if (firstTime == 1) {
              int myerr = errno;
              fprintf(stderr,
                      "writeError cannot open %s (%d)\n%s\n",
                      fName,myerr,strerror(myerr));
          }
          firstTime = 0;
      }
    }
  }
  if (t == errFatal) {
#ifdef USE_MPI
    fprintf(stderr,"FATAL ERROR reported on node %d\n",mpi_rank);
    MPI::COMM_WORLD.Abort(1);
#endif
    exit(1);
  }
}