/* rmem_compat.h — compatibility for R 4.5+ memory API changes */
#ifndef RMEM_COMPAT_H
#define RMEM_COMPAT_H

#include <R.h>
#include <R_ext/Memory.h>
#include <Rversion.h>

/* Map deprecated names to supported ones on R >= 4.5.
   On older R, Calloc/Free/Realloc already exist, so we do nothing. */
#ifndef Calloc
# if defined(R_VERSION) && R_VERSION >= R_Version(4, 5, 0)
#  define Calloc(n, t)       R_Calloc((n), t)
#  define Realloc(p, n, t)   R_Realloc((p), (n), t)
#  define Free(p)            R_Free(p)
# endif
#endif

#endif /* RMEM_COMPAT_H */