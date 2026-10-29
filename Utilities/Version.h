#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "7586ffa25948c6ce3bc6c238cda30ad8d35f60a0"
#define GIT_REFSPEC "refs/heads/master"
#define GIT_LOCAL_STATUS "CLEAN"

#define SPLISHSPLASH_VERSION "2.19.0"

#ifdef DL_OUTPUT

#endif

#endif
