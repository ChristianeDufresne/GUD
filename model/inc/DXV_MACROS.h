C $Header$
C
C     /==========================================================\
C     | DXV_MACROS.h                                             |
C     |==========================================================|
C     | These macros are used to reduce memory requirement and/or|
C     | memory references when variables are fixed along a given |
C     | axis or axes.                                            |
C     \==========================================================/

#ifdef DXV_CONST
#define  _dxV(i,j,bi,bj) dxV(1,1,1,1)
#endif

#ifdef DXV_FX_ONLY
#define  _dxV(i,j,bi,bj) dxV(i,1,bi,1)
#endif

#ifdef DXV_FY_ONLY
#define  _dxV(i,j,bi,bj) dxV(1,j,1,bj)
#endif

#ifndef _dxV
#define  _dxV(i,j,bi,bj) dxV(i,j,bi,bj)
#endif
