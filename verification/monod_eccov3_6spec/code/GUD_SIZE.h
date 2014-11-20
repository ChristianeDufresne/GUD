C $Header$
C $Name$

#ifdef ALLOW_GUD

CBOP
C    !ROUTINE: GUD_SIZE.h
C    !INTERFACE:
C #include GUD_SIZE.h

C    !DESCRIPTION:
C Contains dimensions and index ranges for cell model.

      integer nplank, nGroup, nlam, nopt
      integer iMinPhoto, iMaxPhoto, nPhoto
      integer iMinPrey, iMaxPrey
      integer iMinPred, iMaxPred
      integer nNQuota, nPQuota, nSiQuota, nFeQuota, nChl
      parameter(nlam=1)
      parameter(nopt=1)
      parameter(nplank=8)
      parameter(nGroup=7)
      parameter(iMinPhoto=1, iMaxPhoto=6)
      parameter(nPhoto=iMaxPhoto+1-iMinPhoto)
      parameter(iMinPrey=1, iMaxPrey=iMaxPhoto)
      parameter(iMinPred=iMaxPhoto+1, iMaxPred=nplank)
C     to turn off nitrogen quota, set to zero
      parameter(nNQuota=0)
      parameter(nPQuota=0)
      parameter(nSiQuota=0)
      parameter(nFeQuota=0)
      parameter(nChl=nPhoto)

CEOP
#endif /* ALLOW_GUD */
