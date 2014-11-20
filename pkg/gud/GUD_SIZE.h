C $Header$
C $Name$

#ifdef ALLOW_GUD

CBOP
C    !ROUTINE: GUD_SIZE.h
C    !INTERFACE:
C #include GUD_SIZE.h

C    !DESCRIPTION:
C Contains dimensions and index ranges for gud model.

      integer nplank, nGroup, nlam, nopt
      integer iMinPhoto, iMaxPhoto
      integer iMinPrey, iMaxPrey
      integer iMinPred, iMaxPred
      integer nNQuota, nPQuota, nSiQuota, nFeQuota, nChl
      parameter(nplank=4)
      parameter(nGroup=5)
      parameter(nlam=1)
      parameter(nopt=1)
      parameter(iMinPhoto=1, iMaxPhoto=2)
      parameter(iMinPrey=1, iMaxPrey=iMaxPhoto)
      parameter(iMinPred=iMaxPhoto+1, iMaxPred=nplank)
      parameter(nPhoto=iMaxPhoto+1-iMinPhoto)
C     to turn off nitrogen quota, set to zero
      parameter(nNQuota=nplank)
      parameter(nPQuota=nplank)
      parameter(nSiQuota=nplank)
      parameter(nFeQuota=nplank)
      parameter(nChl=iMaxPhoto)

CEOP
#endif /* ALLOW_GUD */
