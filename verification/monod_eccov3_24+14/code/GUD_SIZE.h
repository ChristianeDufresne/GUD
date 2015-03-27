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
      integer nChl
      integer nPPplank
      integer nGRplank
      parameter(nlam=1)
      parameter(nopt=1)
      parameter(nplank=38)
      parameter(nGroup=8)
      parameter(iMinPhoto=1, iMaxPhoto=24)
      parameter(nPhoto=iMaxPhoto+1-iMinPhoto)
      parameter(iMinPrey=1, iMaxPrey=iMaxPhoto)
      parameter(iMinPred=iMaxPhoto+1, iMaxPred=nplank)
C     to turn off nitrogen quota, set to zero
      parameter(nChl=nPhoto)
      parameter(nPPplank=0)
      parameter(nGRplank=0)

CEOP
#endif /* ALLOW_GUD */
