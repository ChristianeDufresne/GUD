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
      parameter(nlam=13)
      parameter(nopt=4)
      parameter(nplank=80)
      parameter(nGroup=5)
      parameter(iMinPhoto=1, iMaxPhoto=78)
      parameter(iMinPrey=1, iMaxPrey=iMaxPhoto)
      parameter(iMinPred=iMaxPhoto+1, iMaxPred=nplank)
      parameter(nPhoto=iMaxPhoto+1-iMinPhoto)
C     to turn off nitrogen quota, set to zero
      parameter(nChl=0)
      parameter(nPPplank=0)

CEOP
#endif /* ALLOW_GUD */
