C $Header$
C $Name$

#ifdef ALLOW_GUD

CBOP
C    !ROUTINE: GUD_DIAGS.h
C    !INTERFACE:
C #include GUD_DIAGS.h

C    !DESCRIPTION:
C Contains indices into diagnostics array

      integer iPP
      integer iNfix
      integer iDenit
      integer iDenitN
      integer iPPplank
      integer iGRplank
      integer gud_nDiag
      PARAMETER(iPP=   1)
      PARAMETER(iNfix= 2)
      PARAMETER(iDenit=3)
      PARAMETER(iDenitN=4)
      PARAMETER(iPPplank=5)
      PARAMETER(iGRplank=iPPplank+nPPplank)
      PARAMETER(gud_nDiag=iGRplank+nGRplank-1)

CEOP
#endif /* ALLOW_GUD */
