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
      integer gud_nDiag
      PARAMETER(iPP=   1)
      PARAMETER(iNfix= 2)
      PARAMETER(iDenit=3)
      PARAMETER(iDenitN=4)
      PARAMETER(gud_nDiag=iDenitN)

CEOP
#endif /* ALLOW_GUD */
