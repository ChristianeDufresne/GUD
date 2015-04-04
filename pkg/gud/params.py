'''
&locals

_RL pday = 86400.0 _d 0
_RL pyear = 360 _d 0 * pday


&GUD_CONSTANTS-
_RL rad2deg = 180 _d 0 / PI


#ifdef GUD_ALLOW_CARBON

&CARBON_CONSTANTS-

_RL Pa2Atm  = 1.01325 _d 5
_RL ptr2mol = 1.0 _d -3                  ! convert ptracers (in mmol/m3) to mol/m3

! set up coefficients for DIC chemistry
! define Schmidt no. coefficients for CO2
_RL sca1 = 2073.1 _d 0
_RL sca2 = -125.62 _d 0
_RL sca3 =    3.6276 _d 0
_RL sca4 =   -0.043219 _d 0

! define Schmidt no. coefficients for O2
! based on Keeling et al [GBC, 12, 141, (1998)]
_RL sox1 = 1638.0 _d 0
_RL sox2 =  -81.83 _d 0
_RL sox3 =    1.483 _d 0
_RL sox4 =   -0.008004 _d 0

! coefficients for determining saturation O2
_RL oA0 =  2.00907 _d 0
_RL oA1 =  3.22014 _d 0
_RL oA2 =  4.05010 _d 0
_RL oA3 =  4.94457 _d 0
_RL oA4 = -2.56847 _d -1
_RL oA5 =  3.88767 _d 0
_RL oB0 = -6.24523 _d -3
_RL oB1 = -7.37614 _d -3
_RL oB2 = -1.03410 _d -2
_RL oB3 = -8.17083 _d -3
_RL oC0 = -4.88682 _d -7

#endif


&GUD_PARAMS

! general parameters (same for all plankton)

INTEGER gud_seed = 0

INTEGER iDEBUG = 1
INTEGER jDEBUG = 1
INTEGER kDEBUG = 1

LOGICAL gud_linFSConserve = .FALSE.

LOGICAL gud_read_phos = .FALSE.

!_RL phymin = 1 _d -10
!_RL phymin = 1 _d -50
!_RL phymin = 1 _d -20
_RL phymin = 120 _d -20                   ! minimum phyto (below which grazing and mortality doesn't happen)

_RL katten_w   = 4. _d -2                 ! atten coefficient water(m^-1)
_RL katten_chl = 4. _d -2                 ! atten coefficient chl ((mmol chl/m3)-1)

! par parameters
_RL parfrac = 0.4 _d 0                    ! fraction Qsw that is PAR
_RL parconv = 1. _d 0/0.2174 _d 0         ! conversion from W/m2 to uEin/m2/s

_RL tempnorm     = 0.3 _d 0               ! set temperature function (was 1.0)
_RL TempAeArr    = -4000.0 _d 0           ! coefficients for TEMP_VERION 2 (pseudo-Arrhenius)
_RL TemprefArr   = 293.15 _d 0
_RL TempCoeffArr = 0.5882 _d 0

_RL alpfe       = 0.04 _d 0               ! solubility of Fe dust
_RL scav        = 0.4 _d 0/pyear          ! iron chem scavenging rate (s-1)
! scav          = 4.4 _d -3/pday          ! for quota model
_RL ligand_tot  = 1. _d -3                ! total ligand (mol m-3)
_RL ligand_stab = 2. _d 5                 ! ligand stability rate ratio = [FeL]/[Fe'][L'] (m3 mol-1)
_RL freefemax   = 0.4 _d -3               ! max free fe (mol m-3)

_RL scav_rat      = 0.005 _d 0 / pday     ! scavening
_RL scav_inter    = 0.079 _d 0
_RL scav_exp      = 0.58 _d 0
_RL scav_R_POPPOC = 1.1321 _d -4

                                          ! fe sed
_RL depthfesed    = -1.0 _d 0                    ! depth above which to add sediment source (was -1000)
_RL fesedflux     = 1.0 _d 0 * 1.0 _d -3 / pday  ! iron flux (mmol/m2/s)
_RL fesedflux_pcm = 0.68 _d 0 * 1.0 _d -3        ! iron flux (mmol/m3/s) per mmol POC/m3/s
_RL R_CP_fesed    = 106 _d 0

                                          ! oxidation rates for ammonium and nitrite
_RL Knita   = 1.0 _d 0/(.50 _d 0*pday)    ! i.e. Knita ...    NH4 -> NO2
_RL Knitb   = 1.0 _d 0/(10.0 _d 0*pday)   ! i.e. Knitb ...    NO2 -> NO3
_RL PAR_oxi = 10. _d 0                    ! critical light level (muEin/m2/s) after which oxidation starts

_RL Kdoc  = 1.0 _d 0/(100.0 _d 0*pday)    ! DOM remin rates
_RL Kdop  = 1.0 _d 0/(100.0 _d 0*pday)
_RL Kdon  = 1.0 _d 0/(100.0 _d 0*pday)
_RL KdoFe = 1.0 _d 0/(100.0 _d 0*pday)

_RL KPOC  = 1.0 _d 0/(50.0 _d 0*pday)     ! Particulate detritus remin rates
_RL KPOP  = KPOC                          ! z*                             = wx_sink/Kremin_X
_RL KPON  = KPOC                          ! for e-folding length scale, z* = 300 m
_RL KPOFe = KPOC                          ! choose Kremin_X                = 1/30 day-1, wx_sink = 10 m day-1
_RL KPOSi = 1.0 _d 0/(300.0 _d 0*pday)

_RL wC_sink  = 10.0 _d 0/pday             ! sinking rate for particulate matter (m/s)
_RL wP_sink  = wC_sink
_RL wN_sink  = wC_sink
_RL wFe_sink = wC_sink
_RL wSi_sink = wC_sink
_RL wPIC_sink = 15.0 _d 0/pday

_RL Kdissc = 1.0 _d 0/(300.0 _d 0*pday)

#ifdef GUD_ALLOW_CARBON
_RL gud_atmos_pCO2 = 278 _d -6
_RL R_OP           = 170 _d 0
_RL R_OC           = 170 _d 0/120.0 _d 0
_RL m3perkg        = 1. _d 0 / 1024.5 _d 0
#endif

_RL diaz_ini_fac = 1 _d 0                 ! reduce tracer concentrations by this factor on initialization


! denitrification

_RL O2crit    = 6.0 _d 0                  ! critical oxygen for O2/NO3 remineralization (Lipschultz et al 1990, DSR 37, 1513-1541)
_RL denit_NP  = 120.0 _d 0                ! ratio of n to p in denitrification process
_RL denit_NO3 = 104.0 _d 0                ! ratio no3 used relative to all n in denitrification process
_RL NO3crit   = 1. _d -2                  ! critical nitrate below which no denit (or remin) happens



! these should probably be traits
                                          ! from quota package
_RL PARmin    = 0.1 _d 0                  ! minimum light for photosynthesis; for non-Geider: 1.0
_RL chl2nmax  = 3.00 _d 0                 ! Moore 2002 (labelled as GD98)
_RL synthcost = 0.0 _d 0                  ! cost of biosynthesis (mmol C/mmol N)
                                          ! 2.33 for quota model

                                          ! new combined parameters
_RL expPref  = 1.0 _d 0
_RL expPalat = 1.0 _d 0

_RL palat_min = 0 _d 0                    ! was 1D-4 in quota; smaller palat are set to 0

_RL inhib_graz     = 1.0 _d 0             ! for quota-style grazing ((mmol C m-3)-1)
_RL inhib_graz_exp = 0.0 _d 0             ! exponent 0.0 turns it off

_RL hollexp = 1.0 _d 0                    ! grazing exponential 1= holling 2, 2=holling 3

_RL phygrazmin = 120 _d -10               ! minimum total prey conc


! TO BE IMPLEMENTED

_RL depthdenit = 185.0 _d 0               ! depth for denitrification relaxation to start



#ifdef GUD_ALLOW_RADTRANS

&GUD_RADTRANS_PARAMS

CHARACTER*80 gud_waterAbsorbFile = ' '
CHARACTER*80 gud_phytoAbsorbFile = ' '
CHARACTER*80 gud_particleAbsorbFile = ' '

LOGICAL gud_useModelCalendar = .FALSE.

INTEGER gud_radtrans_kmax = Nr

_RL gud_part_size_P = 1 _d -15   ! mmol P per particle

_RL(nlam+1) gud_waveband_edges   = -1.0 _d 0
_RL(nlam)   gud_waveband_centers = -1.0 _d 0          ! representative wavelengths

_RL gud_radmodThresh = 1 _d -4
_RL gud_rmus         = 1.0/0.83 _d 0
_RL gud_rmuu         = 1.0/0.4 _d 0
_RL gud_bbmin        = 0.0002 _d 0
_RL gud_bbw          = 0.5 _d 0

_RL gud_lambda_aCDOM = 450.0 _d 0
_RL gud_Sdom         = 0.014 _d 0
_RL gud_aCDOM_fac    = 0.2 _d 0

#endif



#ifdef GUD_ALLOW_CDOM

&GUD_CDOM_PARAMS

_RL fracCDOM   = 2. _d 0 / 100. _d 0
_RL CDOMdegrd  = 1. _d 0 / (200 _d 0*pday)
_RL CDOMbleach = 1. _d 0 / (15 _d 0*pday)
_RL PARCDOM    = 20. _d 0
_RL R_NP_CDOM  = 16. _d 0
_RL R_FeP_CDOM = 1. _d -3
_RL R_CP_CDOM  = 120. _d 0
_RL CDOMcoeff  = .1 _d -1  / 1. _d -4

#endif



&GUD_DEPENDENT_PARAMS-
! dependent parameters (computed here or in gud_init_fixed)
INTEGER kMinFeSed = 1
INTEGER kMaxFeSed = 0
'''

import parser
globals().update(parser.parse(__doc__))

