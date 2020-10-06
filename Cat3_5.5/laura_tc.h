/*
** svn $Id: basin.h 797 2016-05-11 01:53:51Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for LAURA_TC.
**
** Application flag:   LAURA_TC
** Input script:       ocean_laura_tc.in
*/


#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define SOLVE3D
#define NONLIN_EOS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define ANA_GRID
#define ANA_STFLUX
#define ANA_BTFLUX
#define SOLAR_SOURCE
!#define DIURNAL_SRFLUX
#define GLS_MIXING
#ifdef GLS_MIXING
# define RI_SPLINES
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif
# define BIO_FENNEL
# define DIAGNOSTICS_BIO
# define CARBON
# define OXYGEN
!# define pCO2_RZ
# define TALK_PROGNOSTIC
!#define UV_SMAGORINSKY
!#define TS_SMAGORINSKY
!#define NPZD_POWELL
# define ANA_SPFLUX
# define ANA_BPFLUX
!# define ANA_SRFLUX
!#define CONST_PAR
!#define ANA_DQDSST
!#define ANA_SST
