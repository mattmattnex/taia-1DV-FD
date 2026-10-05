 Module mod_varphy_coupl
 Use comrunmod
 IMPLICIT None

 character(LEN=80):: crapbio,cbioinit,coptique

 real:: rdif,aphymax,bstar,cstar,SRDIR(jpzt,19),&
 &SRDIF(jpzt,19), &
 &RDIR(19,4),XKW(jpzt),CHI(jpzt), &
 &EMO(jpzt),APHY(jpzt),AW(jpzt),ASTAR(jpzt),BW(jpzt)
!
 real:: trpar,SR0P(jpzt),SR0M(jpzt),SRW(jpzt,jpzt)
 integer, parameter:: nfrog = 5

!! TENEUR = concentrations en cours des traceurs
!! BIOAVER = moyenne des teneurs pour sauvegarde
!! FLUXBIO = moyenne des flux biologiques pour sauvegarde
!! DIFFBIO = moyenne des flux diffusifs pour sauvegarde
!! SEDBIO = moyenne des flux de sedimentation pour sauvegarde
!! BIOAIR = flux air - mer des traceurs
!! FLUXAIR = moyenne des flux air - mer pour sauvegarde
!! DEPTHZE = moyenne de la profondeur de la couche euphotique 
!! XNEBUL = id pour nebulosite
!! XPAR0PLUS = id pour xpar0+
!! VPAR = PAR
!! VPUR = PUR
!!
  real(8),allocatable::TENEUR(:,:),DIFFBIO(:,:),&
  &  SEDFLUX(:,:), SEDBIO(:,:)
!
  real,allocatable::EXPORT(:,:)
  real,allocatable::BIOAVER(:,:),BIOAIR(:),FLUXAIR(:,:)
 
 real::  FLUXBIO(jpzt,jpfluxt)
  real(8),allocatable::   VITSED0(:),VITSED(:)
  real :: DEPTHZE(jptemps),XNEBUL(jptemps),XPAR0PLUS(jptemps)
  real :: VPAR(jpzt),VPUR(jpzt)
  integer,allocatable::  MNEG(:,:),MNEGB(:,:) ! indices de vecteurs

! barrier.n: added for computation of optpar etc.
 real(8) :: xdwl,excen,dcln,daylength,ahsun,hsun,azsun,xpar0p,&
  &               xpar0m
  End Module mod_varphy_coupl
