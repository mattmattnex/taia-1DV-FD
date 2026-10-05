
!> Module that contains the parameters read in the "runparam.i" file, i.e general parameters of the simulation
!! (start and end times, time steps, etc).



      MODULE comrunmod
  !-----------------------------------------------------------------
  !
     Implicit None
  !        ----------------------------
  ! =====> DIMENSIONNEMENT DES TABLEAUX
  !        ----------------------------
  !
  ! jpzt = nbr maximum de niveaux
  ! jptract = nbr maximum de traceurs biologiques
  ! jpfluxt = nbr maximum de flux biologiques
  ! jptemps = nbr de sauvegardes maximum dans le run vis a vis du temps
  !
     integer,parameter :: jpzt=42,jpfluxt=50,jptemps=14600
  !simu 5 ans   integer,parameter :: jpzt=42,jpfluxt=50,jptemps=7300
  !simu 10 ans   integer,parameter :: jpzt=42,jpfluxt=50,jptemps=14600
  !simu 30 ans   integer,parameter :: jpzt=42,jpfluxt=50,jptemps=48180
     integer :: jptract
  !
  !        ---------------------------------------
  ! =====> PARAMETRES GENERAUX LUS DANS runparam.i
  !        ---------------------------------------
  !
!!! cnmlog = Directory de sauvegarde des fichiers binaires output.
  !           nmlog = nbr de character de cnmlog
  !
     character*50 cnmlog
      integer::  nmlog
  !
!!! cdir = directory de travail local - adresse absolue
  !           ndir=nbr de character de cdir
  !           ndir0=nbr de character du directory local (dernier dans
  !                 l'adresse)
  !           croot=racine des fichiers output binaire
  !               =cnmlog(1:nmlog)//'/'//cdir(ndir0:ndir)//'_'
  !           nroot=nbr de character de croot
  !
     character*120 cdir,croot
     integer::  ndir,ndir0,nroot
  !
!!! ctit,csubtit = Titre et sous titre du run (character*80)
     character*80 ctit,csbtit
  !
!!! cgrille = nom fichier grille (93)
  !
     character*80 cgrille
  !
!!! dt = pas de temps en MINUTES
!!! tdebut0 = temps initiale en JOURS (par rapport a 1 Janvier 0h)
!!! tfin0 = duree totale du run en JOURS
  !           dts = pas de temps en SECONDES
  !           tdebut = temps initial en MINUTES
  !           tfin = duree total du run en MINUTES
  !           npurday = nbr de pas de temps par jour (1440/dt)
  !           nstep = pas de temps courant
  !           time = temps ecoule durant le run en MINUTES
  !           timep = temps courant en MINUTES (time+tdebut)
  !           timeyr = annee courante
  !           nmonth = mois courant
  !           njour = jour courant du mois courant
  !           timeday = jour courant dans l'annee courante
  !           ntimeday = INT(timeday)
  !           timeminu = minute courante dans la journee courante
  !           ntimeminu = INT(timeminu)
  !
     real:: dt,tdebut0,tfin0
     real:: dts,tdebut,tfin,time,timep,timeyr,timeday,timeminu
     integer:: npurday,nstep,nmonth,njour,ntimeday,ntimeminu
  !
!!! nzt = nbr de niveaux verticaux
  !   nbrbio = nbr de traceurs biogeochimiques
!!! nbrprono = nbr de traceurs biogeochimiques pronostiques
!!! nbrdiag = nbr de traceurs biogeochimiques diagnostiques
!!! CNMTRA = character*3: nom des traceurs biogeochimiques
!!! nbrflux = nbr de flux biogeochimiques
!!! CNMFLX = character*6: nom des flux biogeochimiques
  !
     integer:: nzt,nbrprono,nbrdiag,nbrflux
     character(3),allocatable:: CNMTRA(:)
     character*6 CNMFLX(jpfluxt)
     integer:: nbrbio
  !
!!! timesave0 = periode de sauvegarde em MINUTES
!!! nwrite = Ecritures en nbr de timesave0
  !        timesave = periode de sauvegarde en HEURES
  !        nsave = nbr de pas de temps par periode de sauvegarde
  !        nsa = nbr de pas de temps courant apres derniere sauvegarde
  !        nind = indice de sauvegarde (temporelle) pour les champs 1D
  !        nind0 = indice courant par a l'interieur de nind
  !        nwrite0 = nbr de pas de temps apres derniere ecriture
  !
     real:: timesave0
     integer:: nwrite
     real:: timesave
     integer:: nsave,nsa,nind,nind0,nwrite0
  !
!!! mtke    = flag sauvegarde TKE
!!! mt      = flag sauvegarde T  
!!! ms      = flag sauvegarde S  
!!! mdifft  = flag sauvegarde Diffusion Turbulente
!!! mlength = flag sauvegarde mixing length, Brunt Vaissaila, shear
!!! mtket   = flag tendances equation TKE
!!! mtt     = flag tendances equation T
!!! mst     = flag tendances equation S
  !       nwdyn = nbr de traceurs dynamiques sauvegardes
  !       nwdynfl = nbr de traceurs dynamiques dont tendances sauvegardees
  !
     real:: mtke,mt,ms,mdifft,mlength,mtket,mtt,mst
     integer:: nwdyn,nwdynfl
 !
!!! MSTOCK  = flag sauvegarde des variables biogeochimiques
!!! MFLUX   = flag sauvegarde flux biogeochimiques
!!! MBIODIF = flag sauvegarde des flux diffusifs biogeochimiques pronostiques
!!! mbiosed = flag sauvegarde flux de sedimentation
  !        nwbio = nbr de traceurs biologiques sauvegardes
  !        nwbiofl = nbr de flux biologiques sauvegardes
  !        nwbiodif = nbr de flux diffusifs biologiques sauvegardes
  !        nwbiosed = nbr de flux de sedimentation sauvegardes
  !
     real:: MFLUX(jpfluxt),mbiosed
     real,allocatable:: MSTOCK(:),MBIODIF(:)
     integer:: nwbio,nwbiofl,nwbiodif,nwbiosed
  !
!!! mflagwrite = flag d'ecriture (1 ==> ecriture ; 0 ==> entetes seuls)
  !
     integer:: mflagwrite
  !
  !        --------------------------
  ! =====> PARAMETRES COMPLEMENTAIRES
  !        --------------------------
  !
  ! DEPT = Profondeurs niveau traceur.
  ! DEPW = Profondeurs niveau vitesse verticale
  ! TMASK = Masque traceur
  ! E3T = Distances entre points vitesse verticale
  ! E3W = Distances entre points traceur
  !
     real(8) ::  DEPT(jpzt),DEPW(jpzt),TMASK(jpzt)
     real ::  E3T(jpzt),E3W(jpzt)
  !
  ! nstepneg = nbr de pas de temps avec au moins  une variable  < 0
  ! nstepneg0 = nbr de passage dans boucle (variable biologique < 0)
  !
     integer ::  nstepneg,nstepneg0
  !
  ! nundyn = unite sauvegarde variables dynamiques
  ! nundynfl = unite sauvegarde flux dynamiques
  ! nundyn1d = unite sauvegarde champs dynamiques 1D
  ! nunbio = unite sauvegarde variables biogeochimiques
  ! nunbiofl = unite sauvegarde flux biogeochimiques
  ! nunbio1d = unite sauvegarde champs biogeochimiques 1D
  !
     integer ::  nundyn,nundynfl,nundyn1d,nunbio,nunbiofl,nunbio1d
  !
  ! mflagrun = flag d'erreur dans definition des parametres de runparam.i
  ! mflagtke = flag d'erreur dans definition des parametres de dynparam.i
  ! mflagflux= flag d'erreur dans lecture/definition des forcages
  ! mflagbio = flag d'erreur dans initialisation modele biologique
  ! 
     integer:: mflagrun,mflagtke,mflagflux,mflagbio
  !
  ! Espace de travail
  !
!MB  real:: workb,WORK(jpzt),WORK1(jpzt)
     real(8):: workb,WORK(jpzt),WORK1(jpzt)
  !
  ! Data
     real:: day,daymn,hour,hourm,xyear 
     real::  mday,mdaymn,mhour,mhourm,mxyear 
     real:: temp0 
     real:: valun
     real,parameter::ppi=3.14159265,ppdtr=ppi/180.,pprtd=1./ppdtr
 !
     data day,daymn,hour,hourm,xyear /24.,1440.,3600.,60.,365./
     data mday,mdaymn,mhour,mhourm,mxyear /24,1440,3600,60,365/
     data temp0 /273.15/
     data valun /1./
!------------------------------------------------------------------
     END   MODULE comrunmod
!-----------------------------------------------------------------
