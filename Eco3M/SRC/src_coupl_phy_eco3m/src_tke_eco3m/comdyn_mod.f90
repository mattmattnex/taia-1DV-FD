
!> Module that contains the parameters read in the "dynparam.i" file, i.e parameters
!! that are used in the TKE model.

      Module comdynmod
   
      Use comrunmod

      Implicit None
    !        ----------
    ! =====> PARAMETRES
    !        ----------
    !
    ! jpant = nbr maximum de jours dans le fichier de forcages atmospheriques
    !
    !
     Integer,parameter::jpant=365
    !
    !        ---------------------------------------
    ! =====> PARAMETRES GENERAUX LUS DANS dynparam.i
    !        --------------------------------------- !
    !!! ctsinit = nom fichier profils initiaux variables dynamiques (91)
    !!! cflux = nom fichier forcages atmospheriques (94)
    !!! crappel = nom fichier rappel profils T et S (95)
    !!! crapsurft,crapsurfs = nom fichier rappel surface T et S (92)
    !!! msurft,masirfs = flag de rappel en surface
    !!! xperflux = periode des forcages atmospheriques en HEURES dans cflux
    !!! xperfluxt = periode totale des forcages en JOURS dans cflux
    !!! xperflux0 = date initiale du 1er enregistrement dans cflux en JOURS
    !
      character*80 :: ctsinit,cflux,crappel,crapsurft,crapsurfs
      integer :: msurft,msurfs
      integer :: xperflux,xperfluxt,xperflux0
      integer :: nlengthl,ntketl,nttl,nstl 
    !
    !!! xlat = latitude (calcul de Coriolis dans tke)
    !!! g = cte de gravite
    !!! f = Coriolis
    !!! xdiff, xdiss = constantes utilisees dans tke
    !!! avmb = diffusion background quantite de mouvement
    !!! avtb = diffusion background traceurs
    !!! emin0 = TKE min a la surface (condition aux limites)
    !!! emin = TKE min dans la colonne d'eau
    !!! ecartmx = ecart en profondeur (m) minimale pour calcul de profondeur
    !          de melange
   !!! drag: drag coefficient
   !!! rabs = fraction solaire non penetrante (absorbe sur xsi1 et 1-rabs
   !            absorbe sur xsi2)
    !!! xsi1, xsi2 = longueur d'extinction dans TKE en IR et visible
    !!! raubase = masse volumique de l'eau
    !!! rauair = masse volumique de l'air
    !!! xcp = chaleur specifique de l'eau a pression constante
    !!! tcs = coefficient de conversion de la salinite (flux de surface)
    !!! xpld = nbr de Prandtl
    !!! fave = facteur de diffision TKE
    !!! bb = parametre de calcul de TKE en surface
    !!! tgamma = filtre d'Asselin
    !!! albedo_phy = albedo ocean
    !!! emiss = emissivite ocean
    !!! stephan = coefficient de Stefan Boltzmann
    !
      real :: xlat,g,f,xdiff,xdiss,avmb,avtb,emin0,emin,ecartmx,drag, &
     & rabs,xsi1,xsi2,raubase,&
     & rauair,xcp,tcs,xpdl,fave,bb,tgamma,albedo_phy,emiss,stephan

    !
    !!! rapmax  = constante de rappel maximum en surface (en JOURS)
    !!! rapprop = coefficient sans dimension d'intensification du rappel T,S
    !!! alfcte = Coefficient sans dimension d'amplication du flux solaire
    !           moyen
    !!! alfmin = Coefficient sans dimension d'amplication du flux solaire
    !           en dessous de la moyenne
    !!! alfmax = Coefficient sans dimension d'amplication du flux solaire
    !           au dessus de la moyenne
    !
      real ::  rapmax,rapprop,alfcte,alfmin,alfmax
    !
    !        -------------------------
    ! =====> VARIABLES COMPLEMENTAIRES
    !        -------------------------
    !
    ! -------------------------FORCAGES
    !
    ! taux, tauy = tensions instantannees
    ! q, qsr = flux de chaleur totale et solaire instantannes
    ! ep = E-P instantanne
    ! fse,fle = flux sensible, latent instantannes
    ! rsw, rlw = rayonnement solaire et IR instantannes
    ! ustr, vstr = vitesse de frottement du vent
    !   indices 0 : pas de temps flux avant ; 1 : apres temps courant.
    ! njourflux = duree totale du forcage en JOURS (= INT(xperfluxt))
    ! nfreqflux = nbr de pas de temps forcages par jour (= 24./xperflux)
    ! nfluxt = nbr d'enregistrements total dans fichier forcage 
    !                       (= njourflux*nfreqflux)
    ! nflux0 = pas de temps courant en unites de periodes forcages
    ! npasflux = nbr de pas de temps run par pas de temps forcage 
    !                       (= xperflux/dt)
    ! npasflux0 = nbr de pas de temps run a partir du dernier update forcage
    !
     real ::  taux,tauy,q,qsr,ep
     real ::  fse0,fle0,rsw0,rlw0,ustr0,vstr0,fse1,fle1,rsw1,rlw1,ustr1,vstr1
     integer ::  nfluxt,nflux0,npasflux,npasflux0,njourflux,nfreqflux
    !
    ! windx, windy: vitesse du vent X et Y en m/sec
    ! wind: norme du vent en m/sec
    !
     real:: windx,windy,wind
    !
    ! qsolmean = flux solaire moyen (pour partie > 0)
    ! QSOLDAY = flux solaire journalier moyen (pour partie > 0)
    ! dayflux = jour courant en reel dans le fichier flux (milieu du pas de temps)
    !          entre 0.0 et xperfluxt
    ! dayflux0 = dayflux+xperflux0 = jour courant dans les fichiers rappels dates
    !          en fonction du fichier forcage 
    ! nday0, nday1 = jour-1 dans QSOLDAY des flux avant et apres temps courant
    !
     real ::  dayflux,dayflux0,qsolmean,QSOLDAY(jpant)
     integer::  nday0,nday1
    !
    !
    ! -------------------------RAPPELS 
    !
    !
    ! TRP, SRP = profils T et S de rappel
    ! XTS = profils de constante de rappel
    ! dayts0,dayts1 = jours initial et final inclus dans l'annee du rappel T et S
    !
     real ::  TRP(jpzt),SRP(jpzt),XTS(jpzt)
     real ::  dayts0,dayts1
    !
    ! rapfutt,rapfuts = distance courante en jours par rapport au prochain rappel
    !          en surface en T et en S
    ! rappast,rappass = distance courante en jours par rapport au dernier rappel
    !          en surface en T et en S
    ! tsurfut, ssurfut = valeurs de T et S associees a rapfut
    ! tsurpas, ssurpas = valeurs de T et S associees a rappas
    ! TSURF, SSURF = Vecteurs des valeurs de rappels de surface en T et S
    ! TTSURF,TSSURF = Vecteur des temps de rappel en T et S (en jours de l'annee)
    ! RAPTFUT,RAPSFUT = intervalle pendant lequel le rappel s'exerce apres la date
    !                    de la donnee T et S (en JOURS)
    ! rindfutt,rindfuts = 1 si jour courant < jour donnee + RAP(T-S)FUT ; = 0
    !                      autrement
    ! rindpast,rindpass = 1 si jour courant > jour donnee - RAP(T-S)PAS ; = 0
    !                      autrement
    ! ndsurft,ndsurfs = nbr de rappels TSURF - SSURF
    !
     real :: rapfutt,rappast,rapfuts,rappass,tsurfut,tsurpas,&
    &   ssurfut,ssurpas,rindfutt,rindfuts,rindpast,rindpass

     real ::  TTSURF(jpant),TSURF(jpant),RAPTFUT(jpant),&
    &RAPTPAS(jpant),TSSURF(jpant),SSURF(jpant),&
    & RAPSFUT(jpant),RAPSPAS(jpant)

     integer ::  ndsurft,ndsurfs
    !
    !
    ! -------------------------DYNAMIQUE
    !
    !
    !
    ! Suffixe B = Before ; N = Now ; A = After
    ! U, V = Vitesses zonale et meridienne
   ! E = TKE
    ! T, S = Temperature, Salinite
    ! XMDLF, XMLDS = longueurs de melange quantite de mouvement, traceur
    ! AVM, AVT, AVE = diffusion verticale quantite de mouvement, traceurs,
    !                 TKE
   ! BN = Frequence de Brunt Vaissalla
    ! SH = Cisaillement vertical
    ! RAUT, RAUW = densite aux niveaux traceur et vitesse verticale
    !
     real::  UB(jpzt),UN(jpzt),UA(jpzt),VB(jpzt),VN(jpzt),VA(jpzt), &
    & EB(jpzt),EN(jpzt),EA(jpzt), &
    & XMLDF(jpzt),XMLDS(jpzt),AVM(jpzt),AVT(jpzt),AVE(jpzt),BN(jpzt),&
    & SH(jpzt),RAUT(jpzt),&
    & RAUW(jpzt),TB(jpzt),TN(jpzt),TA(jpzt),SB(jpzt),SN(jpzt),SA(jpzt)

    !
    ! -------------------------SAUVEGARDE
    !
    !
    !
    ! TKEMOY, TMOY, SMOY, DIFFT,  DENSITE = Sommation de TKE, T, S,
    !         Diffusion traceur, densite en vue de moyenne et
    !         sauvegarde
    ! TBN, TSH, TMLDF = Sommation de BN, SH, XMLDF en vue de moyenne et
    !         sauvegarde
    ! TDFTKE = Sommation sur tendance de la diffusion dans equation TKE
    ! TDSTKE = Sommation sur tendance de la dissipation dans equation TKE
    ! TBNTKE = Sommation sur tendance de la flottabilite dans equation TKE
    ! TSHTKE = Sommation sur tendance du cisaillement dans equation TKE
    ! TRAPT, TRAPS = Sommation sur tendance rappels sur T et S
    ! TPENT = Sommation sur tendance energie penetrante sur equation de T
    ! TDIFT = Sommation sur tendance de la diffusion sur equation de T
   ! TDIFS = Sommation sur tendance de la diffusion sur equation de S
    !
     real::TKEMOY(jpzt),TMOY(jpzt),SMOY(jpzt),DIFFT(jpzt),DENSITE(jpzt)
     real::TBN(jpzt),TSH(jpzt),TMLDF(jpzt), &
    &   TDFTKE(jpzt),TDSTKE(jpzt),TBNTKE(jpzt),TSHTKE(jpzt)
     real:: TRAPT(jpzt),TRAPS(jpzt),TPENT(jpzt),TDIFT(jpzt),TDIFS(jpzt)

    !
    ! xdml = profondeur de la ML
    ! DEPTHML = Profondeur de la couche melangee
    ! FLUXSOL = Flux solaire
    ! FLUXNSOL = Flux non solaire (IR, latent, sensible)
    ! FLUXSURT = Flux rappel surface seulement temperature
    ! FLUXADVT = Flux temperature rappel colonne (advection)
    ! FLUXSAL = Flux de salinite
    ! FLUXSURS = Flux rappel surface seulement salinite
    ! FLUXADVS = Flux salinite rappel colonne (advection)
    ! FLUXVENTX = Vent X en m/sec
    ! FLUXVENTY = Vent Y en m/sec
    ! FLUXVENT = Norme du vent  en m/sec
    ! FLUXSURTKE = Forcage surface TKE
    !
     real ::  xdml,DEPTHML(jptemps),FLUXSOL(jptemps), &
    &FLUXNSOL(jptemps),FLUXADVT(jptemps),FLUXSAL(jptemps), &
    &FLUXADVS(jptemps),FLUXVENTX(jptemps),FLUXVENTY(jptemps),&
    &FLUXVENT(jptemps),FLUXSURT(jptemps),FLUXSURS(jptemps),&
    &FLUXSURTKE(jptemps)
    !
    !        ----
    ! =====> Data
    !        ----
    !
    !
    ! omega = rotation de la terre en sec-1
    !
     real:: xlim,xlimd,xlims,omega
     data xlim,xlimd,xlims /1.E-8,1.E-8,1.E-8/
     data omega /0.0000729217/
    !
    !
    ! ------>>> Double precision pour inversion diffusion
    !
    !
     real*8 DZX(jpzt),DZY(jpzt),DZZ(jpzt),DZW(jpzt),DZR(jpzt)
     real*8 DE3T(jpzt),DE3W(jpzt),DAVE(jpzt),DAVM(jpzt),DAVT(jpzt)
     real*8 dtsd,dxdiss,dfact,dtcs
     real*8 dun
     data dun /1.D+00/
     data nlengthl,ntketl,nttl,nstl /3,4,3,2/
   !------------------------------------------------------------------

    End Module comdynmod
