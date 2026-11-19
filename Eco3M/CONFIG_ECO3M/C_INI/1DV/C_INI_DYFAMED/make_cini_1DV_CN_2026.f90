!M Baklouti, Mai 2024, revised Feb 2026
! Fichier de vérification des quotas intracellulaires au sein des conditions initiales
!---------------------------------------------------------------------------------
 implicit real(Q)
 Character(30)::ficell,ficC,ficN,ficP,ficChl,ficparam
 integer::i,j,k,l,il,gf,nzmax
 real::cell,C,N,P,Chl
 real::ChlCmin,ChlCmoy,ChlCmax

 real::quotaC,quotaN,quotaP,ChlC
 real::corrC,corrN,corrP,corrChl 

 character(20)::fic_param
 character(3)::nom
 character(100)::command
 logical :: FORCE,DYFA,MO_REDF,LIV_NP_REDF

 !-- if FORCE = .TRUE. : systematically recalculates the initial conditions in C, N, P, and Chl of living organism
 !-- if FORCE = .FALSE.: verifies that the quotas for living beings are within the authorized limits and recalculates only the initial conditions that do not comply

FORCE = .TRUE.
DYFA = .TRUE.
MO_REDF=.FALSE. !Matière organique selon rapports de Redfield
LIV_NP_REDF = .FALSE. !NP selon Redfield dans le vivant

!-- lecture du fichier de paramètres:
write(*,*) "fich param = ?"
read(*,*) ficparam

!command='python3 extract_quotas_ratios('//trim(fic_param)//')'
!write(*,'(a)') command
!call system(command)

include "quotas_sim.txt"
write(*,*)"Q_MIN_PHYS_CN",Q_MIN_PHYS_CN
! initialization
ChlCmax=1.d0/(20.0/12.0) !Minimun C:Chl ratio in gC/gChl
ChlCmoy=1.d0/(70.0/12.0) 
ChlCmin=1.d0/(150/12.0)  !Max. C:Chl ratio in gC/gChl

IF (DYFA) then !DYFAMED
   nzmax = 42 ! the 42th level is only used in NEMO (level in the seafloor)
ELSE ! LEVANTIN:
   nzmax = 59 ! the 59th level is only used in NEMO (level in the seafloor)
ENDIF

i=1
j=1

Do gf = 1,9
   write(*,*)'gf=',gf
   close(15)
   close(16)
   close(17)
   close(18)
   close(19)
   close(20)
   close(21)
   close(22)
   close(23)

  if (gf==1) THEN
     ficell="zoo_COP_cell.DAT"
     ficC  ="zoo_COP_C.DAT"
     ficN  ="zoo_COP_N.DAT"
     nom=ficell(1:3)
     write(*,*) nom
  elseif (gf==2) THEN
     ficell="phy_PHY_cell.DAT"
     ficC  ="phy_PHY_C.DAT"
     ficChl="phy_PHY_Chl.DAT" 
     ficN  ="phy_PHY_N.DAT"
     nom=ficell(1:3)
     write(*,*) nom
  elseif(gf==3) THEN
     ficell="bac_BAC_cell.DAT"
     ficC  ="bac_BAC_C.DAT"
     ficN  ="bac_BAC_N.DAT"
     nom=ficell(1:3)
     write(*,*) nom
  elseif(gf==4) THEN
     ficC  ="mod_MODLC_C.DAT"
     ficN  ="mod_MODLN_N.DAT"
     ficP  ="mod_MODLP_P.DAT"
     nom=ficell(1:3)
     write(*,*) nom
  elseif(gf==5) THEN
     ficC  ="mop_MOP_C.DAT"
     ficN  ="mop_MOP_N.DAT"
     ficP  ="mop_MOP_P.DAT"
     nom=ficell(1:3)
     write(*,*) nom
  endif

! ouverture du fichier d'abondances  
  if (nom == 'zoo' .or. nom =='phy' .or. nom =='bac') then
   open (15,file=ficell)
  endif

! ouverture du fichier en carbone   
  open (16,file=ficC)
  open (17,file=ficN)
   if (nom =='phy') open (19,file=ficChl)

! read the first line (surface value) of each file
  if (nom == 'zoo' .or. nom =='phy' .or. nom =='bac') then
   read(15,*)i,j,k,cell
  endif
  read(16,*)i,j,k,C
  read(17,*)i,j,k,N
  if (nom =='phy')  read(19,*)i,j,k,Chl

  quotaC = C/cell
  quotaN = N/cell
  if (nom =='phy') ChlC = Chl/C
  
  rewind(15)
  rewind(16)
  rewind(17)
  if (nom =='phy')     rewind(19)


  corrC = -100
  corrN = -100
  corrP = -100
  corrChl = -100

  if (nom=='zoo') then
    if (quotaC>Q_MAX_MZOO_C .or. quotaC<Q_MIN_MZOO_C .or. FORCE) corrC=(Q_MIN_MZOO_C+Q_MAX_MZOO_C)/2.0
    if (quotaN > Q_MAX_MZOO_N .or. quotaN < Q_MIN_MZOO_N .OR. FORCE) corrN = (0.25*Q_MIN_MZOO_N+0.75*Q_MAX_MZOO_N)
  elseif (nom=='phy') then
    if (quotaC > Q_MAX_PHYL_C .or. quotaC < Q_MIN_PHYL_C .OR. FORCE) corrC = (0.90*Q_MIN_PHYL_C+0.10*Q_MAX_PHYL_C)
    if (quotaN > Q_MAX_PHYL_N .or. quotaN < Q_MIN_PHYL_N .OR. FORCE)  corrN = (0.25*Q_MIN_PHYL_N+0.75*Q_MAX_PHYL_N)
    if (ChlC > ChlCmax .or. ChlC < ChlCmin .OR.FORCE)  corrChl = ChlCmoy

  elseif (nom=='bac') then
    if (quotaC > Q_MAX_BAC_C .or. quotaC < Q_MIN_BAC_C .OR. FORCE)   corrC = (0.5*Q_MIN_BAC_C+0.5*Q_MAX_BAC_C)
    if (quotaN > Q_MAX_BAC_N .or. quotaN < Q_MIN_BAC_N .OR. FORCE)   corrN = (0.25*Q_MIN_BAC_N+0.75*Q_MAX_BAC_N)
  endif

  if (nom == 'zoo' .or. nom =='phy' .or. nom =='bac') then
      if (corrC>0) then
         open (20,file=trim(ficC)//"_OK")
         Do l=1,nzmax
           read(15,*)i,j,k,cell
           write(20,'(I2,2X,I2,2X,I2,2X,E13.5)')i,j,k,cell*corrC
         enddo
      endif
      rewind(15)
      if (corrN>0 .and. gf/=1) then
        open (21,file=trim(ficN)//"_OK")
        Do l=1,nzmax
          read(15,*)i,j,k,cell
          write(21,'(I2,2X,I2,2X,I2,2X,E13.5)')i,j,k,cell*corrN
        enddo
      endif
      rewind(15)
      if (corrChl>0 .and. (gf==4 .or. gf==5) ) then
        open (23,file=trim(ficChl)//"_OK")
        Do l=1,nzmax
           read(15,*)i,j,k,cell
           write(23,'(I2,2X,I2,2X,I2,2X,E13.5)')i,j,k,cell*corrC*corrChl
        Enddo
      endif
  else
          if (MO_REDF) then
           open (21,file=trim(ficN)//"_OK")
           open (22,file=trim(ficP)//"_OK")
           Do l=1,nzmax
             read(16,*)i,j,k,C
             write(21,'(I2,2X,I2,2X,I2,2X,E13.5)')i,j,k,C*16./106.
             write(22,'(I2,2X,I2,2X,I2,2X,E13.5)')i,j,k,C*1./106.
           enddo
         endif 
       rewind(16)
 endif
enddo
End

