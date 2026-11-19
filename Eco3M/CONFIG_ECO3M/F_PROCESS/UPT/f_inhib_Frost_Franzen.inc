!--------------------------------------------------------------------
!Derniere modification: 12/07/2011
!
function f_inhib_Frost_Franzen(ili,jcol,inut,kn)
    !
    !--------------------------------------------------------------------
    !
    ! Cette fonction calcule le terme d'inhibition par ammonium au cours de uptake de nitrate 
    ! Frost and Franzen 1992. MEPS 83: 291-303
    !
    ! variables et unites:
    !---------------------
    ! kn : constante d'inhibition de l'uptake de nitrate ou DON par ammonium, ou de DOP par PO4, molX m-3
    !
    ! Romain Mauriac
    !--------------------------------------------------------------------------------------------

   ! Variables globales

   Implicit None

  !-- fonction
  real(8) :: f_inhib_Frost_Franzen(nx_min:nx_max,ny_min:ny_max,1:nz_max)

  !-- arguments
  Real(8) :: kn,inut
  Integer :: ili,jcol,id 

  !-- variables locales:
  character(L_CHAIN)::chain1,chain2
  integer :: i,j,k

  id=nint(inut)

  
  f_inhib_Frost_Franzen = 1.d0

!  do k = 1,nz_max
!    do j = ny_min,ny_max
!       do i = nx_min,nx_max
!       if (var(id)%conc(i,j,k)>1d-10) f_inhib_Frost_Franzen(i,j,k)=1.d0/(1.d0+(var(id)%conc(i,j,k)/(kn+1d-50)))
!     enddo
!   enddo
!  enddo

  !do i = nx_min,nx_max
  !  do j = ny_min,ny_max
  !    do k = 1,nz_max
  !IF (ISNAN(f_inhib_Frost_Franzen(i,j,k))) THEN
  !WRITE(*,*) 'FUNCTION f_inhib_Frost_Franzen contains NaN',i,j,k
  !STOP
  !ENDIF
  !enddo
  !enddo
  !enddo

end function f_inhib_Frost_Franzen
