!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
              function f_out(ili,jcol)
!
!--------------------------------------------------------------------
!
! Modification par Katixa Lajaunie-Salla 07/05/2018 
! Prise en compte de la fonction de quota de croissance et d'uptake dans le calcul de la partie exsudée
!
! Variables globales
USE VAR_GLOBAL
USE VAR_USER
USE COUPLEUR_PHY_BIO

 Implicit None
!-- fonction
Real(8) :: f_out(nx_min:nx_max,ny_min:ny_max,nz_max)
integer :: i,j,k,ili,jcol,iorg
real(8) :: qmin,qmax,k_exsu,f_exsu,var_C
 
do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max

iorg = VAR(ili)%idorg 

if (VAR(ili)%elmt=="C") then
  k_exsu=mu_phy(iorg)%val(i,j,k) !PCnr s-1 = Pmax*f(T)*f(I)
  f_exsu=lim_nut(iorg)%val(i,j,k)
  var_C=VAR(ili)%conc(i,j,k)
elseif (VAR(ili)%elmt=="N") then
  k_exsu=upnh(iorg)%val(i,j,k)+upno(iorg)%val(i,j,k) ! molX/molC/s
  f_exsu=hQnh(iorg)%val(i,j,k)
  var_C=VAR(ili-1)%conc(i,j,k)
elseif (VAR(ili)%elmt=="P") then
  k_exsu=uppo(iorg)%val(i,j,k)
  f_exsu=hQpo(iorg)%val(i,j,k)
  var_C=VAR(ili-2)%conc(i,j,k)
endif
  
f_out(i,j,k)=(1-f_exsu)*k_exsu*var_C 

enddo
enddo
enddo
 Return
 End function f_out
