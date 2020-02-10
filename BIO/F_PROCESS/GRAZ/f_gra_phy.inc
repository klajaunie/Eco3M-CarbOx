!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
              function f_gra_phy(ili,jcol,g1,g2)
!
!--------------------------------------------------------------------
!
! Fonction de calcul du broutage par pourcentage sur la production ou la biomasse
!
! modifier par Katixa lajaunie-Salla 25/04/2018
!
! Paramètres
!  -  g : pourcentage de broutage sur la biomasse
!
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use  VAR_GLOBAL

Implicit None
!-- fonction
Real(8) :: f_gra_phy(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: limG(nx_min:nx_max,ny_min:ny_max,nz_max)
! Variables internes
integer  :: i,j,k,ili,jcol,iorg
real(8)  :: g1,g2,g

iorg = VAR(ili)%idorg

! Calcul du broutage
do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max

!limitation par la temperature
limG=f_lim_graz(ili)

if(VAR(ili)%scomp=="pico") then
  g=g1 
elseif (VAR(ili)%scomp=="dia") then 
  g=g2
endif

f_gra_phy(i,j,k)=g*limG(i,j,k)*VAR(ili)%conc(i,j,k)

enddo
enddo
enddo

 Return
 End function f_gra_phy
