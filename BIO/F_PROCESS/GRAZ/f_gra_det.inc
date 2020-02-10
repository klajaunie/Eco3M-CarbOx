!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
              function f_gra_det(ili,jcol,g1,g2)
!
!--------------------------------------------------------------------
!
! Fonction de calcul du broutage des detritus par le zoo
!
! Paramètres
!
!  - g : pourcentage de broutage sur la biomasse phytoplanctonique
!

Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL

Implicit None
!-- fonction:
Real(8) :: f_gra_det(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: limG(nx_min:nx_max,ny_min:ny_max,nz_max)
! Variables internes
integer :: i,j,k,ili,jcol
real(8) :: g,g1,g2
real(8) :: phyto1,phyto2,tempo

! Calcul du broutage
do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

!limitation par la temperature
limG=f_lim_graz(ili)

tempo = 0.

if(VAR(ili)%elmt=="C") then
  phyto1 = VAR(BIO_C(1))%conc(i,j,k)
  if(VAR(BIO_C(2))%comp=='phy') then
    phyto2 = VAR(BIO_C(2))%conc(i,j,k)
  else
    phyto2 = 0.
  endif
 if((phyto1+phyto2)/=0.) then
   tempo = (g1*phyto1+g2*phyto2)/(phyto1+phyto2)
 else
   tempo =0.
 endif

else if (VAR(ili)%elmt=="N") then
 phyto1 = VAR(BIO_N(1))%conc(i,j,k)
  !debug : verif deux phytos
  if(VAR(BIO_N(2))%comp=='phy') then
    phyto2 = VAR(BIO_N(2))%conc(i,j,k)
  else
    phyto2 = 0.
  endif
 if((phyto1+phyto2)/=0.) then
  tempo = (g1*phyto1+g2*phyto2)/(phyto1+phyto2)
 else
  tempo =0.
 endif
 
else if (VAR(ili)%elmt=="P") then
 phyto1 = VAR(BIO_P(1))%conc(i,j,k)
  if(VAR(BIO_P(2))%comp=='phy') then
    phyto2 = VAR(BIO_P(2))%conc(i,j,k)
  else
    phyto2 = 0.
 endif
 if((phyto1+phyto2)/=0.) then
  tempo = (g1*phyto1+g2*phyto2)/(phyto1+phyto2)
 else
  tempo =0.
 endif

endif

f_gra_det(i,j,k) = tempo*limG(i,j,k)*VAR(ili)%conc(i,j,k)

    enddo
  enddo
enddo

Return

End function f_gra_det
