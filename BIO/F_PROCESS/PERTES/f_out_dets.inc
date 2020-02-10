!--------------------------------------------------------------------
!
 function f_out_dets(ili,jcol,F,g1,g2,choix)
!
!--------------------------------------------------------------------
!
!cette fonction calcule l'excretion en MOP du zooplancton de facon simple : fraction du grazing sur le C_Phy
! Elle calcule aussi l'excretion en NH du zoo...
! 
!
! Paramètres :       
!               F :  Fraction du grazing phyto rejeté sous forme de feces, ou de NH4
!               g : grazing du phyto
!
!
!


! Variables globales
USE VAR_GLOBAL
USE VAR_USER
USE COUPLEUR_PHY_BIO

 Implicit None
!-- fonction
Real(8) :: f_out_dets(nx_min:nx_max,ny_min:ny_max,nz_max)


!--Variables locales:
integer    :: i,j,k,ili,jcol,ivar
real(8)    :: F ,g,choix,g1,g2
real(8)   :: jour_s

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max


! Determination du broutage, en fonction de l'eutrophisation
!if(VAR(1)%conc(i,j,k).GT.2.0) then
!  g = g1
!else
  g = g2
!endif

!jour_s = 1.1574E-5
jour_s=1.

  f_out_dets(i,j,k) = F * jour_s * g * VAR(choix)%conc(i,j,k)

      enddo
  enddo
enddo

 Return
 End function f_out_dets
