!************************************************************************************************************
!***** CALCUL DE LA LUMIERE POUR LA PPB
!************************************************************************************************************
!--------------------------------------------------------------------
!
               function f_lum(kk0,aa0,Chl_tot,dz)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule la limitation en lumière
!


! Variables globales
!Use  DEF_TYPE
!Use  COUPLEUR_PHY_BIO
!Use  VAR_GLOBAL

 Implicit None
!-- fonction
!Real(8) :: f_lum(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: f_lum

!-- arguments
 real(8) :: kk0  !(0.0384)  coefficient d''extinction de l''eau de mer pure (m-1)
 real(8) :: aa0  !(0.01)  coefficient d''absorption de  Chla (m2*mg-1)
 real(8),optional :: Chl_tot,dz


f_lum = exp(-dz*(kk0+aa0*Chl_tot)) 

 Return
 End function f_lum
