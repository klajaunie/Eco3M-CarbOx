!************************************************************************************************************
! --------- Modèles de Production bactérienne------------

!--------------------------------------------------------------------
!
  function f_bp_c(ili,jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
!
!--------------------------------------------------------------------
!
! Cette fonction permet de calculer la production bactérienne
!  sous la forme d'un modèle de Droop : Qmin est le ratio N/C minimal
!  permettant de soutenir la croissance
! on calcule une production nette, d'où le calcul de la respiration bactérienne

! Paramètres
!  - mumaxba : le taux de croissance maximal des bactéries en s-1
!  - qmin : le quota cellulaire N/C minimal pour les bactéries.

! Variables globales
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL

Implicit None
!-- fonction
Real(8)  :: f_bp_c(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8)  :: mu_ba(nx_min:nx_max,ny_min:ny_max,nz_max)
!Real(8)  :: limTB(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8)  :: Q1(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8)  :: Q2(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8)  :: Q3(nx_min:nx_max,ny_min:ny_max,nz_max)
integer  :: i,j,k,ili,jcol
real(8)  :: mumax,qcmin,qnmin,qpmin,var1,var2,var3,Trem,QTEN

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

  Q1(i,j,k) = VAR(var1)%conc(i,j,k)/nb_bact(i,j)
  Q2(i,j,k) = VAR(var2)%conc(i,j,k)/nb_bact(i,j)
  Q3(i,j,k) = VAR(var3)%conc(i,j,k)/nb_bact(i,j)

  limbaC(1)%val(i,j,k)=Qcmin/Q1(i,j,k)
  limbaN(1)%val(i,j,k)=Qnmin/Q2(i,j,k)
  limbaP(1)%val(i,j,k)=Qpmin/Q3(i,j,k)

!-- Droop formulation 1-Qmin/Q --!
 lim_ba(1)%val(i,j,k)=max(limbaC(1)%val(i,j,k),limbaN(1)%val(i,j,k),limbaP(1)%val(i,j,k))

!- Taux de croissance
 if(lim_ba(1)%val(i,j,k).lt.1) then
   mu_ba(i,j,k) = mumax * (1-lim_ba(1)%val(i,j,k))
 else
   mu_ba(i,j,k) = 0.
 endif

!- Production bactérienne
 bp(1)%val(i,j,k)=mu_ba(i,j,k)*nb_bact(i,j)*QTEN**((TEMP_BIO(i,j,k)-Trem)/10.)

 f_bp_c(i,j,k)=bp(1)%val(i,j,k)


    enddo
  enddo
enddo

 Return
 End function f_bp_c
