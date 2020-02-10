!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
    function f_upt_bact_droop(ili,jcol,upmax,ks,QTEN,Trem)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule l'uptake des bactéries sur l'Azote et Phosphore: NOP,NOD,NH4,POP,POD,PO4
! à partir du modèle de Droop
!
! Paramètres :
!       
!  - upmax : uptake max en moC/cell/s ou molN/cell/s
!  - Qmax et Qmin : rapport max et min C/cell ou N/cell
!  - coeff_up : coefficient de partage entre uptake bactérien NH4 et NOD
!  - up_part : pourcentage de bactéries consommant de la MOP (ici du NOP)
!  - ks : constante de demi-saturation pour l'uptake d'ammonium (== uptake NOD)i ou COD

! Variables globales
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use DEF_TYPE

Implicit None
!-- fonction
Real(8)  :: f_upt_bact_droop(nx_min:nx_max,ny_min:ny_max,nz_max)
!-- arguments
integer  :: i,j,k,ili,jcol
Real(8)  :: Q_ba(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8)  :: hQ_ba(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8)  :: up(nx_min:nx_max,ny_min:ny_max,nz_max)
!Real(8)  :: limTB(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8)  :: upmax,ks,Trem,QTEN

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

!-- Uptake suivant DROOP (LEHMAN 1975) --!! (up en molX/cell/s)
up(i,j,k)=upmax*QTEN**((TEMP_BIO(i,j,k)-Trem)/10.)*VAR(jcol)%conc(i,j,k)/(ks+VAR(jcol)%conc(i,j,k))  ! VAR(jcol) = Nutriments

!-- Uptake final (en molX/m3/s) --!!
f_upt_bact_droop(i,j,k)=up(i,j,k)*nb_bact(i,j) ! on multiplie par la biomasse pour revenir en C ou N

    enddo
  enddo
enddo

 Return
 End function  f_upt_bact_droop
