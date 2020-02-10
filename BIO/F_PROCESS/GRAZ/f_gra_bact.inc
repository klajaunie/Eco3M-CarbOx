!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
  function f_gra_bact(ili,jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
!
!--------------------------------------------------------------------
!
! Fonction de calcul du broutage par pourcentage sur la production ou la biomasse
!
! Paramètres : les param sont ceux de la fonction qui calcule la prod bact

Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use  VAR_GLOBAL

Implicit None
!-- fonction:
 Real(8):: f_gra_bact(nx_min:nx_max,ny_min:ny_max,nz_max)
 Real(8):: bacp(nx_min:nx_max,ny_min:ny_max,nz_max)
!arguments
real(8)   :: mumax,qcmin,qnmin,qpmin
real(8)   :: var1,var2,var3,QTEN,Trem
integer    :: i,j,k,ili,jcol
!variables locales
 

!bp(1)%val(i,j,k) = f_bp_c(ili,jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)

! Calcul du broutage
do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

bacp=f_bp_c(ili,jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
f_gra_bact(i,j,k)=bacp(i,j,k)*VAR(ili)%conc(i,j,k)/nb_bact(i,j)


    enddo
  enddo
enddo

 Return
 End function f_gra_bact
