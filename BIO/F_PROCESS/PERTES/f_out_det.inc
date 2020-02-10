!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
  function f_out_det(ili,jcol,D,g1,g2,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
!
!--------------------------------------------------------------------
!
! cette fonction calcule les rejets sous forme de pelotes fécales
! le zoo rejette la partie non assimilée (1-D).grtot sous forme MOP
!
! Modification par Katixa Lajaunie-Salla 01/06/2018 
! Avec deux phytoplanctons
!
! Paramètres :       
!  - D :  coefficient d'assimilation
!  - g1,g2 : grazing du phyto
!
! Variables globales
Use DEF_TYPE
USE VAR_GLOBAL
USE VAR_USER
USE COUPLEUR_PHY_BIO

 Implicit None
!-- fonction
Real(8) :: f_out_det(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grphy(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grbac(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grdet(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grtot(nx_min:nx_max,ny_min:ny_max,nz_max)
!--Variables locales:
integer    :: i,j,k,ili,jcol,ivar
real(8)    :: D,g1,g2
real(8)    :: qcmin,qnmin,qpmin,mumax
real(8)    :: var1,var2,var3,QTEN,Trem

do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max
	
!calcul du broutage totale (carbone, azote ou phosphore)

grphy(i,j,k)=0.

if(VAR(ili)%elmt=="C".OR.VAR(ili)%elmt=="O2") then
    do ivar=1,size(BIO_C)
     if(VAR(BIO_C(ivar))%comp=="phy") then
        grphy=grphy+f_gra_phy(BIO_C(ivar),jcol,g1,g2)
     elseif(VAR(BIO_C(ivar))%comp=="bac") then
         grbac=f_gra_bact(BIO_C(ivar),jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
     elseif(VAR(BIO_C(ivar))%comp=="det") then
         grdet=f_gra_det(BIO_C(ivar),jcol,g1,g2)
     endif
    enddo
grtot(i,j,k)=grphy(i,j,k)+grbac(i,j,k)+grdet(i,j,k)

else if (VAR(ili)%elmt=="N".OR.VAR(ili)%elmt=="nh4") then
    do ivar=1,size(BIO_N)
     if(VAR(BIO_N(ivar))%comp=="phy") then
        grphy=grphy+f_gra_phy(BIO_N(ivar),jcol,g1,g2)
     elseif(VAR(BIO_N(ivar))%comp=="bac") then
         grbac=f_gra_bact(BIO_N(ivar),jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
     elseif(VAR(BIO_N(ivar))%comp=="det") then
         grdet=f_gra_det(BIO_N(ivar),jcol,g1,g2)
     endif 
    enddo
grtot(i,j,k)=grphy(i,j,k)+grbac(i,j,k)+grdet(i,j,k)

 else if (VAR(ili)%elmt=="P".OR.VAR(ili)%elmt=="po4") then
    do ivar=1,size(BIO_P)
     if(VAR(BIO_P(ivar))%comp=="phy") then
        grphy=grphy+f_gra_phy(BIO_P(ivar),jcol,g1,g2)
     elseif(VAR(BIO_P(ivar))%comp=="bac") then
         grbac=f_gra_bact(BIO_P(ivar),jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
     elseif(VAR(BIO_P(ivar))%comp=="det") then
         grdet=f_gra_det(BIO_P(ivar),jcol,g1,g2)
     endif    
    enddo
grtot(i,j,k)=grphy(i,j,k)+grbac(i,j,k)+grdet(i,j,k)

endif

  f_out_det(i,j,k)=(1.-D)*grtot(i,j,k)

    enddo
  enddo
enddo
 Return
 End function f_out_det
