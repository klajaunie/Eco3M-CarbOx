!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
  function f_out_mod(ili,jcol,D,R,KZO,g1,g2,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
!
!--------------------------------------------------------------------
!
!cette fonction calcule l'excrétion en MOD et NH4/PO4 du zooplancton
!le zoo rejette la partie assimilée (D) et celle qui n'est pas utilisée pour la croissance (1-K) sous forme MOD et MID
!
! Modification par Katixa Lajaunie-Salla 01/06/2018 
! Avec deux phytoplanctons
!
! Paramètres :
!  - R : Respiration du zoo (pour COD) ou Exudation NH4 du zoo (pour NOD)
! ATTENTION : pour l'excretion en NH4, il faut donc rentrer pour R : R= 1-exu
!  - D     : coefficient d'assimilation
!  - K     : coefficient d'utilisation de la nourriture assimilée
!  - g1,g2 : grazing du phyto
!
!! Variables globales
USE VAR_GLOBAL
USE VAR_USER
USE COUPLEUR_PHY_BIO
Use  DEF_TYPE

 Implicit None
!-- fonction
Real(8) :: f_out_mod(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grphy(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grbac(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grdet(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: grtot(nx_min:nx_max,ny_min:ny_max,nz_max)
!--Variables locales:
integer    :: i,j,k,ili,jcol,ivar
real(8)    :: D,R,KZO,g1,g2
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

! Parti non refractaire de la MOD 
! tempo*biomasse = ingestion totale du zoo
if(VAR(ili)%comp.EQ."mod") then
  f_out_mod(i,j,k)=(1.-R)*D*(1.-KZO)*grtot(i,j,k)
else
  f_out_mod(i,j,k)=R*D*(1.-KZO)*grtot(i,j,k)
endif


enddo
enddo
enddo

 Return
 End function f_out_mod
