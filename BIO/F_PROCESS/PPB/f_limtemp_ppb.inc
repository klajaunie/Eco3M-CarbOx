!************************************************************************************************************
! --------- Modèles de limiation du grazing par la temperature------------

!--------------------------------------------------------------------
!
              function f_limtemp_ppb(ili,tref,tlet,b)
!
!--------------------------------------------------------------------
!
! Cette fonction permet de calculer la limitation de la production primaire par la temperature
!  
! Fichier crée 04/06/2018 Katixa Lajaunie-Salla

! Paramètres
!  - temp=temperature reference
!  - tlet
!  - b
!  
! Variables globales

Use  DEF_TYPE
Use COUPLEUR_PHY_BIO
USE  VAR_GLOBAL
Implicit None

!-- fonction
Real(8) :: f_limtemp_ppb(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: limtemp(nx_min:nx_max,ny_min:ny_max,nz_max)
!Real(8) :: temp(nx_min:nx_max,ny_min:ny_max,nz_max)
integer :: i,j,k,ili,jcol,iorg
real(8) :: tref,tlet,b,ftemp

iorg = VAR(ili)%idorg

do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max


limtemp(i,j,k)=(2.*(1.-b)*(TEMP_BIO(i,j,k)-tlet)/(tref-tlet))/(((TEMP_BIO(i,j,k)-tlet)/(tref-tlet))**2+2*-b*(TEMP_BIO(i,j,k)-tlet)/(tref-tlet)+1.)

f_limtemp_ppb(i,j,k)=max(limtemp(i,j,k),0.d0)

limT_ppb(iorg)%val(i,j,k)=f_limtemp_ppb(i,j,k)

enddo
enddo
enddo
 
 Return
 End function f_limtemp_ppb
