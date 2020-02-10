!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
   function f_nit(ili,jcol,nitmax0,QTEN,Tnit,KO2)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule la vitesse de remineralisation
! d''un nutriment
!
! variables et unites:
!---------------------
!
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use DEF_TYPE
!USE MOD_FCHAIN
Implicit None
!-- fonction
Real(8) :: f_nit(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: lim_t(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: limO2(nx_min:nx_max,ny_min:ny_max,nz_max)
!-- arguments
Real(8) :: nitmax0,QTEN,KO2,Tnit
Integer :: i,j,k,ili,jcol
character(L_CHAIN) :: scomp,elmt
!-- variables locales:
Integer :: id_nh4,id_o2

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

! Calcul de la fonction de temperature
lim_t(i,j,k)=QTEN**((TEMP_BIO(i,j,k)-Tnit)/10)
! identification de la variable ammonium
scomp = "nh4"
elmt  = "nh4"
id_nh4 = f_scomp2id(scomp,elmt)

! Calcul de la fonction Michaelis-Menten pour l'oxygene
scomp="O2"
elmt="O2"
id_o2=f_scomp2id(scomp,elmt)
limO2(i,j,k)=VAR(id_o2)%conc(i,j,k)/(KO2+VAR(id_o2)%conc(i,j,k))

f_nit(i,j,k) = VAR(id_nh4)%conc(i,j,k)*nitmax0*lim_t(i,j,k)*limO2(i,j,k)


! Remarque importante: vu qu'en rentrant les variables dans "config.ini",
! on commence toujours par la biomasse, le nutriment absorbe
! est toujours dans la colonne jcol du flux (ili --> jcol), et la
! biomasse est toujours dans la ligne ili

    enddo
  enddo
enddo

Return
end function f_nit
