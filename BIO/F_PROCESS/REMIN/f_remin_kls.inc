!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
 function f_remin(ili,jcol,upmax1,upmax2,ks1,ks2,QTEN,Trem,qmin,qmax)
!
!--------------------------------------------------------------------
! Modifié par Katixa Lajaunie-Salla 07/06/2018
! Cette fonction calcule la vitesse de remineralisation de la matière organique dissoute et particulaire par les bactéries
! Degradation de la matière N et P produit de l'ammonium et du phosphate 
! la partie mineralizee est (1-hQba) de tout les uptake par le phyto
!
! variables et unites:
!---------------------
! upmax1 : uptake max de la MOD (s^{-1})
! upmax2 : uptake max de la MOP (s^{-1})
! ks1 : constante de demi-saturation de Michaelis-Menten pour la MO (µM)
! ks2 : constante de demi-saturation de Michaelis-Menten pour la MI (µM)
! QTEN: dans le calcul de la fonction de temperature
! Trem: temperature optimale (°C)
! Qmin : quota min BacN/nbBac ou BacP/nbBac
! Qmax : quota max BacN/nbBac ou BacP/nbBac
! f_remin : vitesse de remineralisation, en s^{-1}
!
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use DEF_TYPE

Implicit None
!-- fonction
Real(8) :: f_remin(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: f_remin_at(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: uppart(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: updiss(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: upmineral(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: Q_ba(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: hQ_ba(nx_min:nx_max,ny_min:ny_max,nz_max)
!-- arguments
Real(8) :: upmax1,upmax2,ks1,ks2,qmax,qmin,QTEN,Trem
Integer :: i,j,k,ili,ili_col,jcol,jcol_tempo

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

if (VAR(ili)%scomp.ne.'at')then

!-- Calcul des variables intermédiaires
Q_ba(i,j,k)=VAR(ili)%conc(i,j,k)/nb_bact(i,j)
! rapport Bc/B ou Bn/B ou Bp/B
! empeche valeurs negatives Katixa Lajaunie-Salla (19/07/2018)
if(Q_ba(i,j,k)>qmax) then
   hQ_ba(i,j,k)=0.
elseif (Q_ba(i,j,k)<qmin) then
   hQ_ba(i,j,k)=1.
else
   hQ_ba(i,j,k)=(qmax-Q_ba(i,j,k))/(qmax-qmin)
endif

jcol_tempo=ili+3 !POM
uppart=f_upt_bact_droop(ili,jcol_tempo,upmax2,ks1,QTEN,Trem) !!(en molX/m3/s)
jcol_tempo=ili+6 !DOM
updiss=f_upt_bact_droop(ili,jcol_tempo,upmax1,ks1,QTEN,Trem) !!(en molX/m3/s)
jcol_tempo=jcol !DIM
upmineral=f_upt_bact_droop(ili,jcol_tempo,upmax1,ks2,QTEN,Trem) !!(en molX/m3/s)

f_remin(i,j,k)=(1-hQ_ba(i,j,k))*(uppart(i,j,k)+updiss(i,j,k)+upmineral(i,j,k))

!!-- Enregistrement de la valeur de remineralisation pour l'alcalinite --!!
if(VAR(jcol)%scomp.eq.'nh4') then
  remnh(1)%val(i,j,k)=f_remin(i,j,k)
elseif(VAR(jcol)%scomp.eq.'po4') then
  rempo(1)%val(i,j,k)=f_remin(i,j,k)
endif

else
f_remin_at(i,j,k)=remnh(1)%val(i,j,k)
f_remin(i,j,k)=f_remin_at(i,j,k)
endif

      enddo
  enddo
enddo

Return

end function f_remin

