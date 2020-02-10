!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
    function f_upt_droop(ili,jcol,Qmax,Qmin,coeff_up,up_part,ks)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule l''absorption d''un nutriment
! par le modele de Droop 
! 
!! Paramètres :
!                      - upmax : uptake max en molN/molC/s
!                      - Qmax et Qmin : rapport max et min N/C
!                      - coeff_up : coefficient de partage entre plusieurs uptake
!                      - up_part : coefficient eventuel pour pondérer l'uptake
!                      - ks : constante de demi-saturation pour l'uptake
!
! Variables globales
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use DEF_TYPE

 Implicit None
!-- fonction
real(8)  :: f_upt_droop(nx_min:nx_max,ny_min:ny_max,nz_max)
!-- arguments
integer  :: i,j,k,ili,jcol,iorg,iphy
real(8)  :: Qmax,Qmin,coeff_up,up_part,ks
real(8)  :: Q(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8)  :: hQ(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8)  :: upmax(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8)  :: up_N(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8)  :: var_C(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8)  :: f_upt_at(nx_min:nx_max,ny_min:ny_max,nz_max)

do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max

iorg = VAR(ili)%idorg

if (VAR(ili)%scomp.ne.'at')then

!-- Calcul des variables intermédiaires
if(VAR(jcol)%scomp.eq.'po4') then
    var_C(i,j,k)=VAR(ili-2)%conc(i,j,k)
else
    var_C(i,j,k)=VAR(ili-1)%conc(i,j,k)
endif
 Q(i,j,k)=VAR(ili)%conc(i,j,k)/var_C(i,j,k) ! rapport N/C ou P/C  

!-- Uptake du nutriments VAR(jcol) suivant LEHMAN (1975)
! empeche valeurs negatives Katixa Lajaunie-Salla (19/07/2018)
if(Q(i,j,k)>Qmax) then
   hQ=0.
elseif (Q(i,j,k)<Qmin) then
   hQ(i,j,k)=1.
else
   hQ(i,j,k)=((Qmax-Q(i,j,k))/(Qmax-Qmin))
endif

!!-- Calcul de la vitesse d'uptake (molX/molC/s) --!!
upmax(i,j,k)=Qmax*mu_phy(iorg)%val(i,j,k) !PCnr s-1 Pmax*f(T)*f(I)
!Uptake pondérée en cas d'uptake sur l'Azote provenant de plusieurs sources (NOD,NH4,NO3...) (coef_up=/1)
up_N=up_part*coeff_up*upmax(i,j,k)*VAR(jcol)%conc(i,j,k)/(ks+VAR(jcol)%conc(i,j,k))  ! VAR(jcol) = Nutriments (up_N en molX/molC/s)
if(VAR(jcol)%scomp.eq.'no3'.and.VAR(ili)%scomp.eq.'pico') then ! inhibition 
  up_N(i,j,k)=up_part*coeff_up*upmax(i,j,k)*(VAR(jcol)%conc(i,j,k)/(ks+VAR(jcol)%conc(i,j,k)))*(1.-0.82*VAR(jcol-1)%conc(i,j,k)/(0.578+VAR(jcol-1)%conc(i,j,k)))
endif

!!-- Uptake final en molX/m3/s --!!
f_upt_droop(i,j,k)=up_N(i,j,k)*var_C(i,j,k) ! on multiplie par le carbone

!!-- Enregistrement de la valeur de l'uptake pour l'exsudation et pour la production d'oxygène --!!
if(VAR(jcol)%scomp.eq.'nh4') then
  upnh(iorg)%val(i,j,k)=up_N(i,j,k)
  hQnh(iorg)%val(i,j,k)=hQ(i,j,k)
  fupnh(iorg)%val(i,j,k)=f_upt_droop(i,j,k)
  Qnh(iorg)%val(i,j,k)=Q(i,j,k)
elseif(VAR(jcol)%scomp.eq.'no3') then
  upno(iorg)%val(i,j,k)=up_N(i,j,k)
  hQno(iorg)%val(i,j,k)=f_upt_droop(i,j,k)
  fupno(iorg)%val(i,j,k)=f_upt_droop(i,j,k)
elseif(VAR(jcol)%scomp.eq.'po4') then
  uppo(iorg)%val(i,j,k)=up_N(i,j,k)
  hQpo(iorg)%val(i,j,k)=hQ(i,j,k)
  fuppo(iorg)%val(i,j,k)=f_upt_droop(i,j,k)
  Qpo(iorg)%val(i,j,k)=Q(i,j,k)
endif


else

f_upt_at(i,j,k)=0.
do iphy=1,nscp_phy
  f_upt_at(i,j,k)=f_upt_at(i,j,k)+(-fupnh(iphy)%val(i,j,k)+fupno(iphy)%val(i,j,k)+fuppo(iphy)%val(i,j,k))
enddo
f_upt_droop(i,j,k)=f_upt_at(i,j,k)
endif

enddo
enddo
enddo

Return
End function f_upt_droop
