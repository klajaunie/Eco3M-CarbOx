 function f_ppb_geider(ili,jcol,mumaxa,tref,tlet,b,qmin,qmax,qpmin,qpmax,qminchl,qmaxchl,alpha)
!--------------------------------------------------------------------
!
! Cette fonction calcule la photosynthèse, et la croissance du phytoplancton
!
! Rq importante : on considère que pour chaque phyto, l'élément qui suit le C est le N
!
! Variables globales
Use  DEF_TYPE
Use COUPLEUR_PHY_BIO
USE  VAR_GLOBAL
 Implicit None
!-- fonction
 Real(8) :: f_ppb_geider(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8)  :: limT(nx_min:nx_max,ny_min:ny_max,nz_max)
 integer :: i,j,k,ili,jcol,iorg,iichl,iin,iip
 Real(8) :: mumaxa,tref,tlet,b,ks,qmin,qmax,qpmin,qpmax,alpha,qmaxchl,qminchl
 Real(8) :: Qn,Qp
 Real(8) :: mumaxt,mumaxC,mu ! Limitation par la temperature, lumière
 Real(8) :: Chl,theta,tempo,N,P

Character(L_VAR)::elmt

iorg = VAR(ili)%idorg

do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max

!! -- Limitation par la temperature (ajout Katixa Lajaunie Salla 15/05/2018) -- !!
limT=f_limtemp_ppb(ili,tref,tlet,b)

!! -- Limitation par les nutriments ("Droop") -- !!
!- ATTENTION : l'élément qui suit le C est le N puis le P
Qn=VAR(ili+1)%conc(i,j,k)/VAR(ili)%conc(i,j,k) ! rapport N/C
Qp=VAR(ili+2)%conc(i,j,k)/VAR(ili)%conc(i,j,k) ! rapport P/C

! Empeche les valeurs negatives Katixa Lajaunie Salla (19/07/2018)
if (VAR(ili)%scomp.eq.'pico') then
limn(iorg)%val(i,j,k)=(Qn-qmin)/(Qn-qmin+0.0072)
limp(iorg)%val(i,j,k)=(Qp-qpmin)/(Qp-qpmin+0.0002)
elseif (VAR(ili)%scomp.eq.'dia') then
limn(iorg)%val(i,j,k)=(Qn-qmin)/(Qn-qmin+0.002)
limp(iorg)%val(i,j,k)=(Qp-qpmin)/(Qp-qpmin+0.0005)
endif

if (Qn>qmax) then
  limn(iorg)%val(i,j,k)=1.
elseif (Qn<qmin) then
  limn(iorg)%val(i,j,k)=0.
endif

if (Qp>qpmax) then
   limp(iorg)%val(i,j,k)=1.
elseif (Qp<qpmin) then
   limp(iorg)%val(i,j,k)=0.
endif

!determination du terme de limitation le plus fort entre le P et le N
lim_nut(iorg)%val(i,j,k)=min(limn(iorg)%val(i,j,k),limp(iorg)%val(i,j,k))
lim_nut(iorg)%val(i,j,k)=max(lim_nut(iorg)%val(i,j,k),1e-05)
!
!Carbon-specific, light-saturated rate of photosynthesis (PCmax dans Geider 1998,pour lumiere max : Pmax*f(T)*f(Nut))
mumaxC=mumaxa*lim_nut(iorg)%val(i,j,k)*limT(i,j,k) ! en s-1

!! -- Calcul de theta le rapport Chl:C -- !!
!-- cas d'un rapport CHL:C constant:
if (.not. CHL_C_BOOL) then
  Chl = CHL_C0 * var(ili)%conc(i,j,k)
  theta = CHL_C0 !-- Rapport Chl:C:
else
  elmt ='Chl' !Indice correspondant a la conc. chlorophylienne du phyto
  iichl = f_idorg2id (iorg,elmt)
  elmt ='N' !Indice correspondant a la conc. en azote du phyto
  iin = f_idorg2id (iorg,elmt)
  N = var(iin)%conc(i,j,k)
  elmt ='P' !Indice correspondant a la conc. en phosphore du phyto
  iip = f_idorg2id (iorg,elmt)
  P = var(iip)%conc(i,j,k)
  Chl = var(iichl)%conc(i,j,k)
  theta = Chl/(var(ili)%conc(i,j,k)+1.d-50)
! ATTENTION AJOUT DU CALCUL DE LA CHLORO PAR VAN DER MEERSCHE: Q*( coorespond a lim_nut)  varie en foction de P
  var(iichl)%conc(i,j,k) = VAR(iin)%conc(i,j,k)/VAR(ili)%conc(i,j,k)&
    *(qmaxchl*lim_nut(iorg)%val(i,j,k)+qminchl*(1.-lim_nut(iorg)%val(i,j,k)))*VAR(ili)%conc(i,j,k)
endif

!! -- Limitation par la lumiere -- !!
lim_lum(iorg)%val(i,j,k)=1.-exp(-alpha*theta*E_PARZ(i,j,k)/(mumaxC+1.d-50))

!! -- Growth Primary Production or Carbon-specific nutrient replete rate of phtosynthesis s-1 (PCnr dans Baklouti 2006 : Pmax*f(T)*f(I))
mu_phy(iorg)%val(i,j,k)=mumaxa*limT(i,j,k)*lim_lum(iorg)%val(i,j,k) ! en s-1

!! -- Croissance molC/m3/s -- !!
f_ppb_geider(i,j,k)=mu_phy(iorg)%val(i,j,k)*VAR(ili)%conc(i,j,k)
if(f_ppb_geider(i,j,k).LT.0.)then
f_ppb_geider(i,j,k)=0.
endif


!! -- Stockage de la production primaire -- !!
 PPB(iorg)%val(i,j,k) = f_ppb_geider(i,j,k)

enddo
enddo
enddo

Return
End function f_ppb_geider
