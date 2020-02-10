 function f_ppb_geider(ili,jcol,mumaxa,mumaxb,qmin,qmax,qpmin,qpmax,qmaxchl,qminchl,alpha)
!--------------------------------------------------------------------
!
! Cette fonction calcule la photosynthèse, et la croissance du phytoplancton
!
! Rq importante : on considère que pour chaque phyto, l'élément qui suit le C est le N
!
!AJOUT P marion 22/04/10
!Marion 11/10/2012 :
!- suppression cout biosynthese,
!-declaration quota pour chloro Van Der Meersche dans fonction et config.ini car avant en dur
!-suppression de +0.5 dans formualtion chloro Van Der Meersche


! Variables globales
Use  DEF_TYPE
Use COUPLEUR_PHY_BIO
USE  VAR_GLOBAL
 Implicit None
!-- fonction
Real(8) :: f_ppb_geider(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
integer    :: i,j,k,ili,jcol
real(8)       :: mumaxa,mumaxb,ks,qmin,qmax,qpmin,qpmax,alpha,qmaxchl,qminchl
!--Variables locales:
real(8)    :: mumaxt         ! Limitation par la lumière,par les nutriments
real(8)    :: Chl, theta,tempo,N
!real(8)    :: f_limi
integer    :: iorg,iichl,iin
Character(L_VAR)::elmt

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max
 iorg = VAR(ili)%idorg
 !-- limitation lumière

if (VAR(ili)%conc(i,j,k).le.0) then
  limn(1)%val(i,j,k) = 0.0
  limp(1)%val(i,j,k) = 0.0
   write(*,*)'pas de calcul de la limitation de l uptake par les nut dans geider' 
 
  f_ppb_geider(i,j,k) = 0.

else
  !-- limitation nutriments ("Droop")
!-- ATTENTION : l'élément qui suit le C est le N !!!!
!! et l'element qui suit le N est P   !!!!!
! lim gedier
limn(1)%val(i,j,k) = (VAR(ili+1)%conc(i,j,k)/VAR(ili)%conc(i,j,k)- qmin)/(qmax-qmin)
limp(1)%val(i,j,k) = (VAR(ili+2)%conc(i,j,k)/VAR(ili)%conc(i,j,k)- qpmin)/(qpmax-qpmin)
!limn droop
!limn(1)%val(i,j,k) = (VAR(ili+1)%conc(i,j,k)/VAR(ili)%conc(i,j,k)- qmin)/(VAR(ili+1)%conc(i,j,k)/VAR(ili)%conc(i,j,k))
!limp(1)%val(i,j,k) = (VAR(ili+2)%conc(i,j,k)/VAR(ili)%conc(i,j,k)- qpmin)/(VAR(ili+2)%conc(i,j,k)/VAR(ili)%conc(i,j,k))
!write(*,*)'limn,limp',limn,limp


!-- cas d'un rapport CHL:C constant:
  if (.not. CHL_C_BOOL) then
	Chl = CHL_C0 * var(ili)%conc(i,j,k)
	theta = CHL_C0 !-- Rapport Chl:C:
  else
    !-- Indice correspondant a la conc. chlorophylienne du phyto:
    elmt ='Chl'
    iichl = f_idorg2id (iorg,elmt)
    !-- Indice correspondant a la conc. en azote du phyto:
    elmt ='N'
    iin = f_idorg2id (iorg,elmt)
    N = var(iin)%conc(i,j,k)
    Chl = var(iichl)%conc(i,j,k)
    theta = Chl/(var(ili)%conc(i,j,k)+1.d-15)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ATTENTION AJOUT DU CALCUL DE LA CHLORO PAR VAN DER MEERSCHE

   var(iichl)%conc(i,j,k) = VAR(iin)%conc(i,j,k)/VAR(ili)%conc(i,j,k)&
    *(qminchl+(qmaxchl-qminchl)*((VAR(iin)%conc(i,j,k)/VAR(ili)%conc(i,j,k)-qmin)/(qmax-qmin))) &
         * VAR(ili)%conc(i,j,k)

  endif
!-- calcul de mumax en fonction de T et FQ (Pcmax dans geider,98)
!determination du terme de limitation le plus fort entre le P et le N
lim_nut(1)%val(i,j,k)=min(limn(1)%val(i,j,k),limp(1)%val(i,j,k))
!write(*,*)'lim_nut',lim_nut(1)%val(i,j,k)
 mumaxt = lim_nut(1)%val(i,j,k)*mumaxa*exp(mumaxb*TEMP_BIO(i,j,k))*1/(24.*60.*60.)

!-- taux de croissance (Pcphot dans geider,98)
 mu_phy(iorg)%val(i,j,k) = mumaxt*(1.-exp(-alpha*theta*E_PARZ(i,j,k)/mumaxt))

 !-- Croissance
 f_ppb_geider(i,j,k) = mu_phy(iorg)%val(i,j,k)*VAR(ili)%conc(i,j,k)
endif


 ! Stockage de la production primaire
PPB(iorg)%val(i,j,k) = f_ppb_geider(i,j,k)


if(f_ppb_geider(i,j,k).LT.0.)then
f_ppb_geider(i,j,k)=0.
endif



    enddo
  enddo
enddo
 Return
 End function f_ppb_geider
