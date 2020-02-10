!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!m.baklouti@univmed.fr; vincent.faure@univmed.fr
!
!This software (Eco3M) is a computer program whose purpose is to perform 
!biogeochemical or coupled physical-biogeochemical modelling.
!
!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software. You can  use, 
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software''s author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability. 
!
!In this respect, the user''s attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software''s suitability as regards their
!requirements in conditions enabling the security of their systems and/or 
!data to be ensured and,  more generally, to use and operate it in the 
!same conditions as regards security. 
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.
!***************************************************************************
!***************************************************************************
!
!----------------------------------------------------------------
!
                SUBROUTINE sub_lec_irrad
!
!-- Routine de lecture de l''irradiance dans le fichier fichirrad
!
! dernière modification: 27/06/07
!----------------------------------------------------------------
! lecture de l''irradiance dans le fichier ad hoc au temps considere
! ou calcul de l''irradiance par une fonction. 
! Une seule valeur de l''irradiance est fournie et affectee a tout le tableau irrad
!
Use COUPLEUR_PHY_BIO
Use VAR_USER
Use VAR_GLOBAL

Implicit None

!-- variables locales:
 Integer :: i,j,k,errlec,lec_irrad
 Integer :: id,nfich1,nfich2,nfich3
 real(8),PARAMETER  :: pi=3.1415d0
 Real(8) :: irradscal ! valeur scalaire de l''irradiance
 real(8) :: tval
 Logical :: present
 Character(10)  :: fich1,fich2,fich3
 Character(70)  :: fichapp


!-- cas d'une irradiance ou d'un PAR fournis par un fichier:
if (fichirrad /= 'IRR_FONCTION' .and. fichirrad /= 'IRR_CODEPHYS') Then
!-- premier appel du module bio: 
  if (nbcallbio == 0 ) then 
! test d''existence du fichier d''irradiance
     inquire (file=fichirrad_long,EXIST=present)
     if ((.NOT. present) ) then
      write(*,*) 'ATTENTION : fichier d irradiance',fichirrad_long,'non trouve'
      write(*,*) 'Le programme va etre arrete'
      stop
     else
      open (20,FILE=fichirrad_long) 
     endif
  endif
!   
!-- Lecture de l''irradiance :
!   Read(20,*,iostat=errlec) irradscal
   Read(20,*) irradscal
!Debug Christel
   if (irradscal < 0.d0) then !irradscal = 0.d0
     irradscal = 0.d0
!   if (errlec /=0) then 
     print *, 'pb de lecture du fichier IRRAD :', fichirrad
     stop
   endif
!--cas d''une irradiance calculee par une fonction (qui peut etre changee si besoin)  
elseif (fichirrad == 'IRR_FONCTION') then 
! irradscal  = irrad_MAX *exp(3.7*(cos(2.0*pi*tps/24/3600.)-1))
!-- variations saisonnieres de l'irradiance en Mediterranee
!-- début conditions hiver
 irradscal = max(1.d0,max(50.,abs(800.*cos(0.45*pi*((tps/3600.-4000.)/4100.)))) * sin(pi*tps/12./3600.));
!-- debut conditions été
! irradscal = max(1.d0,max(50.,abs(800.*sin(0.45*pi*((tps/3600.-4000.)/4100.)))) * sin(pi*tps/12./3600.));
!   irradscal=irrad_MAX
endif
!-- l''irradiance fournie est un scalaire que l'on affecte a tout l'espace horizontal:
  irrad = irradscal
!-- variations saisonnieres de la temp en Mediterranee
!-- début conditions hiver
!TEMP_BIO(:,:,:)= max(13.0,ABS(25.0*cos(0.45*pi*((tps/3600.-5000.)/4380.))))
!TEMP_BIO(:,:,:)= max(13.0,ABS(25.0*cos(0.45*pi*((tps/3600.-4000.)/4100.))))
!TEMP_BIO(:,:,:)= max(18.0,ABS(30.0*cos(0.45*pi*((tps/3600.-5000.)/4380.))))
!-- début condition été
!TEMP_BIO(:,:,:)= max(13.0,ABS(25.0*sin(0.45*pi*((tps/3600.-5000.)/4380.))))
!TEMP_BIO(:,:,:)= max(18.0,ABS(30.0*sin(0.45*pi*((tps/3600.-5000.)/4380.))))
!TEMP_BIO(:,:,:)=TEMP_BIO(:,:,:)+1.5
!-- Lecture temperature dans fichier --!
open (21,FILE='BIO/DATA/Temperature/Planier_TEMPESURF_2017_6ans_10min.txt') 
!open (21,FILE='BIO/DATA/Temperature/TEMPESURF+15_10min.txt') 
  read(21,*)tval
  TEMP_BIO(:,:,:)=tval
!  write(*,*)tval,TEMP_BIO(:,:,:)
open (22,FILE='BIO/DATA/Salinite/Smatch_CARRY_SALSURF_2017_6ans_10min.txt') 
!open (22,FILE='BIO/DATA/Salinite/Sal.txt') 
  read(22,*)tval
  SAL_BIO(:,:,:)=tval
!SAL_BIO(:,:,:)=38.0
open (25,FILE='BIO/DATA/Wind/WRF_WIND_2017_6ans_10min.txt') 
!open (25,FILE='BIO/DATA/Wind/Wind_event.txt')
  read(25,*)tval
  WIND(:,:,:)=tval
!  WIND(:,:,:)=7.
open (26,FILE='BIO/DATA/CO2/CO2_CAV_2017_6ans_10min_atm.txt')
  read(26,*)tval
  CO2atm(:,:,:)=tval
!  CO2atm(:,:,:)=410.
open (27,FILE='BIO/DATA/CO2/Density_sw_2017_6ans_10min.txt')
  read(27,*)tval
  rho_sw(:,:,:)=tval
! Leceture fichier d'apports (Ajout Katixa Lajaunie-Salla 20/03/2019)
do id=1,nbvar
   fich1 = VAR(id)%comp
   nfich1 = len_trim(adjustl(fich1))
   fich2 = VAR(id)%scomp
   nfich2 = len_trim(adjustl(fich2))
   fich3 = VAR(id)%elmt
   nfich3 = len_trim(adjustl(fich3))
   fichapp = trim(repconcapp) // fich1(1:nfich1)//'_'//fich2(1:nfich2)//'_'//fich3(1:nfich3)//'_apport.DAT'
!-- boucle de lecture
   open(400+id,FILE=fichapp)
   read(400+id,*)tval 
   capp(id,:,:,:)= tval
enddo
!stop

#ifdef MODTEST
  Write(*,*) '---  sub_lec_irrad terminée '
#endif
End Subroutine sub_lec_irrad
!----------------------------------------------------------------
!
!
                Subroutine sub_calc_EPARZ
!
! subroutine de calcul de E_PAR(z) sur la colonne d''eau 
! irrad est le rayonnement incident global (W/m2)
! PAR0plus est la portion d''irrad comprise dans la fenetre 400-700 nm
!
! dernière modification: 27/06/07
!----------------------------------------------------------------
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use VAR_USER
Use MOD_FCHAIN
Use F_PROCESS

Implicit None


!-- variables locales:
 Integer    :: iphy_Chl(nscp_phy),iphy_C(nscp_phy)
 Integer    :: i,j,k,nn
 real(8)    :: Chl_tot
 real(8)    :: kextinc
 real(8)    :: dz
 character(L_CHAIN) :: scomp,elmt
 character(3) :: debfichirrad


!-------------------------------------------------------------
!-- Calcul de E_PAR(0+) sauf si irrad est deja un PAR
!-------------------------------------------------------------
 fichirrad = adjustl(fichirrad)
 debfichirrad= fichirrad(1:3)
 if (debfichirrad /= 'PAR') then 
         E_PAR(:,:) = irr2par * irrad !integrale du spectre sur les longueurs d''onde utilisables
 else
         E_PAR(:,:) =  irrad    
 endif	
!-------------------------------------------------------
!-- Calcul de E_PAR(0-) 
!-------------------------------------------------------
!-- Modele Général :  
   E_PAR(:,:) = irr_param * (1-albedo) * E_PAR(:,:) 
write(*,*)'test irrad',E_PAR(:,:),irrad,TEMP_BIO(:,:,:)
!-------------------------------------------------------
!-- Calcul du terme d''extinction (sauf en 0D)
!-------------------------------------------------------
!  
#ifndef M0D-NCOUPL
! initialisation:
iphy_C = 0
iphy_Chl = 0
!
! Le calcul du coefficient d''extinction necessite de connaitre la biomasse chlorophylienne:
!
! 1/ determination des indices de variable contenant des conc en Chl:
  ! 1/ cas ou la conc en Chl est calculee explicitement
  if (CHL_C_BOOL) then
      elmt= 'Chl'
      iphy_Chl =  f_idorg2id_vect(iscp_phy,elmt,nscp_phy)
!      write(11,*) 
!      write(11,*) 'iphy_Chl =',iphy_Chl
!      write(11,*) 
  else
   ! 2/ cas ou la conc en Chl est calculee a partir du Carbone
      elmt='C'
      iphy_C   =  f_idorg2id_vect(iscp_phy,elmt,nscp_phy)
  endif
  
!-- calcul de Chl_tot pour chaque maille:
!----------------------------------------
!- cas ou nz_max est à la surface:
   if (pos_nzmax == 'SURF') then  
!----------------------------------------
    do i=nx_min,nx_max
     do j=ny_min,ny_max
      do k=1,nz_max
! mise a zero de la Chl tot sur chaque maille
         Chl_tot = 0.d0
!--calcul de la conc totale en Chl a la maille i,j,k:
      do nn=1,nscp_phy
       if (CHL_C_BOOL .and. maxval(iphy_Chl) /=0) then
         if (k /= nz_max) then
           Chl_tot = Chl_tot + (VAR(iphy_Chl(nn))%conc(i,j,k) + VAR(iphy_Chl(nn))%conc(i,j,k+1))/2.d0
	 else
           Chl_tot = Chl_tot + VAR(iphy_Chl(nn))%conc(i,j,k) 
	 endif
       elseif(.NOT. (CHL_C_BOOL) .and. maxval(iphy_C) /=0) then
         if (k /= nz_max) then
	   Chl_tot = Chl_tot + (VAR(iphy_C(nn))%conc(i,j,k) + VAR(iphy_C(nn))%conc(i,j,k+1))/2.d0*CHL_C0
	 else
           Chl_tot = Chl_tot + VAR(iphy_Chl(nn))%conc(i,j,k) 
	 endif	        
      endif
      enddo
!--------------------------
!---- Calcul de E_PAR(z) 
!--------------------------
       if (k==nz_max) then ! pres de la surface
!-- Calcul de la distance entre deux noeuds 
	  dz = prof(i,j) - alti_z(i,j,nz_max)
       else
          dz = alti_z(i,j,k+1) - alti_z(i,j,k) 
       endif
!-- lignes de code permettant de calculer le terme d''extinction kextinc
!  (Chl_tot et dz y sont transmis par argument)
!
#ifdef CALC
        include "calc_extinc.inc"
#endif	
	if (k== nz_max) then
           if (size(E_PAR) == 1) then 
             E_PARZ(i,j,k) = E_PAR(1,1) * kextinc
           else
             E_PARZ(i,j,k) = E_PAR(i,j) * kextinc
           endif
        else
           E_PARZ(i,j,k) = E_PARZ(i,j,k+1) * kextinc
	endif
       enddo         
      enddo
     enddo
!- cas ou nz_max est à la surface:
!----------------------------------------
     elseif (pos_nzmax == 'FOND') then   
!----------------------------------------
    do i=nx_min,nx_max
     do j=ny_min,ny_max
      do k=1,nz_max
! mise a zero de la Chl tot sur chaque maille
         Chl_tot = 0.d0
!--calcul de la conc totale en Chl a la maille i,j,k:
      do nn=1,nscp_phy
       if (CHL_C_BOOL .and. maxval(iphy_Chl) /=0) then
         if (k /= 1) then
           Chl_tot = Chl_tot + (VAR(iphy_Chl(nn))%conc(i,j,k) + VAR(iphy_Chl(nn))%conc(i,j,k-1))/2.d0
	 else
           Chl_tot = Chl_tot + VAR(iphy_Chl(nn))%conc(i,j,k) 
	 endif
       elseif(.NOT. (CHL_C_BOOL) .and. maxval(iphy_C) /=0) then
         if (k /= 1) then
	   Chl_tot = Chl_tot + (VAR(iphy_C(nn))%conc(i,j,k) + VAR(iphy_C(nn))%conc(i,j,k-1))/2.d0*CHL_C0
	 else
           Chl_tot = Chl_tot + VAR(iphy_Chl(nn))%conc(i,j,k) 
	 endif	        
      endif
      enddo
      if (k==1) then 
!-- Calcul de la distance entre deux noeuds 
         dz = prof(i,j) - alti_z(i,j,1)
      else 
         dz = alti_z(i,j,k-1) - alti_z(i,j,k) 
      endif
!-- lignes de code permettant de calculer le terme d''extinction kextinc
!  (Chl_tot et dz y sont transmis par argument)
#ifdef CALC
       include "calc_extinc.inc"
#endif
	if (k==1) then
           if (size(E_PAR) == 1) then 
             E_PARZ(i,j,k) = E_PAR(1,1) * kextinc
           else
             E_PARZ(i,j,k) = E_PAR(i,j) * kextinc
           endif
        else
	   E_PARZ(i,j,k) = E_PARZ(i,j,k-1) * kextinc  
	endif

        enddo         
       enddo
     enddo
    endif
#else
  do k=1,nz_max
   E_PARZ(:,:,k) = E_PAR(:,:)
  enddo
#endif
!
#ifdef MODTEST
  Write(*,*) '---  sub_calc_EPARZ terminée '
#endif
!!-- Expressions à transformer en fonctions (très prochainement):
!
! expression David Antoine (1998) sauf le 0.026 qui est une moyenne des valeurs de a*:
!      kk = kk0 + 0.032 * (Chl_tot*1.d3)**0.8  !-- relation pour Chl en mg m-3
!
! A VOIR (morel ???):
! kk = kk0 + 0.0088 * Chl_tot*1.d3  + 0.054 * (Chl_tot*1.d3)**(2./3))

End Subroutine sub_calc_EPARZ
!----------------------------------------------------------------
!
                   SUBROUTINE sub_calc_Chl_C
!
!-- Calcul du rapport Chl:C (à partir des concentrations au 
!   pas de tps precedent)
!  Faire de cette subroutine une fonction plutôt
!
!-- dernière modification: 27/06/07
!----------------------------------------------------------------
!
! parmi les variables d''etat, on recherche le phyto:
! pour chaque classe de phyto, on cherche les conc. en C et en Chl:

Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use MOD_FCHAIN

Implicit None

!-- variables globales:
 Character(L_VAR) :: chaine,chaine2
 integer :: i, ichl, ic ,iorg

if (Allocated (CHL_C)) then  
  ic = 0
  ichl = 0

  do i = 1,nscp_phy 
      iorg = CHL_C(i)%idorg
      chaine = 'C'
      ic = f_idorg2id(iorg,chaine)
      chaine2 = 'Chl'
      ichl = f_idorg2id(iorg,chaine2)
      CHL_C(i)%val(:,:,:) = var(ichl)%conc(:,:,:)/(var(ic)%conc(:,:,:) + 1.d-15)
  enddo
endif

#ifdef MODTEST
  write(*,*) 'sub_calc_Chl_C terminée'
#endif
  
End Subroutine sub_calc_Chl_C
!
!----------------------------------------------------------------
!
                SUBROUTINE sub_calc_mu_phy
		
!
!-- Calcul de la vitesse de croissance (s^-1)
!   pour chaque phytoplancton 
!
!-- dernière modification: 27/06/07
!----------------------------------------------------------------
!
 Use COUPLEUR_PHY_BIO
 Use VAR_GLOBAL
 Use MOD_FCHAIN

 Implicit None

!-- variables globales:
 Character(L_VAR) :: chaine
 integer :: i, ic,  iorg
 
 
 do i = 1,nscp_phy
     iorg = mu_phy(i)%idorg 
     chaine = 'C'
     ic = f_idorg2id(iorg,chaine)
     mu_phy(i)%val(:,:,:) = SELF_VAL(ic)%val / (var(ic)%conc + 1.d-15)
 enddo
End Subroutine sub_calc_mu_phy
!----------------------------------------------------------------
!----------------------------------------------------------------
!
                   SUBROUTINE sub_calc_CO2
!
!-- Calcul de l'equilibre des carbonates 
! derniere modification 02/04/2019 : Katixa Lajaunie-Salla
!----------------------------------------------------------------
!
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use MOD_FCHAIN
Use F_PROCESS

Implicit None
 integer :: ivar,i,j,k
 real    :: pHTol,R,P1atm
 real(8) :: Denom(nx_min:nx_max,ny_min:ny_max,nz_max),CAlk(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: BAlk(nx_min:nx_max,ny_min:ny_max,nz_max),TB(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: TF(nx_min:nx_max,ny_min:ny_max,nz_max),TS(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: KF(nx_min:nx_max,ny_min:ny_max,nz_max),KS(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: IonS(nx_min:nx_max,ny_min:ny_max,nz_max),FREEtoTOT(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: SWStoTOT(nx_min:nx_max,ny_min:ny_max,nz_max),Hfree(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: HSO4(nx_min:nx_max,ny_min:ny_max,nz_max),HF(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: Residual(nx_min:nx_max,ny_min:ny_max,nz_max),Slope(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: deltapH(nx_min:nx_max,ny_min:ny_max,nz_max),FugFac(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: Delta(nx_min:nx_max,ny_min:ny_max,nz_max),b(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: lnK1fac(nx_min:nx_max,ny_min:ny_max,nz_max),lnK2fac(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: lnKBfac(nx_min:nx_max,ny_min:ny_max,nz_max),lnKefac(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: lnKSfac(nx_min:nx_max,ny_min:ny_max,nz_max),lnKFfac(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: lnKcafac(nx_min:nx_max,ny_min:ny_max,nz_max),Kca(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: Pbar(nx_min:nx_max,ny_min:ny_max,nz_max),pHfactor(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: AT_test1(nx_min:nx_max,ny_min:ny_max,nz_max),AT_test2(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8) :: Nut(nx_min:nx_max,ny_min:ny_max,nz_max)

do ivar=1,nbvar
   if (VAR(ivar)%elmt == 'dic') then

VAR(ivar)%conc(:,:,:)=VAR(ivar)%conc(:,:,:)/(rho_sw(:,:,:)*1.d-3)/1.d6 ! DIC en mol/kg
VAR(ivar+3)%conc(:,:,:)=VAR(ivar+3)%conc(:,:,:)/(rho_sw(:,:,:)*1.d-3)/1.d6 ! AT en mol/kg

! Temperature in Kelvin
 TEMP(:,:,:)=TEMP_BIO(:,:,:)+273.15

!!--Calcul of constants--!!

! Calculate TF from Riley, J. P., Deep-Sea Research 12:219-220, 1965:
! 0.000068.*Sali./35. = 0.00000195.*Sali
 TF(:,:,:)=(0.000067/18.998)*(SAL_BIO(:,:,:)/1.80655) ! in mol/kg-SW

! Calculate TS from Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
! 0.02824.*Sali./35. = .0008067.*Sali
 TS(:,:,:)=(0.14/96.062)*(SAL_BIO(:,:,:)/1.80655) ! in mol/kg-SW

! Concentration I : from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
 IonS(:,:,:)=19.924*SAL_BIO(:,:,:)/(1000.-1.005*SAL_BIO(:,:,:));

! Concnetration Tot borate from Uppstrom, L., Deep-Sea Research 21:161-162, 1974
 TB(:,:,:)= 0.000416*SAL_BIO(:,:,:)/35 ! in mol/kg-SW

! KS: constante de dissolution HSO4 from Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
! in mol/kg-H2O ==> convert it to mol/kg-SW : lnKsw=lnK+ln(1-0.001005*Sal) in free pH scale
 KS(:,:,:)=-4276.1/TEMP(:,:,:)+141.328-23.093*log(TEMP(:,:,:))+(324.57-47.986*log(TEMP(:,:,:))-13856./TEMP(:,:,:))*(IonS(:,:,:)**0.5)
 KS(:,:,:)=KS(:,:,:)+(-771.54+114.723*log(TEMP(:,:,:))+35474./TEMP(:,:,:))*IonS(:,:,:)+&
 (-2698./TEMP(:,:,:))*(IonS(:,:,:)**1.5)+(1776./TEMP(:,:,:))*(IonS(:,:,:)**2.)
 KS(:,:,:)=exp(KS(:,:,:))*(1.-0.001005*SAL_BIO(:,:,:))

! KF: Constante formation HF from Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
! in mol/kg-H2O ==> convert it to mol/kg-SW : lnKsw=lnK+ln(1-0.001005*Sal) in free pH scale
 KF(:,:,:)=exp(1590.2/TEMP(:,:,:)-12.641+1.525*IonS(:,:,:)**(0.5))*(1.-0.001005*SAL_BIO(:,:,:))

! Calculate pH Scale Conversion Factors (These are NOT pressure-corrected)
 SWStoTOT(:,:,:)=(1.+TS(:,:,:)/KS(:,:,:))/(1.+TS(:,:,:)/KS(:,:,:)+TF(:,:,:)/KF(:,:,:));
 FREEtoTOT(:,:,:)=1+TS(:,:,:)/KS(:,:,:);

! KB  (Dickson, A. G., Deep-Sea Research 37:755-766,1990) total scale mol/kg-SW
 KB(:,:,:)=(-8966.9-2890.53*(SAL_BIO(:,:,:)**0.5)-77.942*SAL_BIO(:,:,:)+1.728*(SAL_BIO(:,:,:)**(3./2))-0.0996*(SAL_BIO(:,:,:)**2.))/(TEMP(:,:,:))
 KB(:,:,:)=KB(:,:,:)+148.0248+137.1942*(SAL_BIO(:,:,:)**0.5)+1.62142*SAL_BIO(:,:,:)+&
 (-24.4344-25.085*(SAL_BIO(:,:,:)**0.5)-0.2474*SAL_BIO(:,:,:))*LOG(TEMP(:,:,:))+0.053105*(SAL_BIO(:,:,:)**0.5)*TEMP(:,:,:)
 KB(:,:,:)=exp(KB(:,:,:))/SWStoTOT(:,:,:) ! convert to SWS pH scale

! K0 constante de solubilite CO2 dans l'atmosphere from Weiss, R. F., Marine Chemistry 2:203-215, 1974.
! this is in mol/kg-SW/atm
 K0(:,:,:)=exp(-60.2409+93.4517*(100./TEMP(:,:,:))+23.3585*log(TEMP(:,:,:)/100.)+SAL_BIO(:,:,:)*&
 (0.023517-0.023656*(TEMP(:,:,:)/100.)+0.0047036*(TEMP(:,:,:)/100.)**2.)) 

! Ke: produit ionique de l'eau from Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
! this is on the SWS pH scale in (mol/kg-SW)^2
 Ke(:,:,:)=EXP(-13847.26/(TEMP(:,:,:))+148.9802-23.6521*LOG(TEMP(:,:,:))+&
 (-5.977+118.67/(TEMP(:,:,:))+1.0495*LOG(TEMP(:,:,:)))*SAL_BIO(:,:,:)**(0.5)-0.01615*SAL_BIO(:,:,:))

! K1 & K2: From Lueker, Dickson, Keeling, Mar. Chem. 70 (2000) 105-119
! Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work. 
! this is on the Total scale in mol/kg-SW
 K1(:,:,:)=10.**(-3633.86/(TEMP(:,:,:))+61.2172-9.6777*LOG(TEMP(:,:,:))+0.011555*SAL_BIO(:,:,:)-0.0001152*(SAL_BIO(:,:,:)**2.))
 K2(:,:,:)=10.**(-471.78/(TEMP(:,:,:))-25.929+3.16967*LOG(TEMP(:,:,:))+0.01781*SAL_BIO(:,:,:)-0.0001122*SAL_BIO(:,:,:)**2.)
! conversion in SWS pH scale mol/kg-SW
 K1(:,:,:)=K1(:,:,:)/SWStoTOT(:,:,:)
 K2(:,:,:)=K2(:,:,:)/SWStoTOT(:,:,:)

! Ksp for calcite from: Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983
! this is in (mol/kg-SW)^2
 Kca(:,:,:)=10**(-171.9065-0.077993*TEMP(:,:,:)+2839.319/(TEMP(:,:,:))+71.595*LOG10(TEMP(:,:,:))+(-0.77712+0.0028426*TEMP(:,:,:)&
    +178.34/TEMP(:,:,:))*(SAL_BIO(:,:,:)**(0.5))-0.07711*SAL_BIO(:,:,:)+0.0041249*(SAL_BIO(:,:,:)**(1.5)))

! Ca2+ concentration in mol/kg from: Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
! this is .010285.*Sali./35 ==> in mol/kg-SW
 Ca2(:,:,:)=0.02128/40.087*(SAL_BIO(:,:,:)/1.80655)

!if (pHscale==1) then ! SWS pH scale
!  pHfactor=1.
!else ! TOT pH scale
  pHfactor(:,:,:)=SWStoTOT(:,:,:)
!endif

K1(:,:,:) = K1(:,:,:)*pHfactor
K2(:,:,:) = K2(:,:,:)*pHfactor
Ke(:,:,:) = Ke(:,:,:)*pHfactor
KB(:,:,:) = KB(:,:,:)*pHfactor

!!--Corrected by Pressure--!!
R = 83.1451 ! ml/bar/K/mol
Pbar(:,:,:)=1.
lnK1fac(:,:,:)=(25.5-0.1271*TEMP_BIO(:,:,:)+0.5*((-3.08+0.0877*TEMP_BIO(:,:,:))/1000)*Pbar(:,:,:))*Pbar(:,:,:)/(R*TEMP(:,:,:))
lnK2fac(:,:,:)=(15.82-0.0219*TEMP_BIO(:,:,:)+0.5*((1.13-0.1475*TEMP_BIO(:,:,:))/1000)*Pbar(:,:,:))*Pbar(:,:,:)/(R*TEMP(:,:,:))
lnKBfac(:,:,:)=(29.48-0.1622*TEMP_BIO(:,:,:)+0.002608*TEMP_BIO(:,:,:)**2.+0.5*(-2.84/1000)*Pbar(:,:,:))*Pbar(:,:,:)/(R*TEMP(:,:,:))
lnKefac(:,:,:)=(20.02-0.1119*TEMP_BIO(:,:,:)+0.001409*TEMP_BIO(:,:,:)**2.+0.5*((-5.13+0.0794*TEMP_BIO(:,:,:))/1000)*Pbar(:,:,:))*Pbar(:,:,:)/(R*TEMP(:,:,:))
lnKFfac(:,:,:)=(9.78-0.009*TEMP_BIO(:,:,:)+0.000942*TEMP_BIO(:,:,:)**2.+0.5*((-3.91+0.054*TEMP_BIO(:,:,:))/1000)*Pbar(:,:,:))*Pbar(:,:,:)/(R*TEMP(:,:,:))
lnKSfac(:,:,:)=(18.03-0.0466*TEMP_BIO(:,:,:)+0.000316*TEMP_BIO(:,:,:)**2.+0.5*((-4.53+0.009*TEMP_BIO(:,:,:))/1000)*Pbar(:,:,:))*Pbar(:,:,:)/(R*TEMP(:,:,:))
lnKcafac(:,:,:)=(48.76-0.5304*TEMP_BIO(:,:,:)+0.5*((-11.76+0.3692*TEMP_BIO(:,:,:))/1000)*Pbar(:,:,:))*Pbar(:,:,:)/(R*TEMP(:,:,:))

K1(:,:,:)=K1(:,:,:)*exp(lnK1fac(:,:,:))
K2(:,:,:)=K2(:,:,:)*exp(lnK2fac(:,:,:))
KB(:,:,:)=KB(:,:,:)*exp(lnKBfac(:,:,:))
Ke(:,:,:)=Ke(:,:,:)*exp(lnKefac(:,:,:))
KF(:,:,:)=KF(:,:,:)*exp(lnKFfac(:,:,:))
KS(:,:,:)=KS(:,:,:)*exp(lnKSfac(:,:,:))
Kca(:,:,:)=Kca(:,:,:)*exp(lnKcafac(:,:,:))

!!--Calcul Fugacity Constant--!!
! This assumes that the pressure is at one atmosphere, or close to it: Weiss, R. F., Marine Chemistry 2:203-215, 1974.
Delta(:,:,:)=(57.7-0.118*TEMP(:,:,:)) ! in cm3/mol
b(:,:,:)=-1636.75+12.0408*TEMP(:,:,:)-0.0327957*TEMP(:,:,:)**2.+3.16528*0.00001*TEMP(:,:,:)**3 ! in cm3/mol
! For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
P1atm=1.01325 ! in bar
FugFac(:,:,:)=exp((b(:,:,:)+2.*Delta(:,:,:))*P1atm/(R*TEMP(:,:,:)))

! Nut dans AT en mol/kg
Nut(:,:,:)=(VAR(20)%conc(:,:,:)+VAR(19)%conc(:,:,:)-VAR(18)%conc(:,:,:))/(rho_sw(:,:,:)*1.d-3)/1.d6

do i=nx_min,nx_max
  do j=ny_min,ny_max
    do k=1,nz_max
if (nbcallbio <1.) VAR(ivar+2)%conc(i,j,k)=8.
pHTol       = 0.0001  ! tolerance for iterations end
deltapH(i,j,k)    = pHTol+1
do while (abs(deltapH(i,j,k)) > 0.0001)
    H(i,j,k)        = 10.**(-VAR(ivar+2)%conc(i,j,k))
    Denom(i,j,k)    = H(i,j,k)**2.+K1(i,j,k)*H(i,j,k)+K1(i,j,k)*K2(i,j,k)
    CAlk(i,j,k)     = VAR(ivar)%conc(i,j,k)*K1(i,j,k)*(H(i,j,k)+2.*K2(i,j,k))/Denom(i,j,k)
    BAlk(i,j,k)     = TB(i,j,k)*KB(i,j,k)/(KB(i,j,k)+H(i,j,k))
    OH(i,j,k)       = Ke(i,j,k)/H(i,j,k)
    FREEtoTOT(i,j,k)= (1.+TS(i,j,k)/KS(i,j,k)) ! pH scale conversion factor
    Hfree(i,j,k)    = H(i,j,k)/FREEtoTOT(i,j,k) ! for H on the total scale
    HSO4(i,j,k)     = TS(i,j,k)/(1.+KS(i,j,k)/Hfree(i,j,k)) ! since KS is on the free scale
    HF(i,j,k)       = TF(i,j,k)/(1.+KF(i,j,k)/Hfree(i,j,k)) ! since KF is on the free scale
    Residual(i,j,k) = VAR(ivar+3)%conc(i,j,k)-CAlk(i,j,k)-BAlk(i,j,k)-OH(i,j,k)+Hfree(i,j,k)+HSO4(i,j,k)+HF(i,j,k)!-Nut(i,j,k)! - PAlk - SiAlk ;
    ! find Slope dTA/dpH (see Middelburg 2019 table 5.2)
    Slope(i,j,k)    = VAR(ivar)%conc(i,j,k)*H(i,j,k)*K1(i,j,k)*(H(i,j,k)**2.+K1(i,j,k)*K2(i,j,k)+4.*H(i,j,k)*K2(i,j,k))
    Slope(i,j,k)    = (Slope(i,j,k)/(Denom(i,j,k)**2.))+OH(i,j,k)+H(i,j,k)+(BAlk(i,j,k)*H(i,j,k)/(KB(i,j,k)+H(i,j,k)))
    Slope(i,j,k)    = log(10.)*Slope(i,j,k)
    deltapH(i,j,k)  = Residual(i,j,k)/Slope(i,j,k) ! this is Newton's method
    ! to keep the jump from being too big;
    do while (abs(deltapH(i,j,k) > 1.))
        deltapH(i,j,k)=deltapH(i,j,k)/2.
    enddo

VAR(ivar+2)%conc(i,j,k)=VAR(ivar+2)%conc(i,j,k)+deltapH(i,j,k) ! Is on the same scale as K1 and K2 were calculated...

enddo
enddo
enddo
enddo  

VAR(ivar+1)%conc(:,:,:)=(VAR(ivar)%conc(:,:,:)*(H(:,:,:)**2.))/(H(:,:,:)**2.+K1(:,:,:)*H(:,:,:)+K1(:,:,:)*K2(:,:,:))
VAR(ivar+1)%conc(:,:,:)=VAR(ivar+1)%conc(:,:,:)*1.d6/K0(:,:,:)/FugFac(:,:,:) ! in µatm
!write(*,*)'test CO2sys calc pco2',VAR(ivar+1),VAR(ivar),FugFac(:,:,:),K0(:,:,:)

!!--Calcul of others variables of carbonates system--!! in µmol/kg
CO2(:,:,:)   = VAR(ivar)%conc(:,:,:)*1.d6/(1.d0+K1(:,:,:)/H(:,:,:)+K1(:,:,:)*K2(:,:,:)/H(:,:,:)/H(:,:,:))
HCO3(:,:,:)  = K1(:,:,:)*CO2(:,:,:)/H(:,:,:)
CO3(:,:,:)   = K2(:,:,:)*HCO3(:,:,:)/H(:,:,:)
Omega(:,:,:) = Ca2(:,:,:)*CO3(:,:,:)*1.d-6/Kca(:,:,:)

AT_test(:,:,:)=HCO3(:,:,:)+2*CO3(:,:,:)+(Balk(:,:,:)+OH(:,:,:)-Hfree(:,:,:)-HF(:,:,:)-HSO4(:,:,:))*1.d6
AT_test1(:,:,:)=AT_test(:,:,:)+Nut(:,:,:)*1.d6

VAR(ivar)%conc(:,:,:)=VAR(ivar)%conc(:,:,:)*(rho_sw(:,:,:)*1.d-3)*1.d6 ! en µmol/L
VAR(ivar+3)%conc(:,:,:)=VAR(ivar+3)%conc(:,:,:)*(rho_sw(:,:,:)*1.d-3)*1.d6 ! en µmol/L
!write(*,*)'test CO2sys calc var',VAR(ivar),VAR(ivar+1),VAR(ivar+2),VAR(ivar+3)

if (nbcallbio == 1) then
open (144, FILE='./BIO/SORTIES/save_CT.sim')
write(144,'(a)') '%tps CT  pCO2 pH AT CaCO3 CO2 HCO3 CO3 Omega ATtest'
endif
write(144,200) tps/3600,CO2(:,:,:),HCO3(:,:,:),CO3(:,:,:),AT_test(:,:,:),AT_test1(:,:,:),Residual(:,:,:),Nut(:,:,:)*1.d6
200 FORMAT('',F9.3,1X,F7.4,1X,F7.2,1X,F7.2,1X,F9.4,1X,F7.2,1X,F11.10,1X,F7.3)
endif
enddo
end SUBROUTINE sub_calc_CO2
