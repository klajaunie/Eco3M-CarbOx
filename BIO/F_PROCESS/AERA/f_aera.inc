!----------------------------------------------------------------
!
                   function f_aera(ili,jcol,k600)
!----------------------------------------------------------------
!
! Creation Katixa Lajaunie-Salla 11/06/2018
! Calcul de l'aeration: echanges de gaz entre l'eau et l'atmosphère
!
Use COUPLEUR_PHY_BIO
Use  VAR_GLOBAL

Implicit None
!
!-- variables globales:
 REAL(8):: f_aera(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: S_flux_aera(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: O2_0,Uw,Vw,k600,Csat
! REAL(8):: K0(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: Catm(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: O2_sat(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: A(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: sc(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: sc35(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: sc0(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: kw(nx_min:nx_max,ny_min:ny_max,nz_max)
 REAL(8):: wspd(nx_min:nx_max,ny_min:ny_max,nz_max)
! REAL(8):: TEMP(nx_min:nx_max,ny_min:ny_max,nz_max)
integer :: i,j,k,ili,jcol,ivar


do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

! Temperature en Kelvin
!TEMP(i,j,k)=TEMP_BIO(i,j,k)+273.15

! calcul de la composante de vent
!wspd(i,j,k) = (Uw**2.+Vw**2.)**0.5 
wspd(i,j,k) = WIND(i,j,k) ! lu dans fichier sub_calc_prelim

if (VAR(ili)%elmt=='O2') then
! calcul de l'oxygene à saturation
 A(i,j,k)=0.007991-0.0000374*SAL_BIO(i,j,k)-0.000077774*TEMP_BIO(i,j,k)
 O2_sat(i,j,k)=14.652-0.0841*SAL_BIO(i,j,k)+(-0.41022+0.00256*SAL_BIO(i,j,k)+A(i,j,k)*TEMP_BIO(i,j,k))*TEMP_BIO(i,j,k)
 O2_sat(i,j,k)=O2_sat(i,j,k)*1000./32. ! mmol/m3 ou µmol/L 
 Catm(i,j,k)=O2_sat(i,j,k) 
! calcul du nombre de Schmidt number a salinite 35 et salinite 0
 sc35(i,j,k) = 1953.4-128.00*TEMP_BIO(i,j,k)+3.9918*(TEMP_BIO(i,j,k))**2.-0.050091*(TEMP_BIO(i,j,k))**3.
 sc0(i,j,k) = 1800.6-120.1*TEMP_BIO(i,j,k)+3.7818*(TEMP_BIO(i,j,k))**2.-0.047608*(TEMP_BIO(i,j,k))**3.
 sc(i,j,k)=(sc35(i,j,k)/sc0(i,j,k))*SAL_BIO(i,j,k)+sc0(i,j,k)
! Calcul du coefficient de transfert O2 (vitesse de piston en cm h-1) 
 kw(i,j,k)=(k600*wspd(i,j,k)**2.*(660./sc(i,j,k))**0.5)*(0.01/3600.) ! (on divise par 1m) ==> s-1
!Calcul du flux de surface en mmol/m3/s 
 S_flux_aera(i,j,k) = (Catm(i,j,k)-VAR(ili)%conc(i,j,k))*kw(i,j,k)

elseif (VAR(ili)%elmt=='dic') then
! K0 constante de solubilite CO2 dans l'atmosphere en mol/L/atm (Weiss, 1974)
 K0(i,j,k)=exp(-58.0931+90.5069*(100./TEMP(i,j,k))+22.2940*log(TEMP(i,j,k)/100.)+SAL_BIO(i,j,k)*(0.027766-0.025888*(TEMP(i,j,k)/100.)+0.0050578*(TEMP(i,j,k)/100.)**2.))
! Catm(i,j,k)=Csat 
 Catm(i,j,k)=CO2atm(i,j,k) ! en µatm
! calcul du nombre de Schmidt number a salinite 35 et salinite 0
 sc35(i,j,k)=2073.1-125.62*TEMP_BIO(i,j,k)+3.6276*(TEMP_BIO(i,j,k))**2.-0.043219*(TEMP_BIO(i,j,k))**3.
 sc0(i,j,k)=1911.1-118.11*TEMP_BIO(i,j,k)+3.4527*(TEMP_BIO(i,j,k))**2.-0.041320*(TEMP_BIO(i,j,k))**3.
 sc(i,j,k)=(sc35(i,j,k)/sc0(i,j,k))*SAL_BIO(i,j,k)+sc0(i,j,k)
! Calcul du coefficient de transfert O2 (vitesse de piston en cm h-1) 
 kw(i,j,k)=(k600*wspd(i,j,k)**2.*(660./sc(i,j,k))**0.5)*(0.01/3600.) ! (on divise par 1m) ==> s-1
 kw(i,j,k)=kw(i,j,k)*K0(i,j,k) ! mol/L/atm/s
!Calcul du flux de surface en µmol/L/s 
S_flux_aera(i,j,k) = (Catm(i,j,k)-VAR(ili+1)%conc(i,j,k))*kw(i,j,k) !en µmol/L/s

endif

if (nbcallbio <=1 ) then
  S_flux_aera(i,j,k)=0.d0
endif

f_aera(i,j,k)=S_flux_aera(i,j,k)


    enddo
  enddo
enddo

Return
end function f_aera


