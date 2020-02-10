!************************************************************************************************************
! --------- Modèles de limiation du grazing par la temperature------------
! ------------------------------------------------------------------------
!
              function f_lim_graz(ili)
!
! ------------------------------------------------------------------------
!
! Cette fonction permet de calculer la limitation du grazing par la temperature (via les concentration de phy et det)
!  g=(Proie)i/somme (Phyto+Det) ; ex graz du Pico ==> g=Pico/(Pico+Dia+Det)
! 

Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use  VAR_GLOBAL

Implicit None
!-- fonction
Real(8) :: f_lim_graz(nx_min:nx_max,ny_min:ny_max,nz_max)
integer :: i,j,k,ili,jcol,iorg
real(8) :: b,tref,tot
Real(8) :: var1(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: var2(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: var3(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: limgraz(nx_min:nx_max,ny_min:ny_max,nz_max)

iorg = VAR(ili)%idorg


do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max

if(VAR(ili)%elmt=="C") then
  var1(i,j,k)=VAR(1)%conc(i,j,k)
  var2(i,j,k)=VAR(5)%conc(i,j,k)
  var3(i,j,k)=VAR(12)%conc(i,j,k)
elseif(VAR(ili)%elmt=="N") then
  var1(i,j,k)=VAR(2)%conc(i,j,k)
  var2(i,j,k)=VAR(6)%conc(i,j,k)
  var3(i,j,k)=VAR(13)%conc(i,j,k)
elseif(VAR(ili)%elmt=="P") then
  var1(i,j,k)=VAR(3)%conc(i,j,k)
  var2(i,j,k)=VAR(7)%conc(i,j,k)
  var3(i,j,k)=VAR(14)%conc(i,j,k)
endif

tot=var1(i,j,k)+var2(i,j,k)+var3(i,j,k)
if (tot.ge.1.e-10) then
  if(VAR(ili)%scomp=="pico") then
    limgraz(i,j,k)=var1(i,j,k)/(var1(i,j,k)+var2(i,j,k)+var3(i,j,k))
    limT_graz(iorg)%val(i,j,k)=limgraz(i,j,k)
  elseif (VAR(ili)%scomp=="dia") then
    limgraz(i,j,k)=var2(i,j,k)/(var1(i,j,k)+var2(i,j,k)+var3(i,j,k))
    limT_graz(iorg)%val(i,j,k)=limgraz(i,j,k)
  elseif (VAR(ili)%scomp=='detritus') then
    limgraz(i,j,k)=var3(i,j,k)/(var1(i,j,k)+var2(i,j,k)+var3(i,j,k))
  endif
else
   limgraz(i,j,k)=0.
endif
f_lim_graz(i,j,k)=limgraz(i,j,k)


enddo
enddo
enddo

Return
End function f_lim_graz
