!----------------------------------------------------------------
!
                   function f_disscarbo(ili,jcol,kcarb)
!----------------------------------------------------------------
!
! Creation Katixa Lajaunie-Salla 11/06/2018
! Calcul du flux de dossolution des carbonates: CaCO3 ==> Ca2+ + CO32- 
!
Use COUPLEUR_PHY_BIO
Use  VAR_GLOBAL

Implicit None

REAL(8):: f_disscarbo(nx_min:nx_max,ny_min:ny_max,nz_max)
REAL(8):: diss(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: kcarb,varCaCO3
Integer :: i,j,k,ili,jcol,id_diaC,id_picoC
character(L_CHAIN) :: scomp,elmt

do i=nx_min,nx_max
do j=ny_min, ny_max
do k=1,nz_max

if (Omega(i,j,k)-1>0.) then
 diss(i,j,k)=0.
else
 diss(i,j,k)=(1-Omega(i,j,k))*kcarb
endif

if(VAR(ili)%scomp=='dic')then
  varCaCO3=VAR(ili+4)%conc(i,j,k)
elseif(VAR(ili)%scomp=='at')then
  varCaCO3=VAR(ili+1)%conc(i,j,k)
endif

f_disscarbo(i,j,k)=diss(i,j,k)*varCaCO3


enddo
enddo
enddo

Return
end function f_disscarbo
