!----------------------------------------------------------------
!
                   function f_precipcarbo(ili,jcol,kprecip,Kc)
!----------------------------------------------------------------
!
! Creation Katixa Lajaunie-Salla 11/06/2018
! Calcul du flux de dossolution des carbonates: Ca2+ + CO32- ==> CaCO3
!
Use COUPLEUR_PHY_BIO
Use  VAR_GLOBAL

Implicit None

REAL(8):: f_precipcarbo(nx_min:nx_max,ny_min:ny_max,nz_max)
REAL(8):: precip(nx_min:nx_max,ny_min:ny_max,nz_max)
Integer :: i,j,k,ili,jcol,iphy
Real(8) :: kprecip, Kc !fraction PIC/POCmax & Michaelis-Menten constant
 real(8) :: npptot

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

 npptot = 0.0

do iphy=1,nscp_phy
  npptot= npptot+(PPB(iphy)%val(i,j,k)-pr(iphy)%val(i,j,k))
enddo

if (Omega(i,j,k)-1>0.) then
  precip(i,j,k)=kprecip*((Omega(i,j,k)-1)/(Kc+Omega(i,j,k)-1))
  if (npptot<0.) npptot=0.
else
 precip(i,j,k)=0.
endif

! une fraction de la nette production primaire (en molC/m3/s)
 f_precipcarbo(i,j,k)=precip(i,j,k)*npptot


    enddo
  enddo
enddo

Return
end function f_precipcarbo
