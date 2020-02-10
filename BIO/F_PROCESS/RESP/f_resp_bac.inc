!--------------------------------------------------------------------
!
    function f_resp_bac(ili,jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)
!
!--------------------------------------------------------------------
!
!Cette fonction calcule la respiration bactérienne, en fonction de la production calculée
!
! Variables globales
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL

 Implicit None
!-- fonction
Real(8) :: f_resp_bac(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: qcba(nx_min:nx_max,ny_min:ny_max,nz_max)
Real(8) :: bacp(nx_min:nx_max,ny_min:ny_max,nz_max)
!-- arguments
integer    :: i,j,k,ili,jcol
real(8)    :: qcmin,qnmin,qpmin,QTEN,Trem,mumax,var1,var2,var3

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

if (VAR(ili)%scomp.eq.'at')then
qcba(i,j,k) = VAR(var1)%conc(i,j,k)/nb_bact(i,j)
else
qcba(i,j,k) = VAR(ili)%conc(i,j,k)/nb_bact(i,j)
endif

!Respiration bactérienne
bacp=f_bp_c(ili,jcol,mumax,qcmin,qnmin,qpmin,QTEN,Trem,var1,var2,var3)

br(1)%val(i,j,k)=0.6*qcba(i,j,k)*bacp(i,j,k)+0.01/(60*60*24)*(qcba(i,j,k)-qcmin)*nb_bact(i,j)

f_resp_bac(i,j,k)=br(1)%val(i,j,k)

    enddo
  enddo
enddo
   
 Return
 End function f_resp_bac
