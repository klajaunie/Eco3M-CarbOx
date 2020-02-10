!--------------------------------------------------------------------
!
   function f_resp(ili,jcol,r)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule la respiration, fraction de la biomasse carbonée
!
!
! Paramètres :
!  r : respiration par rapport à la biomasse

USE VAR_GLOBAL
USE COUPLEUR_PHY_BIO
Implicit None
!-- fonction
Real(8) :: f_resp(nx_min:nx_max,ny_min:ny_max,nz_max)
!-- arguments
Real(8) :: r
Integer :: i,j,k,ili,jcol,iorg

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

iorg=VAR(ili)%idorg

!13/04/06 respiration sur la production brute
 f_resp(i,j,k)=r*VAR(ili)%conc(i,j,k)

!Enregistrement de la respiration
pr(iorg)%val(i,j,k) = f_resp(i,j,k)
     

    enddo
  enddo
enddo

Return
end function f_resp
