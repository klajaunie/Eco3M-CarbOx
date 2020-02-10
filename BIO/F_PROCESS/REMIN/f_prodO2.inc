!--------------------------------------------------------------------
!!------------ MODIF Lajaunie-Salla 12/06/2019 --------------------!!
!--------------------------------------------------------------------
!
    function f_prodO2(ili,jcol)
!
!--------------------------------------------------------------------
!
! Calcul de la production d'oxygène provenant de la production phyto., de l'uptake en nitrates
! et on soustrait la respiration bactérienne
!
! variables et unites:
!---------------------
!
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use DEF_TYPE

Implicit None
!-- fonction
Real(8) :: f_prodO2(nx_min:nx_max,ny_min:ny_max,nz_max)
!-- arguments
 Integer :: i,j,k,ili,jcol,iphy
!-- variables locales:
real(8) :: mutot(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8) :: prtot(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8) :: upnotot(nx_min:nx_max,ny_min:ny_max,nz_max)
real(8) :: brtot(nx_min:nx_max,ny_min:ny_max,nz_max)

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

 mutot(i,j,k) = 0.0
 prtot(i,j,k) = 0.0
 upnotot(i,j,k) = 0.0
 brtot(i,j,k) = 0.0

 do iphy=1,nscp_phy
   mutot(i,j,k)= mutot(i,j,k)+PPB(iphy)%val(i,j,k)
   upnotot(i,j,k)= upnotot(i,j,k)+hQno(iphy)%val(i,j,k)
   prtot(i,j,k)=prtot(i,j,k)+pr(iphy)%val(i,j,k)
   brtot(i,j,k)=br(1)%val(i,j,k)
 enddo

 f_prodO2(i,j,k)=mutot(i,j,k)+2.*upnotot(i,j,k)-prtot(i,j,k)-brtot(i,j,k)


! Remarque importante: vu qu'en rentrant les variables dans "config.ini",
! on commence toujours par la biomasse, le nutriment absorbe
! est toujours dans la colonne jcol du flux (ili --> jcol), et la
! biomasse est toujours dans la ligne ili

    enddo
  enddo
enddo

Return
end function f_prodO2
