!--------------------------------------------------------------------
! Dernière modification: 18/03/2005
!
! Cette fonction calcule le facteur adimensionnel (fonction
! des quotats) qui limite l'UPTAKE du fait du nutriment X
! dont la conc. interne dans le phyto concerne est donnee par 
! l''indice indX
! ref: article de Geider, 1998
! NB: cette fonction est plutot adaptee au Nitrate mais pas
!     trop aux autres nutriments (cf article Geider)
!
       function f_fQupt_geid(ili,jcol,indX,Qmin,Qmax,n)
!
!--------------------------------------------------------------------
 USE COUPLEUR_PHY_BIO
 USE VAR_GLOBAL
 USE MOD_FCHAIN

 Implicit None

!-- fonction:
 Real(8):: f_fQupt_geid(nx,ny,nz)

!-- arguments
 Real(8) :: Qmin, Qmax,indX,n
 Integer :: ili,jcol,iindX

!--variables locales:
 Integer :: iiC
 Real(8) :: Q(nx,ny,nz)
 Integer:: identite
 Integer:: iphyC
 Character(L_VAR)::chain

 identite = var(ili)%idorg
 chain = 'C'
 iphyC = f_idorg2id(identite,chain)
 iindX = nint(indX)
 
!-- Calcul du quota cellulaire:
Q = var(iindX)%conc / (var(iphyC)%conc)! quota cellulaire X:C

if (maxval(Q ) >= Qmax) then
  !write(*,*) 'Q > Qmax dans f_fQupt_geid à t(h) =', tps/3600.
  f_fQupt_geid =0.d0
else
  f_fQupt_geid =  ((Qmax-Q) / (Qmax - Qmin))**n
endif


end function f_fQupt_geid
!--------------------------------------------------------------------
