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
!--------------------------------------------------------------------
!
       function f_extinc_morel(kk0,aa0,Chl_tot,dz)
!    
!-------------------------------------------------------------------------
! Cette fonction calcule le terme d'extinction lumineuse
! dans la colonne d'eau par le modele de Morel
! 
! Variables et parametres et leurs unites:
!-----------------------------------------
!  kk0  : coefficient d''extinction de l''eau de mer pure (m-1), val. typique 0.0384
!  aa0  : coefficient d''absorption de  Chla (m2*mg-1), val. typique 0.029

! Ref:  Morel et al. (1993) J. MAR. Res., 51:617-
!
! Melika Baklouti 15/06/2007
!-------------------------------------------------------------------------
!
! USE COUPLEUR_PHY_BIO
! USE VAR_GLOBAL
! USE MOD_FCHAIN

 Implicit None

!-- fonction:
 Real(8):: f_extinc_morel

!-- arguments
 real(8) :: kk0  !(0.0384)  coefficient d''extinction de l''eau de mer pure (m-1)
 real(8) :: aa0  !(0.029)  coefficient d''absorption de  Chla (m2*mg-1)
 real(8),optional :: Chl_tot,dz

 
!-- relation pour Chl exprime en mg m-3
   f_extinc_morel = exp(-dz*(kk0 + 0.0088 * Chl_tot*1.d3  + 0.054 * (Chl_tot*1.d3)**(2./3))) 
 
end function f_extinc_morel
!--------------------------------------------------------------------
