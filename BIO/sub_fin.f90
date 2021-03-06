!***************************************************************************
!***************************************************************************
!Copyright or � or Copr. CNRS/IRD/Universit� de la M�diterran�e
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
               Subroutine Sub_fin_INI
!
!-- Routine de fin de l'execution du mod�le en mode INI
!
!   * lib�ration de la m�moire de tous les tableaux dynamiques
!
! Derniere modif: 05/12/2005
!
!---------------------------------------------------------------------
Use DEF_TYPE
Use VAR_USER
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL

Implicit None

! Variables locales
 
 Integer::i,j,k
 Integer :: ivar

!-- Desallocation des tableaux dynamiques:
deallocate(VAR)
if (allocated(nb_elmt)) deallocate(nb_elmt)
if (allocated(tab_noms_fich)) deallocate(tab_noms_fich)
!----lib�ration de la m�moire
 
 do  i=1,nbvar 
   do j=1,nbvar
    if(associated(FLUX_VAL(i,j)%idproc)) nullify(FLUX_VAL(i,j)%idproc)
   enddo
    if(associated(SELF_VAL(i)%idproc)) nullify(SELF_VAL(i)%idproc)
  enddo

 do i=1,nbproc_bib
   if(associated(PROC_MOD(i)%nompar)) nullify(PROC_MOD(i)%nompar)
 enddo

 do i=1,nbproc_flux
   if(associated(FLUX_PAR(i)%valpar)) nullify(FLUX_PAR(i)%valpar)
 enddo
 do i=1,nbproc_self
   if(associated(SELF_PAR(i)%valpar)) nullify(SELF_PAR(i)%valpar)
 enddo

   if(associated(IRR_PAR%valpar)) nullify(IRR_PAR%valpar)

if (allocated(FLUX_VAL)) deallocate(FLUX_VAL)
if (allocated(SELF_VAL)) deallocate(SELF_VAL)

if (allocated(FLUX_PAR)) deallocate(FLUX_PAR)
if (allocated(SELF_PAR)) deallocate(SELF_PAR)

if (allocated(PROC_MOD)) deallocate(PROC_MOD)

if (allocated(CHL_C)) deallocate(CHL_C)


 if (allocated (iscp_phy)) deallocate(iscp_phy)
 if (allocated (iscp_zoo)) deallocate(iscp_zoo)
 if (allocated (iscp_bac)) deallocate(iscp_bac)

if (allocated (coord_flux)) deallocate (coord_flux)
if (allocated (noms_fich_flux)) deallocate (noms_fich_flux)
 
 write(*,*) " -- Sub_fin_INI Termin� ! -- "

End Subroutine Sub_fin_INI

!-------------------------------------------------------------------------------
!
               Subroutine Sub_fin
!	       
!-- Routine de fin du mod�le
!
!   * lib�ration de la m�moire de tous les tableaux dynamiques
!
!---------------------------------------------------------------------
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Use VAR_USER

 Implicit None

! Variables locales
 
 Integer::i,j,k,istat
 Integer :: ivar


!----lib�ration de la m�moire
 
 do  i=1,nbvar 
   do j=1,nbvar
     if(associated(FLUX_VAL(i,j)%val)) nullify (FLUX_VAL(i,j)%val)
     if(associated(FLUX_VAL(i,j)%idproc)) nullify(FLUX_VAL(i,j)%idproc)
   enddo
    if(associated(SELF_VAL(i)%val)) nullify(SELF_VAL(i)%val)
    if(associated(SELF_VAL(i)%idproc)) nullify(SELF_VAL(i)%idproc)
    if(associated(VAR(i)%conc)) nullify(VAR(i)%conc)
  enddo

 do i=1,size(PROC_MOD)
   if(associated(PROC_MOD(i)%nompar)) nullify(PROC_MOD(i)%nompar)
 enddo

do i=1,size(FLUX_PAR)
   if(associated(FLUX_PAR(i)%valpar)) nullify(FLUX_PAR(i)%valpar)
enddo

if (allocated(VAR)) then
  deallocate(VAR)
endif
if (allocated(FLUX_VAL)) then
 deallocate(FLUX_VAL)
endif
if (allocated(SELF_VAL)) then
 deallocate(SELF_VAL)
endif
if (allocated(FLUX_PAR)) then
 deallocate(FLUX_PAR)
endif
if (allocated(PROC_MOD)) then
 deallocate(PROC_MOD)
endif
if (allocated(mu_phy)) then
 deallocate(mu_phy)
endif
! deallocate(IRRAD)
if (allocated(E_PARZ)) then
 deallocate(E_PARZ)
endif
if (allocated(TEMP_BIO)) then
 deallocate(TEMP_BIO)
endif
if (allocated(SAL_BIO)) then
 deallocate(SAL_BIO, STAT=istat)
if (istat /=0) write(*,*) 'pb de desallocation de SAL_BIO'
endif
if (allocated(WIND)) then
 deallocate(WIND)
endif
if (allocated(CO2atm)) then
 deallocate(CO2atm)
endif
if (allocated(rho_sw)) then
 deallocate(rho_sw)
endif
!if (allocated(capp)) then
! deallocate(capp)
!endif
if (allocated(CHL_C)) then
 deallocate(CHL_C)
endif
if (allocated(coord_flux)) then
 deallocate(coord_flux, STAT = istat)
if (istat /=0) write(*,*) 'probleme de desallocation de coord_flux'
endif

   print * , " -- Sub_fin Termin� ! -- "


End Subroutine Sub_fin
