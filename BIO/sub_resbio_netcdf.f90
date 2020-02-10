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
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
            SUBROUTINE sub_resbio_netcdf (t_arg)

! Subroutine de sauvegarde des resultats en format netcdf
! 
! A MODIFIER par utilisateur pour zone géographique, niveau sigma, 
!  et date de reference. Cette version est compatible à MARS3D
!
! derniers modifs : 27/06/07
!********************************************************************
!     
!
!     appelant         : main_bio
!                      :
!     appeles          : ionc_init, ionc_creer, ionc_torigine
!     appeles          : ionc_creer8_zxyt, ionc_vatt_missing
!     appeles          : ionc_ecritr8_zxyt
!
!     variables in     : VAR,t_arg
!                      :
!     variables out    :
!                      :
!     variables inout  :                     :
!                      :
!***********************************************************
!     
      Use COUPLEUR_PHY_BIO
      Use PHY_BIO
      Use VAR_USER
      
implicit none

!-- variables locales
      integer ivar,i,j,k,l,compt
      character(20) nom_var_long,nom_var_court
      real*8 , allocatable :: varecrit(:,:,:)

!-- Variables geographiques pour le fichier netcdf             
       real*8 dgac  !pas de la maille en longitude
       real*8 gwac  !borne West de la zone en latitude
       real*8 dfiac !pas de la maille en latitude
       real*8 fisac !borne Sud zone en longitude
       
!-- Variables physiques
       real*8,allocatable :: sig(:)    
       
!-- Parametre temporel
       character*20 dateref
       parameter (dateref='01-Jan-0000 00:00:00')   
       
!-- Parametre de valeur NaN
      real*8 valmanq
      parameter (valmanq = -999.)
!-- Arguments     
      real*8 t_arg
      
!-- Fonctions:
      logical :: isnan


!-- Au premier appel au module bio -> Creation du fichier

  if(nbcallbio.eq.0)then
  
   !Definition des parametres geographiques
   !ex de la caledonie
   dgac = 0.009
   gwac = -166.5
   dfiac = 0.0009
   fisac = -21.8
   
   !Params physiques et temporels
   allocate(sig(nz_max))
   !ATTENTION : DOIT PROVENIR DE LA PHYSIQUE
   !niveau sigma homogene
   if(nz_max.GT.1) then
     sig(1) = 1/sig(k)/2
   else
     sig(1) = 1
   endif
   do k=2,nz_max
     sig(k) = 1/sig(k) + sig(k-1)
   enddo
  
   !Definition du nom du fichier de sortie
    resu_bio = trim(rep_sortie) // 'sortie_' // trim(suffbio) // '.nc'

    write(*,*) 'creation fichier de sortie bio : ',resu_bio
    !-- Création du fichier         
    call ionc_creer(trim(resu_bio),dgac,-gwac,dfiac,fisac,nx_min,&
&                       nx_max,1,ny_min,ny_max,1,1,nz_max,1,nz_max,0,sig)
    !-- Date à l'origine
    call ionc_torigine(trim(resu_bio),dateref) 


    !-- Ajout de la définition de chacune des variables d'état         
    do ivar=1,nbvar      
      nom_var_court = trim(VAR(ivar)%scomp) // '_' // trim(VAR(ivar)%elmt)
      nom_var_long  =trim(VAR(ivar)%comp) // '_' // trim(VAR(ivar)%scomp) // '_' // trim(VAR(ivar)%elmt)
      call ionc_creer8_zxyt(trim(resu_bio),trim(nom_var_court),'micromol/l',trim(nom_var_long),0,0)
      call ionc_vatt_missing(trim(resu_bio),trim(nom_var_court),valmanq)
    enddo   

  else


!-- Ecriture dans le fichier des variables d'etats

  !-- Ajout d'un pas de temps
  call ionc_ecrit_t(trim(resu_bio),0,t_arg)
  allocate(varecrit(nz_max,nx_min:nx_max,ny_min:ny_max))
  varecrit = 0.0
  do ivar = 1, nbvar  
    !-- permutation des indices (i,j,k)->(k,i,j) pour s'adapter avec IONETCDF
    
    do i=nx_min,nx_max
    do j=ny_min,ny_max
    do k=1,nz_max
    !-- Conversion des NaN en -1
      if(isnan(VAR(ivar)%conc(i,j,k))) then
        varecrit(k,i,j) = -1.0
      else
        if(VAR(ivar)%conc(i,j,k).GT.500.0) then
          varecrit(k,i,j) = -2.0
        else
          varecrit(k,i,j) = VAR(ivar)%conc(i,j,k)
        endif
      endif
    enddo
    enddo
    enddo
    !-- ecriture dans le fichier de sortie de la variable ivar   
    nom_var_court = trim(VAR(ivar)%scomp) // '_' // trim(VAR(ivar)%elmt)  
    call ionc_ecritr8_zxyt(trim(resu_bio),trim(nom_var_court), &
&                           varecrit,nx_min,nx_max,ny_min,ny_max,nz_max,0)    
    
  enddo
  deallocate(varecrit)
  
  endif
      
end subroutine sub_resbio_netcdf
