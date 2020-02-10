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
!------------------------------------------------------------------------
! 
                  SUBROUTINE MAIN_BIO
!
!
! Programme Eco3M (Ecological Mechanistic and Modular Modelling)
!
! cette entité de programmation est le programme principal de la 
! modélisation biogéochimique 	

!-- Ce programme est destine à la modelisation des dynamiques des
!   premiers echelons du reseau trophique à l aide de 
!   modèles de processus biogeochimiques. Il peut etre utilise sous 2 modes:
!
!   1/ en mode INITIALISATION :  ce mode sert à l implementation du  
!      modèle choisi par l utilisateur et à la creation d''un fichier
!      call.inc qui regroupe l appel aux fonctions et routines pour
!      le calcul des flux de masse entre variables d etat
!
!   2/ en mode CALCUL : ce mode ne peut être active que si le mode
!      precedent a ete lance au moins une fois au prealable, de façon à 
!      ce que le fichier call.inc ait ete cree. Ce mode permet la 
!      simulation dynamique des processus mis en jeu (et notamment 
!      le calcul des termes sources et puits pour chaque variable d''etat
!      afin de les communiquer au modèle physique !
!
!
! Concepteurs :
!--------------
!
! Melika BAKLOUTI : Laboratoire d Oceanographie et de Biogeochimie (LOB-UMR 6535)
!                    détachée à l Institut de Recherche pour le Developpement (IRD) entre 2002-2006
!                   
! Vincent FAURE :   Université de la Méditerranée  
!                   Laboratoire d Océanographie et de Biogéochimie (LOB-UMR 6535)
!                   Institut de Recherche pour le Développement (IRD)
!
!-- Dernière modification: 27/06/2007
!------------------------------------------------------------------------
!
! Variables globales
 Use  DEF_TYPE
 Use  COUPLEUR_PHY_BIO
 Use  VAR_USER
 Use  VAR_GLOBAL
 Use  F_PROCESS

 Implicit None

!-- Variables locales
 integer,PARAMETER :: nb_call_exec = 1 !-- nb d appels successifs de la subroutine de calcul des flux
 integer :: i,j,k,ili,jcol,it,ii
 integer :: ivar,jvar,kvar
 integer :: numfich
 integer :: choix
 logical :: present
 character(60) nom_fic1,nom_fic2,nom_fic3,nom_fic4

!----------------------------------------------------------------------
!-- Lecture des donnees en vue de la creation du fichier call.inc  ----
!----------------------------------------------------------------------
 
!       ----------------------------------------------
!       -------     I.  MODE INITIALISATION  ---------
!       ----------------------------------------------
#ifdef INI
      write(*,*) 'MODE BIO  =  INITIALISATION '
!--Initialisation des processus du modèle:
      Call Sub_init_mod

!--Initialisation du profil des Variables:
      Call Sub_init

!-- Creation du fichier call.inc:
      Call Sub_call_tab

!#ifndef NCOUPL-M0D 
#ifndef M0D-NCOUPL 
!-- Creation du fichier calc_extinc.inc
      Call Sub_calc_extinc
#endif

!-- Creation du fichier save_flux.inc:
if (allocated(coord_flux)) then
      Call Sub_save_flux     
endif
!
      write(*,*) '--------- FIN CREATION DES APPELS --------------'
      write(*,*) 
#endif

!       ----------------------------------------------
!       -------      II.  MODE CALCUL        ---------
!       ----------------------------------------------
#ifdef CALC
!-- Premier appel de main_bio en mode calcul
  if (nbcallbio ==-1) then
 !--Initialisation des processus du modèle:
      Call Sub_init_mod

!--Initialisation du profil des Variables:
      Call Sub_init
      
#ifdef COUPL
!-- creation de la liste de fichiers pour       
!-- lire les conc initiales ds le programme phy
      Call Sub_var_ini
#endif        

!--Creation des noms des fichiers de sorties
      nom_fic1 = trim(rep_sortie)//'CR_'//trim(suffbio)//'.sim'
      nom_fic2 = trim(rep_sortie)//'bilan_'//trim(suffbio)//'.sim'
      nom_fic3 = trim(rep_sortie)//'var_phy_'//trim(suffbio)//'.sim'

!-- Reouverture des fichiers 
      Open(22,FILE=nom_fic1) 
      Open(23,FILE=nom_fic2)
      Open(24,FILE=nom_fic3)
!   

      Call sub_affiche(1) ! Ecriture du fichier compte-rendu (fichier 22)
      Call Sub_affiche(2) ! Ecriture des entetes fichiers 23,24,25

!--------------------------------------------------------------------------
  elseif (nbcallbio > -1) then 
!--------------------------------------------------------------------------      
#ifdef MODTEST
   !  write(*,*) '------ debut main_bio (cas general) ------'
#endif

!-- Second appel de main_bio en mode calcul
  if (nbcallbio ==0) then

#ifdef NCOUPL      
!-- Allocation et initialisation des champs de concentrations:
      Call Sub_var_ini
#endif     
!-- Calcul des rapports Chl_C :
     if (CHL_C_BOOL) then
         Call sub_calc_Chl_C
      endif 
!--Initialisation et allocation des variables physiques
!      Call sub_phys_ini  !inutile a priori et supprimé ds la prochaine version

!-- Creation du fichier de resultat et sauvegarde de l''instant initial
    #ifdef RESBIO_A
       call sub_resbio_ascii(tps)
    #endif
    #ifdef RESBIO_N
       call sub_resbio_netcdf(tps)
    #endif
 #ifdef MODTEST
  !  write(*,*) 'retour de sub_resbio'
 #endif       
  endif ! fin des instructions specifiques au second appel de main_bio


!--Gestion de l irradiance
!--dt_irrad : pas de temps (converti de min en sec dans sub_init) pour lequel on a une 
!  nouvelle valeur d irradiance 
!  tps = tps en s      

!--lecture ou calcul de l irradiance (ou du PAR) en surface (sauf si elle est transmise
!  par le code physique):
if (fichirrad /= 'IRR_CODEPHYS') then
   if (nbcallbio== 0 .or. mod(tps,dt_irrad)==0) then 
   write(*,*)'pas de tps irrad',tps,dt_irrad,mod(tps,dt_irrad)
        call sub_lec_irrad 
!--Calcul de l irradiance sur la verticale (1DV ou 3D)
       call sub_calc_EPARZ
    endif
endif
!stop
!-- fermeture du fichier de compte-rendu
if (nbcallbio==0) close(22)

!--------------------------------------------------------------------------
!----------- calcul des flux biogeochimiques     --------------------------
!--------------------------------------------------------------------------

!-- Appel des fonctions calculant les flux
!marion
!write(*,*)'nb_bact',nb_bact(:,:)
do i=1,nb_call_exec
! DEBUG Christel
write(*,*) 'nb call',i
     CALL SUB_EXEC_CALL
enddo
!
#ifdef MODTEST
    write(*,*) '------ Calcul des flux termine ------'
#endif

!-- Sauvegarde des flux
if (nbcallbio==0 .or. mod(tps,dt_save_bio*60)==0) then
     if (allocated (noms_fich_flux))then
       Call Sub_affiche(4) ! affichage des flux en un point donne
                                  ! de coordonnees (i,j,k) dans le
     endif                             ! fichier numfich (numfich >= 30)
endif
!debug marion
!write(*,*)'Phyto C, N , P et Chl',VAR(1)%conc(i,j,k),VAR(2)%conc(i,j,k),VAR(3)%conc(i,j,k),VAR(4)%conc(i,j,k) 
!write(*,*)'bac C, N et P',VAR(5)%conc(i,j,k),VAR(6)%conc(i,j,k),VAR(7)%conc(i,j,k)



!--------------------------------------------------------------------------
!----------- calcul des tendances sur le maillage --------------------------
!--------------------------------------------------------------------------
!
#ifdef MODTEST
        write(*,*) 'Calcul des Tendances.....'
#endif
        TEND = 0.d0
!-- Somme des termes de tendance (sources - puits) pour chaque variable
     do ivar=1,nbvar
!-- termes puits/sources ne faisant intervenir qu une variable d etat
        if (associated (SELF_VAL(ivar)%idproc)) then
!            if (ivar==18) write(*,*)'Tend DIC self',ivar,TEND(ivar,:,:,:)
           TEND(ivar,:,:,:) = SELF_VAL(ivar)%val(:,:,:)
!            if (ivar==18) write(*,*)'Tend DIC self',ivar,TEND(ivar,:,:,:),SELF_VAL(ivar)%val(:,:,:)
        endif
!-- Rappel: FLUX_VAL(i,j) est un flux de i vers j : 
!           donc c est un puit pour i et une source pour j
       if (ivar < nbvar ) then
        do jvar = ivar+1,nbvar
          if (associated (FLUX_VAL(ivar,jvar)%idproc)) then
!            if (ivar==18.or.jvar==18) write(*,*)'Tend DIC flux',ivar,TEND(ivar,:,:,:)
           TEND(ivar,:,:,:) = TEND(ivar,:,:,:) -  FLUX_VAL(ivar,jvar)%val(:,:,:) 
!            if (ivar==18.or.jvar==18) write(*,*)'Tend DIC flux',ivar,jvar,TEND(ivar,:,:,:),FLUX_VAL(ivar,jvar)%val(:,:,:)	
          endif
        enddo
       endif
!-- On prend aussi en compte les termes sources pour la variable numero ivar
      if (ivar > 1) then 
         do kvar = 1,ivar-1
          if (associated (FLUX_VAL(kvar,ivar)%idproc)) then
!              if (ivar==18.or.kvar==18) write(*,*)'Tend DIC flux source',ivar,TEND(ivar,:,:,:)
            TEND(ivar,:,:,:) = TEND(ivar,:,:,:) + FLUX_VAL(kvar,ivar)%val(:,:,:)
!              if (ivar==18.or.kvar==18) write(*,*)'Tend DIC flux sources',ivar,kvar,TEND(ivar,:,:,:),FLUX_VAL(kvar,ivar)%val(:,:,:)	
          endif
         enddo
      endif
    ! write(*,*) 'ii,max(TEND)',ivar,maxval(TEND(ivar,:,:,:))
     enddo
    
#ifdef MODTEST
!     write(*,*)'TEMP_BIO main271=',TEMP_BIO
     write(*,*) '------ Calcul des tendances terminé ------'
#endif

#ifdef NCOUPL
!-- Integration methode d''Euler explicite (à améliorer !):
     do ivar=1,nbvar !var_new=var_old+(dvar/dt)*dt
       if (ivar==19) write(*,*)'Var at',ivar,TEND(ivar,:,:,:),VAR(ivar)%conc(i,j,k),dt_bio
        do i=nx_min,nx_max
        do j=ny_min,ny_max
        do k=1,nz_max
         VAR(ivar)%conc(i,j,k) = max (1.d-50 ,VAR(ivar)%conc(i,j,k) + TEND(ivar,i,j,k) * dt_bio) 
!            if (SAL_BIO(i,j,k)<=37.0d0) then
!            if (ivar==25) then 
!            VAR(ivar)%conc(i,j,k) = capp(ivar,i,j,k) 
!            VAR(25)%conc(i,j,k) = (2.5*SAL_BIO(i,j,k)+2500)*rho_sw(i,j,k)*1.d-3 ! µmol/L eq. avec toutes les valeurs SOMLIT
!            VAR(25)%conc(i,j,k) = (-21*SAL_BIO(i,j,k)+3400)*rho_sw(i,j,k)*1.d-3 ! µmol/L eq. avec les valeurs seulement d intrusion
!            VAR(25)%conc(i,j,k) = 2600.*rho_sw(i,j,k)*1.d-3 ! µmol/L
!            endif
!            endif
        enddo
        enddo
        enddo
!       open (145, FILE='./BIO/SORTIES/tend.txt')
!       write(145,*)ivar,TEND(ivar,1,1,1),VAR(ivar)%conc(1,1,1)
     enddo
#endif
!
!
if (nbcallbio==0 .or. mod(tps,dt_save_bio*60)==0) then
   Call Sub_affiche(3) ! Bilan de matiere, conservativite 
! 


   Call Sub_affiche(5) ! affichage des variables physiques
endif
!
#ifdef MODTEST
 write(*,*) '---------- FIN CALCUL BIO --------------'
#endif

!-- Calcul du pCO2 et pH a partir de DIC et AT
 CALL sub_calc_CO2
 write(*,*)'calcul systeme des carbonates'

!Ecriture dans le fichier de resultat
if (nbcallbio > 0 .and. mod(tps,dt_save_bio*60)==0) then !dt_save_bio en mn
#ifdef MODTEST
 write(*,*) '---------- DEBUT SAUVEGARDE BIO --------------'
#endif

  

 #ifdef RESBIO_A
      call sub_resbio_ascii(tps)
 #endif
 #ifdef RESBIO_N
     call sub_resbio_netcdf(tps)
 #endif
#ifdef MODTEST
 write(*,*) '---------- FIN SAUVEGARDE BIO --------------'
#endif
endif

#ifdef NCOUPL
if (tps>duree*journee) then
      write(*,*) 'temps final', duree,'(j) atteint'
      Call sub_fin
      stop
endif
#endif
endif !-- fin du if nbcallbio >-1
!---------------------------------------------------------------------
!           Nombre d''appels du module bio (debute a -1)
!---------------------------------------------------------------------
nbcallbio = nbcallbio + 1 


#endif !-- va avec le ifdef CALC

End SUBROUTINE MAIN_BIO
!---------------------------------------------------------------------
      SUBROUTINE SUB_EXEC_CALL
!--------------------------------------------------------------------
! Variables globales
 Use  DEF_TYPE
 Use  COUPLEUR_PHY_BIO
 Use  VAR_USER
 Use  VAR_GLOBAL
 Use  F_PROCESS

Implicit None
!
#ifdef MODTEST
 write(*,*) 'Debut de SUB_EXEC_CALL'
#endif


#ifdef CALC
  include "call.inc"
#endif 
write(*,*) 'END SUBROUTINE SUB_EXEC_CALL'

#ifdef MODTEST
 write(*,*) 'Fin de SUB_EXEC_CALL'
#endif
END SUBROUTINE SUB_EXEC_CALL
!------------------------------------------------------------------
