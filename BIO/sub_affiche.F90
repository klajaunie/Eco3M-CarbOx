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
! Dernière modification: 26/01/05
!
!-- Routine d''affichage dans le fichier compte-rendu 22, et le fichier de resulat 23!
!
!  choix=1 -> Remplissage du fichier de compte rendu (n 22)
!  choix=2 -> en tete des autres fichiers de resultats 
!  choix=3 -> Bilans de matiere par element (pour étudier la conservativite ou au contraire l''export)
!  choix=4 -> calcul des flux (vers un numero de fichier et pour une station donnée)
!  choix=5 -> variables physiques
!
! numfich : numero de l''unité du fichier (0=ecran) dans lequel on veut écrire les flux (choix=4)
! i_st,j_st,k_st : station pour laquelle on veut écrire l''analyse (choix=4)
!-------------------------------------------------------------------------------

      Subroutine Sub_affiche(choix)                 

!-------------------------------------------------------------------------------						
Use DEF_TYPE
Use VAR_USER
Use VAR_GLOBAL
Use COUPLEUR_PHY_BIO
Use MOD_FCHAIN
Use  F_PROCESS

Implicit None
 
! Variables locales
!------------------- 
 integer          :: i,j,k,ili,jcol,imat,i_st,j_st,k_st,numfich
 integer          :: iii,dimfonc
 real(8)          :: tempor
 character(10)    :: valeur
 character(80)    :: numfichtxt
 character(40)    :: nomfic,namo
 Character(200)   :: chaine
 character(15)    :: chaine_petite
 Character(15)    :: nom_param
 Character(50)    :: fich
 Character(3)     :: debfichirrad
 integer          :: choix,choix_calcul
 real             :: som,somC, somN 
 real(8),Allocatable:: ss_flux(:,:,:,:)


select case (choix)

!----------------------------------------------------------------------
  case(1) !-- Ecriture des caracteristiques du modele et de la simulation (fichier CR)
!----------------------------------------------------------------------	
!	 
      write(22,*) '------------- COMPTE RENDU DE SIMULATION -----------'
      write(22,*) '----------------------------------------------------'
      write(22,*) '---------------------- INITIALISATION --------------'
      write(22,*) ' '
!-- Ecriture de la matrice PROC_MOD dans le fichier de compte-rendu
      write(22,*) " "
      write(22,*) "---PROC_MOD"
      do i=1,size(PROC_MOD)
        write(22,'(I3,$)') PROC_MOD(i)%idproc
        write(22,*)' : ',trim(adjustl(PROC_MOD(i)%nomproc)),' : ',&
                       trim(adjustl(PROC_MOD(i)%nomsub)),       &
                       PROC_MOD(i)%nbpar,' :',                 &
                       (trim(adjustl(PROC_MOD(i)%nompar(k))),'|', k=1,size(PROC_MOD(i)%nompar))
      enddo 
      write(22,*) '-------------------------------------------------------------------------'   
      write(22,'(a40,F5.1,a25,F7.1)') ' paramètres temporels : duree simul (j)=',duree,' - pas de temps bio(s)=',dt_bio
      write(22,*) 'Dimensions spatiales : ','nx,ny,nz',nx_min,':',nx_max,'/',ny_min,':',ny_max,'/',nz_max
      Write(22,*) 'Nombre de variables : ','nbcomp,nbscomp,nbvar',nbcomp,nbscomp,nbvar
      write(22,*) '-------------------------------------------------------------------------'

! -- Affichage de la matrice VAR
      write(22,*) "---VAR"
      do imat=1,nbvar
        write(22,*) imat,' - ',VAR(imat)%comp,' - ',VAR(imat)%scomp,' - ',&
       VAR(imat)%elmt
      enddo
      write(22,*) '-----------------'

!-- Affichage de la matrice FLUX

      write(22,*) "---MATRICE DES FLUX"
      do ili=1,nbvar
        do jcol=1,nbvar
          if(ASSOCIATED(FLUX_VAL(ili,jcol)%idproc)) then 
           write(22, *) ili,';',jcol,':',(FLUX_VAL(ili,jcol)%idproc(k),'/',&
                           k=1,size(FLUX_VAL(ili,jcol)%idproc))
          endif
        end do
     end do
     write(22,*) "-----------------"  

!-- Affichage de la matrice SELF

      write(22,*) "---MATRICE DES PUITS et SOURCES"
      do ili=1,nbvar
       if(ASSOCIATED(SELF_VAL(ili)%idproc)) then   
         write(22, *) ili,':',(SELF_VAL(ili)%idproc(k),'/',&
                             k=1,size(SELF_VAL(ili)%idproc))
       endif
      end do
!
      write(22,*) "-----------------"   

!-- Affichage des parametres lies aux processus du modele et
!   impliques dans des flux entre 2 variables d''etat:

      write(22,*)  '-----PARAMETRES des PROCESSUS du MODELE-----'
      write(22,*)  '      (Flux entre 2 variables d''etat)      '
      write(22,*) ''
      do i=1,size(FLUX_PAR)
        write(22,'(I3,I3,A3,A1,A10,A3,$)')  FLUX_PAR(i)%ipos ,' > ',                 &
          FLUX_PAR(i)%signe,PROC_MOD(FLUX_PAR(i)%idproc)%nomproc , ' : '              
        write(22,'(A7,A1,G10.3,A1,$)') (PROC_MOD(FLUX_PAR(i)%idproc)%nompar(k),'=',   &
       FLUX_PAR(i)%valpar(k),';',k=1,size(FLUX_PAR(i)%valpar))
        write(22,*) 
      enddo
      write(22,*) "-----------------"   

!-- Affichage des parametres lies aux processus du modele et
!   impliques par une seule variable d''etat::

      write(22,*) '-----PARAMETRES des PROCESSUS du MODELE-----'
      write(22,*) '       (Flux avec 1 variable d''etat)       '
      write(22,*) ''
      do i=1,size(SELF_PAR)
        write(22,'(I3,A3,A1,A10,A3,$)') SELF_PAR(i)%ipos(1) ,' > ', SELF_PAR(i)%signe, &
          PROC_MOD(SELF_PAR(i)%idproc)%nomproc , ' : '                             
        write(22,'(A7,A1,G10.3,A1,$)') (PROC_MOD(SELF_PAR(i)%idproc)%nompar(k),'=',     &
                SELF_PAR(i)%valpar(k),';',k=1,size(SELF_PAR(i)%valpar))
        write(22,*) 
      enddo

!                   --------------------------------
!                   --  Gestion de l''irradiance:  --
!                   --------------------------------

        fichirrad = adjustl(fichirrad)
        debfichirrad= fichirrad(1:3)        
 if (fichirrad /= 'IRR_FONCTION' .and.  fichirrad /= 'IRR_CODEPHYS') then
      if (debfichirrad == 'IRR' .or. debfichirrad == 'Irr'.or. debfichirrad == 'irr') then
	write(22,*) '-------------------- IRRADIANCES ----------------------'
        write(22,*) 'fichier d''irradiances E(0+):',  fichirrad
      elseif(debfichirrad == 'PAR' .or. debfichirrad == 'Par'.or. debfichirrad == 'par') then
        write(22,*) '--------------  PAR (Photosynthetic Active Radiation) -------------------'
        write(22,*) 'fichier donnant les valeurs de PAR E_PAR(0+):',  fichirrad
      endif
 elseif (fichirrad == 'IRR_FONCTION') then
	write(22,*) '-------------------- IRRADIANCES ----------------------'
        write(22,*) 'Irradiance E(0+) calculee par une fonction'   
 elseif (fichirrad == 'IRR_CODEPHYS') then
	write(22,*) '-------------------- IRRADIANCES ----------------------'
        write(22,*) 'Irradiance E(0+) fournie par le code physique'   
 endif
 write(22,*) '-----------------  FIN DU COMPTE-RENDU  --------------------------'
! 
!----------------------------------------------------------------------
  case(2) !-- Entêtes des autres fichiers de sauvegarde
!----------------------------------------------------------------------		 
       
!-- Entete du fichier calculant la biomasse par element (C,N,P,Si,...)  
      write(23,*) '%----------- COMPTE RENDU DE SIMULATION ------------'
      write(23,*) '%---------------------------------------------------'
      write(23,*) '%----- CALCUL  DE lA BIOMASSE PAR ELEMENT ----------'
      write(23,*) '%-----   (en mol/m3 sauf Chl en g/m3)      ---------'
      write(23,*) '% '
!     write(23,*) '%--------------------------------  VAR (temps initial) '
!   write(23,'(a11,$)') '% time(h)||'
!   do i=1,nbvar
!       write(23,'($a5,$a1,$a3,$a2)') VAR(i)%scomp, '/', VAR(i)%elmt, '||'
!   enddo
!   write(23,'(/)')
!-- Ecriture des variables d''état en (1,1,1) au temps initial:
!   write(23,'(G10.5,$)')  tps/3600.d0
!   write(23,'($1X,f9.4)') (VAR(i)%conc(1,1,1),i=1 , nbvar)
   
!   write(23,*) 
      write(23,*) 
      write(23,*) 'temps (h) || TOTAL ||  C  ||  N  ||  P || Si ||  Chl  ||  Fe  ||'
      write(23,*) 
   
!-- En tete du fichier donnant l''irradiance au cours du temps 
      write(24,'(a18)') '% time(h) || EPAR || TEMP || SAL || WIND || CO2atm  '
!
!----------------------------------------------------------------------
  case(3) ! Bilan de matiere par element
!----------------------------------------------------------------------		 
!
!-- Conservativité de la matiere totale
   som=0.d0
   do i=1,nbvar
      som=som+VAR(i)%conc(1,1,1)
   enddo
   
!-- recherche de toute la biomasse carbonée

somC=0.d0

do i=1,nbvar
   if (var(i)%elmt=='C') then 
      somC = somC + VAR(i)%conc(1,1,1)
   endif
enddo
!-- recherche de toute la biomasse azotée
somN=0.d0

do i=1,nbvar
   if (var(i)%elmt=='N') then 
      somN = somN + VAR(i)%conc(1,1,1)
   endif
enddo
   
write(23,'(g12.3,1X,g12.3,1X,g12.3,1X,g12.3)') tps/3600.d0,som, somC,somN

!-------------------------------------------------------------------  
  case(4) ! Affichage de la matrice flux en un point donne de coordonnees (i_st,j_st,k_st)
          !  dans le fichier de numero logique numfich
!-------------------------------------------------------------------  
 if (nbcallbio ==0) then
  do iii=1,size(noms_fich_flux)
    fich = trim(adjustl(rep_sortie))// trim(adjustl(noms_fich_flux(iii)))//'.sim'
    numfich = 200+coord_flux(iii,1)+coord_flux(iii,2)+coord_flux(iii,3)
    open(numfich, File = fich)
    write(numfich,*) '# Fichier sauvant les flux entre variables au point d''espace'
    write(numfich,*) '# de coordonnées selon x, y et z respectivement:'
    write(numfich,*) '#',coord_flux(iii,1),coord_flux(iii,2),coord_flux(iii,3)
    write(numfich,*) '#'
    write(numfich,*) '# la 1ere colonne du fichier donne le temps en heure'
    write(numfich,*) '# les colonnes 2 et 3 donnent les indices des variables entre lesquelles'
    write(numfich,*) '# il existe un flux non nul'
    write(numfich,*) '# les autres colonnes donnent les valeurs des différentes FONCTIONS mises '
    write(numfich,*) '# en jeu par ce flux'
    write(numfich,*) '# la dernière colonne donne la valeur totale du flux entre ces deux variables'
    write(numfich,*) '# '
    write(numfich,*) '# ATTENTION: la dernière colonne n''est pas nécessairement la somme algébrique '
    write(numfich,*) '# des colonnes précédentes car les fonctions peuvent intervenir sous forme de '
    write(numfich,*) '# produit et pas uniquement sous forme de sommes algébriques '
    write(numfich,*) '# '
  enddo
 endif
 
!-- sauvegarde des flux: 
#ifdef CALC
    if (allocated(coord_flux)) then
       include './BIO/call_save_flux.inc'
    endif
#endif
!-------------------------------------------------------------------  
 case(5) ! Affichage des variables physiques (irradiance, temp, salinit�)
         ! a chaque pas de temps:
      write(24,'(F8.2,1X,F8.4,1X,F7.4,1X,F7.4,1X,F6.2,1X,F6.2)')tps/3600.d0,E_PARZ(1,1,1),TEMP_BIO(1,1,1),SAL_BIO(1,1,1),WIND(1,1,1),CO2atm(1,1,1)

!-------------------------------------------------------------------
  ! A FAIRE SI NECESSAIRE !

end select

End Subroutine Sub_affiche
