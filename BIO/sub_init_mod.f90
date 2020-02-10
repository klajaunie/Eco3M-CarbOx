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
                       SUBROUTINE Sub_init_mod
!				
!-- Routine définissant les differents modeles de processus disponibles
!   après lecture du fichier fichconfmod (de nom "modele.def" par defaut)
!
!-- ATTENTION: tous les processus apparaissant dans le fichier  fichconfmod 
!            (de nom "modele.def" par defaut) doivent avoir été déclarés au
!             préalable dans le fichier mod_process.f90  
!
! - Les différents modèles de processus sont ensuite référencés dans la matrice PROC_MOD
!   qui permet l association entre les noms des processus et les identifiants des
!   processus (nombres entiers)
!
!-- Dernière modif: 06/07/07
!--------------------------------------------------------------------------------
				
Use DEF_TYPE
Use VAR_USER
Use MOD_FCHAIN
Use VAR_GLOBAL

 Implicit None
 
 !--Variables locales
 integer            :: errlec,istat
 integer            :: nb_proc
 integer            :: i,j,itemp
 character(L_CHAIN) :: chaine,chaine_tempo

!-------------------------------------------------------
! Initialisation de la matrice PROC_MOD
! Ce tableau contient les associations entre 
! nom_proc et id_proc
!
!-------------------------------------------------------
 !-- Ouverture du fichier de données:

 Open(116,FILE =fichconfmod)
 write(*,*) 'fichier contenant les modeles de processus =',fichconfmod
 
!-- initialisation du nb de fonctions de processus:
 nb_proc =0

!-- Calcul du nombre de processus
 do
   read(116,*,iostat=errlec) chaine
   if (errlec /=0) stop 'pb de lecture du fichier de config sub_init_mod nb'
   if (chaine(1:1)=='#') then
        cycle
    elseif (chaine(1:4)=='!fin' .OR. chaine(1:4)=='!FIN'.OR. chaine(1:4)=='!Fin') then
        exit
    else
        nb_proc = nb_proc + 1
    endif
 end do

 !--Creation de la matrice de processus, liant un nom avec un numéro
 !--  et le nombre et les noms des parametres associes
 if (nb_proc /=0)  Allocate(PROC_MOD(nb_proc))

!-- sauvegarde du nombre de modeles disponibles dans la bibliotheque:
nbproc_bib= nb_proc

 rewind(116)

 i=0
 do while (i < nb_proc)   
   read(116,*,iostat=errlec) chaine
    if (errlec /=0) stop 'pb de lecture du fichier de config dans sub_init_mod '
    if (chaine(1:1)=='#') cycle
    i = i+1
    PROC_MOD(i)%idproc = i
      chaine_tempo=f_chain(chaine,1,':')
      read(chaine_tempo,*)PROC_MOD(i)%nomproc
      chaine_tempo=f_chain(chaine,2,':')
      read(chaine_tempo,*)PROC_MOD(i)%nomsub
      chaine_tempo=f_chain(chaine,3,':')
      read(chaine_tempo,*)PROC_MOD(i)%nbpar
      itemp = PROC_MOD(i)%nbpar
!-- Allocation de l'element i du tableau de pointeurs PROC_MOD%nompar:
      if (associated (PROC_MOD(i)%nompar)) NULLIFY(PROC_MOD(i)%nompar)
      Allocate(PROC_MOD(i)%nompar(itemp),STAT=istat)
      if (istat /= 0) write(*,*) 'pd d allocation de PROC_MOD(',i,')%nompar'
      write(11,*) 'alloc de PROC_MOD(',i,')%nompar'

!
      do j=1,PROC_MOD(i)%nbpar
         chaine_tempo=f_chain(chaine,3+j,':')
         read(chaine_tempo,*)PROC_MOD(i)%nompar(j)
      enddo
 end do 

!-- Fermeture du fichier
   close(116)


   write(*,*)  " -- Sub_init_mod Terminé ! -- "
   write(*,*)  'nb_proc =',nb_proc
   write(*,*)  
   
End Subroutine Sub_init_mod
