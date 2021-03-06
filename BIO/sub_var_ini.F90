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
!
                    SUBROUTINE  Sub_var_ini
!
!  * Cette routine permet de lire dans les fichiers ad hoc dont les
! noms doivent etre de la forme COMP_SCOMP_ELMT.DAT (tout en majuscules),
! les concentrations initiales pour l''ensemble des variables sous reserve
! que le fichier correspondant a chaque variable existe. En cas 
! contraire, la conc. initiale est d''office mise a une valeur
! cini standart
! 
! Derni�re modification: 27/06/07
!
!-------------------------------------------------------------						  
!
 Use DEF_TYPE
 Use VAR_USER
 Use COUPLEUR_PHY_BIO

 Implicit None
 
 ! Variables locales
 Integer        :: iphy,id,nfich1,nfich2,nfich3,i,j,k
 Integer        :: imat,jmat
 integer        :: errlec,istat
 Character(10)  :: fich1,fich2,fich3
 Character(70)  :: fich
 Character(20)  :: tempo
 Logical        :: present
 
!--------------------------------------------------------------------
! ALLOCATION DU TABLEAU DE NOMS de FICHIERS des CONCENTRATIONS POUR 
! CHAQUE VARIABLE d''ETAT
!--------------------------------------------------------------------
write(*,*)'passage dans sub_var_ini'
#ifdef COUPL   
  Allocate(tab_noms_fich(nbvar), STAT=istat)
  if (istat /= 0) write(*,*) 'probleme d allocation de tab_noms_fich'
  write(11,*) 'alloc de tab_noms'
#endif
!--------------------------------------------------------------------
!   FICHIERS DE CONCENTRATIONS INITIALES
!--------------------------------------------------------------------

do id=1,nbvar
   fich1 = VAR(id)%comp
   nfich1 = len_trim(adjustl(fich1))
!
   fich2 = VAR(id)%scomp
   nfich2 = len_trim(adjustl(fich2))
!
   fich3 = VAR(id)%elmt
   nfich3 = len_trim(adjustl(fich3))
!
!-- Creation du nom du fichier de conc. initiale pour la variable
!-- numero id sous la forme comp_scomp_element.DAT
   fich = trim(repconcinit) // fich1(1:nfich1)//'_'//fich2(1:nfich2)//'_'//fich3(1:nfich3)//'.DAT'

#ifdef COUPL   
     tab_noms_fich(id) = fich
#endif
!
#ifdef NCOUPL
!-- test d''existence du fichier
   INQUIRE(file=fich,EXIST=present)
!
   if (.NOT. present) then
        VAR(id)%conc = cini
        WRITE(22,*) fich ,'non trouve. => conc. prise par defaut a',cini
        WRITE(*,*) fich ,'non trouve. => conc.(',id,') prise par defaut a',cini
        read(*,*)
   else
        CALL sub_lec_cini(id,fich)
   endif
   
#endif   

 enddo
!--------------------------------------------------------------------
!-- Initialisation de la temperature et la salinite en 0D: 
!--------------------------------------------------------------------
#ifdef NCOUPL
     SAL_BIO=36.d0
     TEMP_BIO=20.0d0
#endif

#ifdef MODTEST
  write(*,*)  " -- Sub_var_ini Termin�e ! -- "
#endif

END Subroutine Sub_var_ini
!-----------------------------------------------------------------------
!
                 SUBROUTINE Sub_lec_cini(id,fichier)
!
! Cette sous-routine charge les fichiers de concentration initiale
! Donn�es sous forme ASCII
! FORMAT du fichier : i,j,k,Cini
! 
!-----------------------------------------------------------------------

 Use DEF_TYPE
 Use COUPLEUR_PHY_BIO
 Use VAR_USER
 Use VAR_GLOBAL
 
 !-- arguments
 Character(*)    ::fichier
 Integer         :: id
 
 !-- Variables locales:
 integer         :: i,j,k
 real(8)         :: conc
 integer         :: errlec
 character(200)  :: tempo

!-- boucle de lecture
   open(135,FILE=fichier)
 do
   Read(135,'(a)',iostat=errlec) tempo
   if (tempo==' ') exit
   if (errlec /=0) then 
      write(*,*)'pb de lecture du fichier de concentration initiale :', fichier
      stop
   endif
   Read(tempo,*) i,j,k,conc
!print *, i,j,k,conc
   VAR(id)%conc(i,j,k) = conc
   write(*,*)'conc ini',id,i,j,k,conc
 enddo


 close(135)
END Subroutine  Sub_lec_cini
