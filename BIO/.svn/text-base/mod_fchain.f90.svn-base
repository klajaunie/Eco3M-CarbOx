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
!

                 MODULE MOD_FCHAIN
!
! Module contenant les differentes fonctions
! de traitement de chaines de caracteres et celles de
! recherche d'elements particuliers dans un tableau
!
! Derniere modification 27/06/07
!--------------------------------------------------------------------
USE LONG_CHAIN

Implicit None

 CONTAINS

!--------------------------------------------------------------------
!
              Character(L_CHAIN) function f_chain(chain,ipos,separateur)
!
!--------------------------------------------------------------------
!
! Cette fonction permet de retourner la sous chaine de 
! caracteres situee entre le ipos-ieme-1 : et le (ipos-ieme) separateur
! de la chaine mere chain
!
 Implicit None
 Character(L_CHAIN) :: chain
 Character(1)   :: separateur
 Integer        :: ipos

 character(L_CHAIN)  :: ch_tempo
 
 !--Variables locales:
   Integer::ii,ideb,ifin,ip

 ch_tempo = adjustl(chain)
 ch_tempo = trim(ch_tempo) // separateur
  
 ideb = 1
 ifin = len_trim(ch_tempo)
 ip = 0
 do
  ii=index(ch_tempo(ideb:ifin),separateur)+ideb-1
  ip = ip+1
  if (ip == ipos) then
    f_chain=ch_tempo(ideb:ii-1)
    f_chain=adjustl(f_chain)
    f_chain=trim(f_chain)
    exit
  endif
  ideb=ii+1
 enddo
 Return
 End 	function f_chain
!--------------------------------------------------------------------
!
             Integer function f_nschain(chain,separateur)
!
!--------------------------------------------------------------------
! Cette fonction permet de retourner le nb de sous chaine de 
! caracteres (separees par des :) contenues dans la chaine mere chain
!
 Implicit None
 Character(L_CHAIN)  :: chain
 Character(1)        :: separateur
 character(L_CHAIN)  :: ch_tempo2

 !--Variables locales:
   Integer::ii,ideb,ifin,ip

if(len_trim(chain)/=0) then
 ch_tempo2 = adjustl(chain)
 ch_tempo2 = trim(ch_tempo2) // separateur

 ideb = 1
 ifin = len_trim(ch_tempo2)
 ip = 0 ! nb d'occurence du caractere "separateur"
 do
  ii=index(ch_tempo2(ideb:ifin),separateur)+ideb-1
  ip = ip+1
  ideb=ii+1
  if (ideb >= ifin) exit
 enddo
 
else
  ip=0
endif
 f_nschain = ip
 
 Return
 End 	function f_nschain
!--------------------------------------------------------------------
!
             Integer function f_nschain_seq(chain,sequence,lseq)
!
!--------------------------------------------------------------------
! Cette fonction est identique a f_nschain mais recherche une
! sequence de plusieurs (i.e. lseq ) caracteres dans la chaine mere
!
 Implicit None
 integer             :: lseq
 Character(L_CHAIN)  :: chain
 Character(lseq)     :: sequence
 character(L_CHAIN)  :: ch_tempo2
 
 !--Variables locales:
   Integer::ii,ideb,ifin,ip

if(len_trim(chain)/=0) then
 ch_tempo2 = adjustl(chain)
 ch_tempo2 = trim(ch_tempo2) // sequence

 ideb = 1
 ifin = len_trim(ch_tempo2)
 ip = 0 ! nb d'occurence du caractere ":"
 do
  ii=index(ch_tempo2(ideb:ifin),sequence)+ideb-1
  ip = ip+1
  ideb=ii+lseq
  if (ideb >= ifin) exit
 enddo
 
else
  ip=0
endif
 f_nschain_seq = ip
 
 Return
 End 	function f_nschain_seq
!--------------------------------------------------------------------
!
             Integer function f_proc2id(nom_proc)
!
!--------------------------------------------------------------------
! Cette fonction renvoie le numéro du processus associé au nom
! Si le nom="0", le numéro renvoyé est 0

use DEF_TYPE
Implicit none
 Character(L_VAR)     :: nom_proc
 character(L_VAR)     :: tempo

integer          :: i

f_proc2id=0
  do i=1,size(PROC_MOD)  
    tempo = PROC_MOD(i)%nomproc
    if(tempo==trim(adjustl(nom_proc))) then
      f_proc2id = PROC_MOD(i)%idproc
      exit

    end if
  end do
  if(f_proc2id==0) then
   write(*,*)'vous avez utilisé un processus qui n''existe pas &
&dans le modèle :',nom_proc
    write(*,*) 'Ce programme va être arreté'
    read(*,*)
    stop
  endif

End function f_proc2id

!-------------------------------------------------------------------------------------------
!
     character (L_CHAIN)  Function f_Int2chain (int,longchain)
!
!-------------------------------------------------------------------------------------------
! Fonction qui convertit un entier "int" en chaine de caracteres 
! "chain" calee a gauche 
! et dont la longueur (sans les blancs) est longchain

   Implicit None
!	character(200) :: f_Int2chain
character(L_CHAIN) :: chain
Integer            :: int,longchain
  
write(chain,*) int
chain = trim(adjustl(chain))
longchain = len_trim(chain)
f_Int2chain = chain(1:longchain)
  
End Function f_Int2chain

!-------------------------------------------------------------------------------
!
      Integer  Function f_scomp2id (scomp_in,elmt_in)
!
!-------------------------------------------------------------------------------
! Fonction qui permet de trouver l'identifiant id du tableau VAR, 
! correspondant à l'élément elmt (rentré par l'utilisateur),
! du sous-compartiment scomp (rentré par l'utilisateur)
!
! Exemple : on recherche l'Azote N du sous-compartiment diatomés...
   
 USE COUPLEUR_PHY_BIO

 Implicit None
 character(L_CHAIN)   :: elmt_in
 character(L_CHAIN)   :: scomp_in 
 character(L_CHAIN)   :: tempo,tempo2,tempo3,tempo4
 integer              :: i,j,k

 f_scomp2id = 0
 i=1
 do while (i<=size(VAR))
    if( trim(adjustl(VAR(i)%scomp)) == trim(adjustl(scomp_in)).AND. &
         trim(adjustl(VAR(i)%elmt)) == trim(adjustl(elmt_in)) ) then
           f_scomp2id = i
           exit
    endif
    i=i+1
 Enddo

 Return 
 End Function f_scomp2id
!-------------------------------------------------------------------------------	
!
      Integer  Function f_idorg2id (idorg_in,elmt_in)
!
!-------------------------------------------------------------------------------
! Fonction qui permet de trouver l'identifiant id du tableau VAR 
! correspondant à la variable dont la biomasse est exprimee
! en elmt_in (input),du numero de sous-compartiment idorg (input)
!
! Exemple : on recherche l''identifiant id de la variable correspondant
!           au sous-compartiment diatomés (idorg) exprime 
!           en termes d'azote (elmt_in)
   
 USE COUPLEUR_PHY_BIO
 Implicit None
 
!-- arguments:
 character(L_VAR)     :: elmt_in
 integer              :: idorg_in 
!-- variables locales
 integer              :: i, long, long_in

 f_idorg2id = 0
 
 i=1

 do while (i<=size(VAR))  
    long = len_trim(trim(adjustl(VAR(i)%elmt)))
    long_in = len_trim(trim(elmt_in))
 
    if ( VAR(i)%idorg == idorg_in .AND. &
       (VAR(i)%elmt(1:long) == elmt_in(1:long_in) )) then
         f_idorg2id = i
         exit
    endif
    i=i+1
 enddo
 
 return 
 End Function f_idorg2id
!-------------------------------------------------------------------------------	
!
       Function f_idorg2id_vect (idorg_in,elmt_in,size_idorg)
!
!-------------------------------------------------------------------------------
! Fonction identique a f_idorg2id mais qui s'applique pour une entree
! vectorielle
   
 USE COUPLEUR_PHY_BIO
 Implicit None
 
!-- arguments:
 character(L_VAR)     :: elmt_in
 integer              :: size_idorg,idorg_in(size_idorg) 
 Integer              :: f_idorg2id_vect(size_idorg)
!-- variables locales
 integer              :: i, long, long_in,ll

 f_idorg2id_vect = 0

 do ll=1,size_idorg
 i=1
   do while (i<=size(VAR))  
      long = len_trim(trim(adjustl(VAR(i)%elmt)))
      long_in = len_trim(trim(elmt_in))
 
      if ( VAR(i)%idorg == idorg_in(ll) .AND. &
        (VAR(i)%elmt(1:long) == elmt_in(1:long_in) )) then
          f_idorg2id_vect(ll) = i
          exit
    endif
    i=i+1
   enddo
 enddo
 return 
 End Function f_idorg2id_vect
!-------------------------------------------------------------------------------
 END MODULE MOD_FCHAIN
