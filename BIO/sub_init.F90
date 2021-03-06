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

                          SUBROUTINE Sub_init

!-- Routine d''initialisation des choix de l''utilisateur
!   mentionn�s dans le fichier ''config.ini''
!
!   Initialisation des matrices  :
!
! - VAR(nbvar)            : Variables d �tat du mod�le
! - FLUX_VAL(nbvar,nbvar) : interactions entre variables d �tat
! - SELF_VAL(nbvar)       : termes sources ou puits ne faisant intervenir
!                           qu''une seule var. d''etat
!
! derni�re modification: 27/06/07
!-------------------------------------------------------------------------------
!
! Variables globales
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_USER
Use VAR_GLOBAL
Use MOD_FCHAIN
 
Implicit None
 
! Variables locales
 Integer::i,j,k,istat
 Integer::ili,jcol,iii
 Integer::icomp, iscomp,iscomp1,ivar,inbelmt,inbelmt2,nborg,inbmult
 Integer::idprocmat,idtemp,nbproctot
 Integer::imat,jmat
 Integer::iposini,iposfin
 Integer::errlec
 Integer            :: CHL_C_VAL
 integer            :: i_st,j_st,k_st   
 Integer            :: ii,nbpoints
 integer            :: lon_cist, lon_cjst, lon_ckst
 Character(L_VAR)   :: nomcomp
 Character(L_VAR)   :: signe
 Character(L_CHAIN) :: chaine,fichtemp
 Character(L_CHAIN) :: tempo,tempo2,tempopo,nomproc_tempo,chaine_proc
 Character(L_VAR)   :: processus
 character(L_CHAIN) :: c_ist,c_jst,c_kst


!-- Ouverture du fichier de donn�es:
 Open(15,FILE =fichconf)
      write(*,*)'fichier de configuration =',fichconf
 
!-------------------------------------------------------
!  Lecture des donn�es dans le fichier de configuration
!-------------------------------------------------------
!
!-- Lecture de parametres divers
!-------------------------------------------------------
 do 
  Read(15,*,iostat=errlec)chaine
  if (errlec /=0) stop 'pb de lecture du fichier de config, param divers'
  if (chaine(1:1)=='#') cycle
   suffbio=chaine !-- lecture du suffixe associe a la simulation
   exit
 enddo
!-------------------------------------------------------
!-- Lecture des param�tres temporels
!-------------------------------------------------------
do 
     Read(15,*,iostat=errlec)chaine
#ifdef MODTEST
     write(*,*) chaine
     
#endif     
     if (errlec /=0) stop 'pb de lecture du fichier de config,param temporels'
     if (chaine(1:1)=='#') cycle
     tempo=f_chain(chaine,1,':')
     read(tempo,*) duree !-- duree simul (jours)
     tempo=f_chain(chaine,2,':')
     read(tempo,*) dt_bio !-- pas de tps calcul biogeo (sec)
     tempo=f_chain(chaine,3,':')
     read(tempo,*) dt_save_bio !-- pas de tps de sauvegarde bio(mn)
     exit
enddo 
!-------------------------------------------------------
!--   Lecture des dimensions spatiales
!-------------------------------------------------------
do 
     Read(15,*,iostat=errlec)chaine
#ifdef MODTEST
     write(*,*) chaine
#endif     
     if (errlec /=0) stop 'pb de lecture du fichier de config,param spatiaux'
     if (chaine(1:1)=='#') cycle
   #ifdef NCOUPL
     tempo=f_chain(chaine,1,':')
     read(tempo,*) nx_min
     tempo=f_chain(chaine,2,':')
     read(tempo,*)nx_max
     tempo=f_chain(chaine,3,':')
     read(tempo,*)ny_min
     tempo=f_chain(chaine,4,':')
     read(tempo,*)ny_max
     tempo=f_chain(chaine,5,':')
     read(tempo,*)nz_max
   #endif
     exit
 enddo
!-------------------------------------------------------
!-- Lecture du nom du fichier d''irradiance
!-------------------------------------------------------
 do
     Read(15,*,iostat=errlec)chaine
#ifdef MODTEST
     write(*,*) chaine
#endif
     if (errlec /=0) stop 'pb de lecture du nom du fichier d irradiance'
     if (chaine(1:1)=='#') cycle
     tempo=f_chain(chaine,1,':')
     read(tempo,*) fichirrad
     tempo=f_chain(chaine,2,':')
     if (fichirrad == 'IRR_FONCTION') then
     write(*,*)'irrad fonction'
       read(tempo,*) irrad_MAX
!-- pas de temps de calcul ou de lecture d irrad
       tempo=f_chain(chaine,3,':')
       read(tempo,*)dt_irrad
       tempo=f_chain(chaine,4,':')
       read(tempo,*)irr2par
       tempo=f_chain(chaine,5,':')
       read(tempo,*)albedo
       tempo=f_chain(chaine,6,':')
       read(tempo,*)irr_param
     elseif (fichirrad == 'IRR_CODEPHYS') then
     write(*,*)'irrad code physique'
       read(tempo,*)irr2par
       tempo=f_chain(chaine,3,':')
       read(tempo,*)albedo
       tempo=f_chain(chaine,4,':')
       read(tempo,*)irr_param
     else
     write(*,*)'irrad fichier'
!-- pas de temps de calcul ou de lecture d irrad
       read(tempo,*)dt_irrad
       tempo=f_chain(chaine,3,':')
       read(tempo,*)irr2par
       tempo=f_chain(chaine,4,':')
       read(tempo,*)albedo
       tempo=f_chain(chaine,5,':')
       read(tempo,*)irr_param
     endif
 !-- Ajout du chemin
     fichirrad_long = './BIO/DATA/'//fichirrad
 !-- conversion en secondes
     dt_irrad = dt_irrad * 60
     exit
 enddo

!------------------------------------------------------------------
!-- Lecture du nom de la fonction calculant l''extinction lumineuse
!   valable pour les differents modes sauf le mode 0D
!------------------------------------------------------------------
#ifndef M0D-NCOUPL
 do 
     write(*,*) 'pas M0D-NCOUPL'
     Read(15,*,iostat=errlec)chaine
     if (errlec /=0) stop 'pb de lecture du nom du fichier d''irradiance'
     if (chaine(1:1)=='#') cycle

!-- numero du processus
     tempo= f_chain(chaine,1,'(')
     tempo = trim(adjustl(tempo))
     write(*,*) 'tempo =',tempo
     tempo2=chaine(len_trim(tempo)+2:len_trim(chaine)-1)
     write(*,*) 'tempo2 =',tempo2

!-- affectation a id_extinc du numero de processus:
        idtemp = f_proc2id(tempo)
	IRR_PAR%idproc = idtemp 
        inbelmt= PROC_MOD(idtemp)%nbpar
        Allocate(IRR_PAR%valpar(inbelmt),STAT=istat)
        if (istat /= 0) write(*,*) 'pb d''allocation de IRR_PAR%valpar'
        write(11,*) 'alloc de IRR_PAR%valpar a ',inbelmt
         do k= 1, inbelmt 
            tempopo=f_chain(tempo2,k,'>')
            Read(tempopo,*)IRR_PAR%valpar(k)
         enddo
       exit
 enddo
#endif
!-------------------------------------------------------
! Specificite du maillage sur la verticale:
! indique le sens (fond--> surf ou surf-->fond) adopte dans
! le code physique pour la numerotation des mailles sur la verticale
!-------------------------------------------------------
!
#ifndef M0D-NCOUPL
errlec = 0
  do 
write(*,*) 'pas M0D-NCOUPL'

     Read(15,*,iostat=errlec) chaine
     if (errlec /=0) stop 'pb de lecture de la position de nz_max'
     if (chaine(1:1)=='#') cycle
     read(chaine,*) pos_nzmax
     write(*,*) 'pos =',pos_nzmax
     exit
 enddo 
#endif
!
!
!-------------------------------------------------------
! Determination de l'�chelle de pH
! pHscale=1 ==> SWS pH scale; else TOT pH scale
!-------------------------------------------------------
!
!#ifndef M0D-NCOUPL
!errlec = 0
!  do 
!write(*,*) 'pas M0D-NCOUPL'
!
!     Read(15,*,iostat=errlec) chaine
!     write(*,*)'test',chaine
!     if (errlec /=0) stop 'pb de lecture de la position de pH scale'
!     if (chaine(1:1)=='#') cycle
!     read(chaine,*) pHscale
!     write(*,*) 'pHscale =',pHscale
!     exit
! enddo 
!#endif
!
!
!-------------------------------------------------------
!-- Lecture du nb de compartiments du modele
!-------------------------------------------------------
 errlec = 0
  do 
     Read(15,*,iostat=errlec)chaine
     write(*,*)chaine
     if (errlec /=0) stop 'pb de lecture du fichier de config, nbre compartiments'
     if (chaine(1:1)=='#') cycle
     tempo= f_chain(chaine,1,':')
     read(tempo,*) nbcomp
     tempo=f_chain(chaine,2,':')
     read(tempo,*) nbscomp
     tempo=f_chain(chaine,3,':')
     read(tempo,*) nbvar
     exit
  enddo
!
!-------------------------------------------------------
! Initialisation de la matrice VAR
!-------------------------------------------------------
!
!-- Allocation du tableau contenant les variables d �tat
ALLOCATE(VAR(nbvar),STAT=istat)
if (istat /=0) write(*,*) 'pb d''allocation de VAR dans sub_init'
   write(11,*) 'alloc de VAR  a ',nbvar
#ifdef COUPL
   ALLOCATE(nb_elmt(nbscomp),STAT=istat)
   if (istat /=0) write(*,*) 'pb d''allocation de nb_elmt dans sub_init'
   write(11,*) 'alloc de nb_elmt  a ',nbscomp
#endif

!--Identification des compartiments et ss compartiments:
!-- ET remplissage de la matrice VAR(nbvar)
icomp  = 0
iscomp = 0
iscomp1= 0
ivar   = 0
nborg  = 0

 do while (icomp < nbcomp)
     Read(15,*,iostat=errlec)chaine
     if (errlec /=0) stop 'pb de lecture du fichier de config, matrice VAR'
     if (chaine(1:1)=='#') cycle
     icomp = icomp+1
     nomcomp = f_chain(chaine,1,':')
     nomcomp = adjustl(nomcomp)
     tempo=f_chain(chaine,2,':')
     Read(tempo,*)iscomp

   !-- pour chaque sous-compartiment:
     do while (iscomp1 < iscomp)
     Read(15,*,iostat=errlec)chaine
       if (errlec /=0) stop 'pb de lecture du fichier de config,matrice VAR,2'
       if (chaine(1:1)=='#') cycle
       iscomp1=iscomp1+1
       inbelmt=f_nschain(chaine,':')
       nborg = nborg + 1 ! compte le nb d''organismes differents 

#ifdef COUPL
       nb_elmt(nborg)=inbelmt-1 !-- compte le nb d elements associes a
                               !   chaque organisme
#endif

!-- pour chaque element:
       do k=1,inbelmt-1
         ivar=ivar+1
         VAR(ivar)%comp=nomcomp
         VAR(ivar)%scomp=f_chain(chaine,1,':')
         VAR(ivar)%elmt =f_chain(chaine,k+1,':')
         VAR(ivar)%idorg = nborg 
       enddo
    enddo
iscomp1=0
enddo

!-- Verifications:	
 if (ivar /= nbvar) stop 'probleme dans le nb de variables totales'

!-------------------------------------------------------
! Initialisation de la matrice FLUX_VAL
! Ce tableau contient les interactions entre les variables
! d �tat du mod�le, ie les processus
!
! Cette matrice est diagonale superieure (sans elements sur la diag) 
!
!-------------------------------------------------------

!-- Allocation du tableau contenant les interactions entre variables
if (.not. allocated (FLUX_VAL)) ALLOCATE(FLUX_VAL(nbvar,nbvar),STAT=istat)
write(11,*) 'alloc de FLUX_VAL  a ',nbvar,nbvar
if (istat /= 0) write(*,*) 'pb d''allocation de FLUX_VAL'
!--Lecture de la matrice de l utilisateur
!-- ET remplissage de la matrice FLUX_VAL(nbvar,nbvar)
inbelmt = 0
inbelmt2 = 0
idprocmat=0
    
!-- Desassociation des pointeurs:

do ili=1,nbvar
  do jcol=1,nbvar
     NULLIFY(FLUX_VAL(ili,jcol)%idproc)
!     write(11,*) 'nullify de FLUX_VAL(',ili,jcol,')'
  enddo
enddo

 do  
     Read(15,*,iostat=errlec) chaine
      if (errlec /=0) stop 'pb de lecture du fichier de config, matrice FLUX_VAL'
      if (chaine(1:1)=='#') cycle
      if (trim(adjustl(chaine))=='finflux') exit
!-- lecture du nb cumule de processus utilises dans l ensemble des  
!-- cellules de la matrice de flux :
      Read(chaine,*)nbproctot 
      nbproc_flux=nbproctot
      if (.not. allocated (FLUX_PAR)) Allocate(FLUX_PAR(nbproctot),STAT=istat)
      write(11,*) 'alloc de FLUX_PAR  a ',nbproctot
      if (istat /= 0) write(*,*) 'pb d''allocation de FLUX_PAR'
      exit
 enddo

 do 
!-- Lit en lignes chaque processus de la matrice flux ds le fichier config.ini
      Read(15,*,iostat=errlec) chaine
#ifdef MODTEST
  write(*,*) 'chaine = ',chaine
#endif  
      if (errlec /=0) stop 'pb de lecture du fichier de config, matrice FLUX_VAL,2'
      if (chaine(1:1)=='#') cycle
      if (trim(adjustl(chaine))=='finflux') exit 

  !-- lecture: ili,jcol
      tempo=f_chain(chaine,1,':')
      tempo2=f_chain(tempo,1,';')
      read(tempo2,*)ili !numero de ligne
      tempo2=f_chain(tempo,2,';')
      read(tempo2,*)jcol !numero de col	   	  
!------------------------------------------------------------------	  	  	    
!   calcul du nombre de fonctions et de processus dans la cellule
!------------------------------------------------------------------
!-- on isole la portion de chaine ne contenant que les processus et leur signe: 
      chaine_proc=chaine(len_trim(tempo)+2:len_trim(chaine))

!-- calcul du nb de sous chaines (et donc de processus differents)
      inbelmt=f_nschain(chaine_proc,':')
      
 ! nb de signes multiplie ds la chaine  (on cherche une chaine de 3 caracteres au cas 
 ! ou il y aurait la lettre 'x' dans un nom de fonction   
       inbmult=f_nschain_seq(chaine_proc,':x;',3) - 1
      
!-- calcul du nb de "sous-flux" constituant le flux de ili vers jcol       
      FLUX_VAL(ili,jcol)%nb_sflux = inbelmt - inbmult
      
!-- Allocation dynamique en fonction du nb de processus
      Allocate (FLUX_VAL(ili,jcol)%idproc(inbelmt),STAT=istat)
      if (istat /= 0) write(*,*) 'pb d''allocation de FLUX_VAL%idproc'
      write(11,*) 'alloc de FLUX_VAL(',ili,jcol,')%idproc  a ',inbelmt

!-- Boucle sur le nb de processus de la cellule (ili,jcol) de la matrice flux
      do i=1,inbelmt
!-- Incrementation du nb de processus cumules:
       idprocmat = idprocmat + 1
  
!-- test sur le nb cumule de processus
      if (idprocmat > nbproctot) stop 'nb de processus cumules de la matrice > a celui indique dans matrice FLUX'
 
!-- Stockage des coordonnees (ili,jcol,signe) relatifs au processus n�idprocmat
       !-- lecture du signe
       tempo=f_chain(f_chain(chaine,i+1,':'),1,';')
       read(tempo,*)signe !signe du flux
       FLUX_PAR(idprocmat)%ipos(1)=ili
       FLUX_PAR(idprocmat)%ipos(2)=jcol
!
       FLUX_PAR(idprocmat)%signe= signe  
  
!--Associe le nom du proc � son num�ro
       nomproc_tempo = f_chain(f_chain(chaine,i+1,':'),2,';')
       tempo = f_chain(nomproc_tempo,1,'(') !coupure a la parenthese
  !-- numero du processus
       idtemp = f_proc2id(tempo)
       FLUX_VAL(ili,jcol)%idproc(i)=idtemp
       FLUX_PAR(idprocmat)%idproc=idtemp

  !-- recherche du nb de parametres
       iposini=index(nomproc_tempo,'(') !--position d ouverture de parenthese
       iposfin=index(nomproc_tempo,')') !--position de fermeture de parenthese
       tempo=nomproc_tempo(iposini+1:iposfin-1)
       inbelmt2=f_nschain(tempo,'>')!--nb de parametres pour le proc. en cours
!-- on verifie que ce nb de param. est correct
       if (inbelmt2/=PROC_MOD(idtemp)%nbpar) then
         write(*,*) 'nb de parametres du processus ',PROC_MOD(idtemp)%nomproc,&
&'dans la cellule ',ili,jcol,' est incoherent'
         write(22,*)'nb de parametres du processus ',PROC_MOD(idtemp)%nomproc,&
&'dans la cellule ',ili,jcol,' est incoherent'
#ifdef MODTEST
       read(*,*)
#endif      
 endif
!         STOP 
!       else
         if (inbelmt2 >=1) then
            if (associated (FLUX_PAR(idprocmat)%valpar)) NULLIFY(FLUX_PAR(idprocmat)%valpar)
            Allocate(FLUX_PAR(idprocmat)%valpar(inbelmt2),STAT=istat)
            if (istat /= 0) write(*,*) 'pb d allocation de FLUX_PAR%valpar'
            do k= 1, inbelmt2 
              tempopo=f_chain(tempo,k,'>')
              Read(tempopo,*)FLUX_PAR(idprocmat)%valpar(k)
            enddo
         else
            if (associated (FLUX_PAR(idprocmat)%valpar)) NULLIFY(FLUX_PAR(idprocmat)%valpar)
         endif
!       endif
    enddo    
 enddo
!
!-------------------------------------------------------
! Initialisation du vecteur SELF_VAL
! Ce vecteur contient les termes source ou puit supplementaires 
! dont l origine ou la destination ne fait pas partie des variables
! d''etat du modele.
!-------------------------------------------------------


!-- Allocation du vecteur SELF_VAL
      Allocate(SELF_VAL(nbvar),STAT=istat)
      write(11,*) 'alloc de SELF_VAL  a ',nbvar
      if (istat /= 0) write(*,*) 'pb d''allocation de SELF_VAL'

!-- Desassociation des pointeurs:
      do iii=1,nbvar
       NULLIFY(SELF_VAL(iii)%idproc)
      write(11,*) 'nullify de SELF_VAL(',iii,')%idproc  '
      enddo

!--Lecture du vecteur SELF(nbvar) d�fini par l utilisateur
!-- ET remplissagede ce vecteur

inbelmt = 0
inbelmt2 = 0
idprocmat=0
errlec = 0

do 
!-- 1ere lecture de 'config.ini' pour allocation de la matrice SELF :
     Read(15,*,iostat=errlec) chaine
      if (errlec /=0) stop 'pb de lecture du fichier de config dans SELF_VAL'
      if (chaine(1:1)=='#') cycle
      if (trim(adjustl(chaine))=='finself') exit 
!---  nb cumule de processus utilises sur l ensemble des cellules de la matrice:
      Read(chaine,*) nbproctot  
      nbproc_self=nbproctot
      Allocate(SELF_PAR(nbproctot),STAT=istat) 
      if (istat /= 0) write(*,*) 'pb d''allocation de SELF_PAR'
      write(11,*) 'alloc de SELF_PAR ',nbproctot
      exit
enddo
 do
!-- 2nde lecture de 'config.ini' pour remplissage de la matrice SELF :
      Read(15,*,iostat=errlec) chaine
#ifdef MODTEST
     write(*,*) chaine
#endif 
      if (errlec /=0) stop 'pb de lecture du fichier de config dans SELF_VAL'
      if (chaine(1:1)=='#') cycle
      if (trim(adjustl(chaine))=='finself') exit 
  !-- lecture: ili
      tempo=f_chain(chaine,1,':')
      tempo2=f_chain(tempo,1,';')
      read(tempo2,*)ili !numero de ligne
      
!------------------------------------------------------------------
!-- calcul du nombre de fonctions appelees et de processus dans la cellule:
!------------------------------------------------------------------
!  on isole la chaine contenant que les processus et leur signe (sans i;j)
      chaine_proc=chaine(len_trim(tempo)+2:len_trim(chaine))
!-- calcul du nb de fonctions utilisees      
      inbelmt=f_nschain(chaine_proc,':')
! nb de signes multiplie ds la chaine  (on cherche une chaine de 3 caracteres au cas 
 ! ou il y aurait la lettre 'x' dans un nom de fonction   
       inbmult=f_nschain_seq(chaine_proc,':x;',3) - 1

!-- calcul du nb de "sous-flux" constituant le flux de ili vers une variable       
      SELF_VAL(ili)%nb_sflux = inbelmt - inbmult

      Allocate (SELF_VAL(ili)%idproc(inbelmt),STAT=istat)
      if (istat /= 0) write(*,*) 'pb d allocation de SELF_VAL%idproc'
      write(11,*) 'alloc de SELF_VAL(',ili,')%idproc  a ',inbelmt
!
    do i=1,inbelmt
!-- Incrementation du nb de processus cumules:
       idprocmat = idprocmat + 1
!--   test sur le nb cumule de processus
       if (idprocmat > nbproctot) STOP 'nb de processus cumules de la matrice > a celui indique dans SELF'

!--   Stockage de la coordonnee (ili,signe) correspondant au processus n idprocmat
!     lecture du signe
       tempo=f_chain(f_chain(chaine,i+1,':'),1,';')
       read(tempo,*)signe !signe du flux

       SELF_PAR(idprocmat)%ipos(1)=ili
       SELF_PAR(idprocmat)%ipos(2)=0
       SELF_PAR(idprocmat)%signe= signe 

   !--Associe le nom du proc � son num�ro
       nomproc_tempo = f_chain(f_chain(chaine,i+1,':'),2,';')
       tempo = f_chain(nomproc_tempo,1,'(') !coupure a l ouverture de parenthese
   !-- numero du processus
       idtemp = f_proc2id(tempo)
       SELF_VAL(ili)%idproc(i)=idtemp
       SELF_PAR(idprocmat)%idproc=idtemp
!	
!-- Recherche des parametres  associes au processus en cours 
       iposini=index(nomproc_tempo,'(') !--position ouverture parenthese
       iposfin=index(nomproc_tempo,')') !--position fermeture parenthese
       tempo=nomproc_tempo(iposini+1:iposfin-1)
       inbelmt2=f_nschain(tempo,'>')!--nb de param. pour le proc. en cours
!
!-- on verifie que ce nb de param. est correct
       if (inbelmt2/=PROC_MOD(idtemp)%nbpar) then
          write(*,*) 'nb de parametres du processus ',PROC_MOD(idtemp)%nomproc,&
&     'dans la cellule ',ili,' est incoherent'
          write(22,*)'nb de parametres du processus ',PROC_MOD(idtemp)%nomproc,&
&     'dans la cellule ',ili,' est incoherent'
#ifdef MODTEST
       read(*,*)
#endif      
        endif       
        If (associated(SELF_PAR(idprocmat)%valpar)) NULLIFY(SELF_PAR(idprocmat)%valpar)
        Allocate(SELF_PAR(idprocmat)%valpar(inbelmt2),STAT=istat)
        write(11,*) 'alloc de SELF_PAR(',idprocmat,')%valpar a', inbelmt2
        if (istat /= 0) write(*,*) 'pb d''allocation de SELF_PAR%valpar'
          do k= 1, inbelmt2
             tempopo=f_chain(tempo,k,'>')
             Read(tempopo,*)SELF_PAR(idprocmat)%valpar(k)
          enddo
    enddo    
 enddo
!	  
!-- on verifie que le nb cumule de processus est correct
 if (idprocmat /= nbproctot) STOP 'nb de processus cumules de la matrice < a celui indique'

!-- Lecture du rapport Chlorophylle/Carbone:
errlec = 0
do   
      Read(15,*,iostat=errlec) chaine
#ifdef MODTEST
     write(*,*) chaine     
#endif 
 if (errlec /=0) stop 'pb de lecture du fichier de config,rapport Chl/C'
     if (chaine(1:1)=='#') cycle
     tempo=f_chain(chaine,1,':')
     Read(tempo,*) CHL_C_VAL
     if ( CHL_C_VAL == 0) then  ! rapport Chl/C est cst
        CHL_C_BOOL = .FALSE.
        tempo=f_chain(chaine,2,':')
        Read(tempo,*)CHL_C0 
     elseif ( CHL_C_VAL == 1) then
       CHL_C_BOOL = .TRUE.
      else
       stop 'Probleme avec le rapport CHL:C'
      endif
      exit
 enddo
 
!-- lecture du nb de points de sauvegarde des flux
!   au cours du temps le cas echeant

!--initialisation � zero par d�faut:
 errlec = 0
 nbpoints = 0
 do   
      Read(15,*,iostat=errlec) chaine
      if (chaine(1:1)=='#' .and. errlec == 0) cycle
      if (errlec /= 0 ) exit
      read(chaine,*) nbpoints
#ifdef MODTEST
     write(*,*) 'nb de points de sauvegarde des flux =',nbpoints
#endif  
     exit	 
 enddo 
 
!-- lecture des coordonnes de sauvegarde des flux
!   au cours du temps le cas echeant

!-- CAS OU IL Y A DES FLUX A SAUVEGARDER
if (nbpoints > 0) then
 Allocate (coord_flux(nbpoints,3))  ! tableau des coordonnees des points
      write(11,*) 'alloc de coord_flux a',nbpoints,3
 Allocate (noms_fich_flux(nbpoints))! tableau des noms de fichier
      write(11,*) 'alloc de noms_fich_flux a',nbpoints
 
 do ii=1,nbpoints
 do  
      Read(15,*,iostat=errlec) chaine
#ifdef MODTEST
     write(*,*) chaine
#endif 
      if (chaine(1:1)=='#') cycle
      tempo=f_chain(chaine,1,':')
#ifdef MODTEST
      write(*,*) 'tempo =',tempo
#endif 
      Read(tempo,*) i_st
!
      tempo=f_chain(chaine,2,':')
#ifdef MODTEST
      write(*,*) 'tempo =',tempo
#endif 
      Read(tempo,*) j_st
!
      tempo=f_chain(chaine,3,':')
#ifdef MODTEST
      write(*,*) 'tempo =',tempo
#endif 
      Read(tempo,*) k_st
      
!-- Affectation des coordonn�es:
      coord_flux(ii,1)=i_st
      coord_flux(ii,2)=j_st
      coord_flux(ii,3)=k_st

!-- Creation des noms de fichier de sauvegarde :
      c_ist = f_Int2chain(i_st,lon_cist)
      c_jst = f_Int2chain(j_st,lon_cjst)
      c_kst = f_Int2chain(k_st,lon_ckst)


      noms_fich_flux(ii)='flux_'//c_ist(1:lon_cist)//'_'//c_jst(1:lon_cjst)//'_'//c_kst(1:lon_ckst)

!-- Test sur les coordonnees:
      
#ifdef NCOUPL
      if (i_st<nx_min ) then 
        write(*,*) 'l indice de maille selon x pour la sauvegarde des flux est inferieur a',nx_min
	stop
      elseif (i_st>nx_max) then
        write(*,*) 'l indice de maille selon x pour la sauvegarde des flux est superieur a',nx_max
	stop
      elseif (j_st<ny_min ) then 
        write(*,*) 'l indice de maille selon y pour la sauvegarde des flux est inferieur a',ny_min
        stop  
      elseif (j_st>ny_max) then
        write(*,*)  'l indice de maille selon y pour la sauvegarde des flux est superieur a',ny_max
        stop
      elseif (k_st<1 )then
        write(*,*)   'l indice de maille selon z pour la sauvegarde des flux est inferieur a 1'
        stop
      elseif (k_st>nz_max) then
        write(*,*)    'l indice de maille selon z pour la sauvegarde des flux est superieur a',nz_max
        stop
      endif
#endif
      exit
   enddo
 enddo
endif 

!-- Fermeture du fichier
 close(15)

 Call Sub_init_misc
 Call Sub_init_zoo

#ifdef MODTEST
 Write(*,*) ' -- Sub_init termin�e ! -- '
 Write(*,*) ' -- nbcallbio = ',nbcallbio
! Write(*,*) 
#endif

End Subroutine Sub_init
!------------------------------------------------------------------------------
!
                       SUBROUTINE Sub_init_misc
!		       
!-- identifie le nb d organismes ds les compartiments zoo, phy et bac et les
!   indices de variable correspondants
!------------------------------------------------------------------------------
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_USER
Use VAR_GLOBAL

Implicit None

!-- variables locales
Integer :: ivar, iorg, iorgold, iorgold_z,iorgold_b
Integer :: ii,jj,kk,i


!-- Determination du nb de sous-compartiments de phyto et de celui de
!   zoo :
nscp_phy = 0
iorgold = 0
nscp_zoo = 0
iorgold_z = 0
nscp_bac = 0
iorgold_b = 0

do ivar = 1,nbvar
      iorg = VAR(ivar)%idorg
      if (VAR(ivar)%comp == 'phy' .AND. iorg /= iorgold ) then
        nscp_phy = nscp_phy + 1
        iorgold = iorg
      elseif (VAR(ivar)%comp == 'zoo' .AND. iorg /= iorgold_z ) then
        nscp_zoo = nscp_zoo + 1
        iorgold_z = iorg
      elseif (VAR(ivar)%comp == 'bac' .AND. iorg /= iorgold_b ) then
        nscp_bac = nscp_bac + 1
        iorgold_b = iorg
      endif
enddo

!-- Allocation du vecteur contenant les organismes phyto:
if (nscp_phy /=0)  then
    Allocate(iscp_phy(nscp_phy))
      write(11,*) 'alloc de iscp_phy a',nscp_phy
    iscp_phy = 0
    iorgold = 0
endif

if (nscp_zoo /=0)    then
    Allocate(iscp_zoo(nscp_zoo))
      write(11,*) 'alloc de iscp_zoo a',nscp_zoo
    iscp_zoo = 0
    iorgold_z = 0
endif

if (nscp_bac /=0)   then
   Allocate(iscp_bac(nscp_bac))
      write(11,*) 'alloc de iscp_bac a',nscp_bac
   iscp_bac = 0
   iorgold_b = 0
endif
      
!-- Enregistrement des indices d organismes phyto et zoo
ii=1     
jj=1  
kk=1    

      
do ivar = 1,nbvar
      iorg = VAR(ivar)%idorg
      if (VAR(ivar)%comp == 'phy' .AND. iorg /= iorgold ) then
        iscp_phy(ii) = iorg
        ii=ii+1
        iorgold = iorg
      elseif (VAR(ivar)%comp == 'zoo' .AND. iorg /= iorgold_z ) then
        iscp_zoo(jj) = iorg
        jj=jj+1
        iorgold_z = iorg
      elseif (VAR(ivar)%comp == 'bac' .AND. iorg /= iorgold_b ) then
        iscp_bac(kk) = iorg
        kk=kk+1
        iorgold_b = iorg
      endif
enddo

#ifdef MODTEST
 Write(*,*) ' -- Sub_init_misc termin�e  ! -- '
#endif
End Subroutine Sub_init_misc
!------------------------------------------------------------------------------
! Routine d initialisation des variables utilis�es pour simuler le zoo
!------------------------------------------------------------------------------
 subroutine sub_init_zoo()

!************************************
!tableau contenant les indices des biomasses carbon�es dans VAR
!tableau contenant les indices des biomasses azot�es dans VAR

Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use VAR_USER
Use VAR_GLOBAL

implicit none
!!!!!!!!!!!!!!!!modif marion 22/04/10 AJOUT P
!-- variables locales
Integer :: nb_C, nb_N,nb_P,ivar
integer, allocatable :: tempo_C(:)        !tableau contenant les indices des biomasses carbon�es dans VAR
integer, allocatable :: tempo_N(:)
integer, allocatable :: tempo_P(:)

!Allocation des tableaux temporaires referecant les biomasses en C et N
! afin de pouvoir allouer les tableaux finaux sans surdimensionner
allocate(tempo_C(nbvar))
allocate(tempo_N(nbvar))
allocate(tempo_P(nbvar))

nb_C=0
nb_N=0
nb_P=0

do ivar=1, nbvar
 
  if( VAR(ivar)%comp=='phy'.OR.VAR(ivar)%comp=='bac'.OR.VAR(ivar)%comp=='det') then
   
     select case (VAR(ivar)%elmt)
	   case("C")
	     nb_C=nb_C+1
		 tempo_C(nb_C) = ivar
       case("N")
	     nb_N=nb_N+1
		 tempo_N(nb_N) = ivar
      case("P")
             nb_P=nb_P+1
                tempo_P(nb_P) = ivar
	 end select

  end if   
end do
! Allocation et remplissage des tableaux finaux
allocate(BIO_C(nb_C))
allocate(BIO_N(nb_N))
allocate(BIO_P(nb_P))
allocate(CBIO_C(nb_C))
allocate(CBIO_N(nb_N))
allocate(CBIO_P(nb_P))

do ivar=1,nb_C
  BIO_C(ivar) = tempo_C(ivar) 
end do
do ivar=1,nb_N
BIO_N(ivar) = tempo_N(ivar)
end do
do ivar=1,nb_P
BIO_P(ivar) = tempo_P(ivar)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ATTENTION A DECOMMENTER POUR LE MODELE THESE DE VINCENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CBIO_C(:) = 1.
!CBIO_C(3) = 1.
CBIO_N(:) = 1.
!CBIO_N(3) = 1.
CBIO_P(:) = 1.
!!!!!!!!!!!!fin marion

 end subroutine sub_init_zoo
