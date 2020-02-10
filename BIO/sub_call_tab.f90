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
            Subroutine Sub_call_tab ()
!	    
!
!-- Routine de création du fichier "call.inc"
!
!   Création du fichier "call.inc" à partir de la matrice FLUX_PAR.
!   Ce fichier contient les appels aux fonctions modélisant les processus.
!   Il sera inclus dans la boucle temporelle du programme principal.
!   
!   Ces appels permettent de calculer les flux de FLUX_VAL(i,j)
!   et SELF_VAL(i)
!
! Dernière modification: 27/06/07
!-------------------------------------------------------------------------------
 Use DEF_TYPE
 Use MOD_FCHAIN
 Use COUPLEUR_PHY_BIO
 Use VAR_GLOBAL
 Implicit None
   
!--Variables locales
 
 Integer             :: errlec,nb_fonc_cell,idtemp
 Integer             :: ili,jcol,ii,iorg,j,iitest
 integer(2)          :: kk,dim_param
 character(L_CHAIN)  :: c_ili,c_jcol,ichain,ichain2,c_iorg
 character           :: c_jcol_fict
 character(L_CHAIN)  :: nomsub
 character(L_VAR)    :: signe_facteur
 character(L_VAR)    :: signe2
 Integer             :: lon_c_ili,lon_c_jcol,lon_ichain,lon_ichain2,lon_c_iorg
 Integer             :: size_sign2
! 
!-- Ouverture (ou création) du fichier destine a l''appel des fonctions 
!   de processus pour le calcul des flux

 Open(19,FILE="call.inc")
!
 write(19,*) '!-- Ce fichier est généré automatiquement lors de l execution du programme'
 write(19,*) '!-- a partir des informations du fichier config.ini   --'
 write(19,*) '!-- Il contient les expressions necessaires au calcul des termes de --'
 write(19,*) '!-- production primaire (utiles pour le calcul des flux) --'
 write(19,*) '!-- ainsi que l''appel des fonctions permettant le calcul des flux entre --'
 write(19,*) '!-- variables d''etat, et  celui des flux n''impliquant qu''une seule  --'
 write(19,*) '!-- variable d''etat                                                   --'

!------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------ 
 write(19,*)
 write(19,*) '!--- 1/ Expression des flux (non nuls) entre variables d''état' 
 write(19,*) 
!------------------------------------------------------------------------------ 
 
!-- Boucle sur la matrice FLUX_PAR pour les appels des fonctions
!-- permettant de calculer les éléments de la matrice FLUX_VAL(ili,jcol)

ii=1 ! initialisation
do while (ii <= size(FLUX_PAR))
!-- Transformation en chaines de caractères des entiers impliques 
!-- dans les appels des processus (ie les entiers ii, ili, jcol)
! -- entier ii <--> chaine ichain  de longueur lon_ichain
    ichain=f_Int2chain(ii,lon_ichain)

!-- entier ili <--> chaine c_ili de longueur utile lon_c_ili
    ili=FLUX_PAR(ii)%ipos(1)
    c_ili = f_Int2chain(ili,lon_c_ili)

!-- entier jcol <--> chaine c_jcol de longueur utile lon_c_jcol
    jcol=FLUX_PAR(ii)%ipos(2)
    c_jcol=f_Int2chain(jcol,lon_c_jcol)
  
!--Appel pour la matrice FLUX_VAL (NB: le triplet (i,j,k) correspond aux
!--indices du maillage spatial)
 
  write(19,*)'FLUX_VAL(',c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol),')%val = &'
                 
  do nb_fonc_cell=0,size(FLUX_VAL(ili,jcol)%idproc)-1

!-- recuperation du nom de la fction appelée par le processus n°(ii+nb_fonc_cell)
      write(nomsub,*)PROC_MOD(FLUX_PAR(ii+nb_fonc_cell)%idproc)%nomsub
      nomsub=trim(adjustl(nomsub))

!-----------------------------------------------
!---- cas d''une fonction avec parametres ------
!-----------------------------------------------
    if (size(FLUX_PAR(ii+nb_fonc_cell)%valpar )/=0 ) then
! dimension du vecteur des parametres :      
         dim_param=size(FLUX_PAR(ii+nb_fonc_cell)%valpar)
         signe_facteur = TRIM(ADJUSTL(FLUX_PAR(ii+nb_fonc_cell)%signe))

         if (signe_facteur(1:1) /= 'x') then  
           if (len_trim(signe_facteur)>1) then   
               write(19,'(a1,a8,a5)') '(',signe_facteur, ') * &'
           else
               write(19,'(a,a1,a6)') '(',signe_facteur(1:1), '1) * &'
           endif 
         endif
         write(19,*) '  ',nomsub(1:len_trim(nomsub)),'(', &
         c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ', &'
         
!suite appel : un argument par ligne, à l''aide de l''opérateur &
!-- chaine de caracteres pour l''entier ii + nb_fonc_cell:
         ichain2=f_Int2chain(ii+nb_fonc_cell,lon_ichain2)
         do kk=1,dim_param-1
            write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
           ')%valpar(',  kk  ,')', ', & '
         enddo
!-- Test: y-a-t-il d''autres processus dans la meme cellule?
         if (nb_fonc_cell < size(FLUX_VAL(ili,jcol)%idproc)-1) then
!
            signe2 = TRIM(ADJUSTL(FLUX_PAR(ii + nb_fonc_cell + 1 )%signe)) 
            size_sign2 = LEN_TRIM(signe2)
            if (signe2(1:1) =='x') then 
               write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
               ')%valpar(',  dim_param  ,') ) &'	
              if (size_sign2 > 1) then	
                write(19,'(a1,a8,a5)') '*',signe2(2:size_sign2),' * &'	
              else
                write(19,'(A3)') '* &'
              endif
!		
            else
               write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
              ')%valpar(',  dim_param  ,') ) &'
               write(19,'(A3)') '+ &'
            endif
         else
            write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
          ')%valpar(',  dim_param  ,') )'
         endif
 !-----------------------------------------------
!---- cas d''une fonction sans parametres ------
!-----------------------------------------------
   else ! cas d''une fction de processus sans arguments autres que ili,jcol
         dim_param=0
         signe_facteur = TRIM(ADJUSTL(FLUX_PAR(ii+nb_fonc_cell)%signe))
         if (signe_facteur(1:1)/='x') then 
            if (len_trim(signe_facteur)>1) then
              write(19,'(a1,a8,a5)') '(',signe_facteur, ') * &'
            else
              write(19,'(a,a1,a6)') '(',FLUX_PAR(ii+nb_fonc_cell)%signe, '1) * &'
            endif
         endif
    if (nb_fonc_cell < size(FLUX_VAL(ili,jcol)%idproc)-1) then
         signe_facteur = FLUX_PAR(ii+1+nb_fonc_cell)%signe
         if  (signe_facteur(1:1) =='x') then 
              write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(', &
              c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ') &'
               if (len_trim(signe_facteur) > 1) then
                 write(19,'(A1,A8,A3)') '*',signe_facteur(2:len_trim(signe_facteur)),'* &'
               else
                 write(19,'(A3)') '* &'
               endif
         else
             write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(', &
              c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ') &'
            write(19,'(A3)') '+ &'
         endif
    else
            write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(', &
           c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ') '
    endif
endif
  enddo
  write(19,*)' '
  write(19,*)' '
  ii=ii+ size(FLUX_VAL(ili,jcol)%idproc) 
enddo
!------------------------------------------------------------------------------
Write(19,*)
write(19,*) '!--- 2/ Expressions des flux (non nuls) n''impliquant qu''une ' 
write(19,*) '!       seule variable d''état' 
write(19,*)
 
!-- Boucle sur la matrice SELF_PAR pour les appels des fonctions
!-- permettant de calculer les éléments de la matrice SELF_VAL(ii)

ii=1 ! initialisation

do while (ii <= size(SELF_PAR))

!-- Mise sous forme de chaine de caractères des entiers impliques dans les
!--  appels des processus(variable ii, ili, jcol)
  ichain=f_Int2chain(ii,lon_ichain)

!-- transformation de ili en chaine de caracteres c_ili de longueur 
!-- utile lon_c_ili  
  ili=SELF_PAR(ii)%ipos(1) 
  c_ili = f_Int2chain(ili,lon_c_ili)

!-- utilisation fictive de jcol
  c_jcol_fict = '0'
  
!Appel pour la matrice SELF_VAL
  write(19,*) 'SELF_VAL(',c_ili(:lon_c_ili),')%val = &'
                 
  do nb_fonc_cell=0,size(SELF_VAL(ili)%idproc)-1

!-- recuperation du nom de la fction appelée par le processus n°(i+nb_fonc_cell)
    write(nomsub,*)PROC_MOD(SELF_PAR(ii+nb_fonc_cell)%idproc)%nomsub
    nomsub=trim(adjustl(nomsub))
!
!-- Si la fonction impliquee par le processus a des parametres:
!
    if( size(SELF_PAR(ii+nb_fonc_cell)%valpar )/=0 ) then
! dimension du vecteur des parametres	
       dim_param=size(SELF_PAR(ii+nb_fonc_cell)%valpar) 
       signe_facteur = TRIM(ADJUSTL(SELF_PAR(ii+nb_fonc_cell)%signe))
        if (signe_facteur(1:1) /= 'x') then 
            if (LEN_TRIM(signe_facteur)>1) then
                write(19,'(a,a8,a5)') '(',signe_facteur, ') * &'
            else
                write(19,'(a,a1,a6)') '(',SELF_PAR(ii+nb_fonc_cell)%signe, '1) * &'
            endif
       endif  
       write(19,*) '  ',nomsub(1:len_trim(nomsub)),'(',c_ili(:lon_c_ili),',',c_jcol_fict,', &'
!
!suite appel : un argument par ligne, à l'aide de l'opérateur &
! entier ii + nb_fonc_cell <--> ichain2 de long. utile lon_ichain2
       ichain2=f_Int2chain(ii+nb_fonc_cell,lon_ichain2)

       do kk=1,dim_param-1
          write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
         ')%valpar(',  kk  ,')', ', & '
       enddo
!Autres processus ds la cellule?
       if(nb_fonc_cell < size(SELF_VAL(ili)%idproc)-1) then 
!
          signe2 = TRIM(ADJUSTL(SELF_PAR(ii + nb_fonc_cell + 1 )%signe))
          size_sign2 = LEN_TRIM(signe2)
!          if (SELF_PAR(ii + nb_fonc_cell + 1 )%signe =='x') then 
          if (signe2(1:1)== 'x') then
             write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
             ')%valpar(',  dim_param  ,') ) &'	
             if (size_sign2 > 1 ) then
                write(19,'(A1,A8,A3)') '*',signe2(2:size_sign2),'* &'	
             else
                write(19,'(A4)') '* &'	
             endif
          else
              write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
            ')%valpar(',  dim_param  ,')', ') & '
             write(19,'(A3)') '+ &'
          endif
       else
          write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
          ')%valpar(',  dim_param  ,')', ') '
       endif
    else ! cas d''une fonction ou sub de processus sans arguments
       dim_param=0
       signe2 = TRIM(ADJUSTL(SELF_PAR(ii+nb_fonc_cell)%signe))
       size_sign2 = LEN_TRIM(signe2)
       if (signe2(1:1)/= 'x') then
           if (size_sign2 > 1 ) then
               write(19,'(a1,a8,a6)') '(',signe2, ') * &'
           else
               write(19,'(a,a1,a6)') '(',SELF_PAR(ii+nb_fonc_cell)%signe, '1) * &'
           endif
       endif
       if(nb_fonc_cell < size(SELF_VAL(ili)%idproc)-1) then ! d''autres processus ds la cellule
        signe2 = TRIM(ADJUSTL(SELF_PAR(ii+1+nb_fonc_cell)%signe))
        size_sign2 = LEN_TRIM(signe2)
         if (signe2(1:1) /= 'x') then  
            write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(',&
            c_ili(:lon_c_ili),',',c_jcol_fict,') &'
            write(19,'(A3)') '+ &'
         else
            write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(',&
            c_ili(:lon_c_ili),',',c_jcol_fict,') &'
             if (size_sign2 > 1 ) then
                write(19,'(A1,A8,A3)') '*',signe2(2:size_sign2),'* &'	
             else
                write(19,'(A4)') '* &'	
             endif           
         endif
        else !pas d''autres processus ds la cellule
          write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(',&
          c_ili(:lon_c_ili),',',c_jcol_fict,')'
        endif 
    endif   
  enddo
  write(19,*)' '
  write(19,*)' '
  ii=ii+ size(SELF_VAL(ili)%idproc) 
enddo

 !-- Fermeture du fichier "call.inc"
 close(19)
      
         print * , " -- Sub_call_tab Terminé ! -- "
 
End Subroutine Sub_call_tab
!-------------------------------------------------------------------------------
!
            Subroutine Sub_calc_extinc ()
!	    
!-------------------------------------------------------------------------------
 Use DEF_TYPE
 Use MOD_FCHAIN
 Use COUPLEUR_PHY_BIO
 Use VAR_GLOBAL

 Implicit None
!--Variables locales
 
 Integer             :: errlec,nb_fonc_cell,idtemp
 Integer             :: ili,jcol,ii,iorg,j,iitest
 integer(2)          :: kk,dim_param
 character(L_CHAIN)  :: c_ili,c_jcol,ichain,ichain2,c_iorg
 character           :: c_jcol_fict
 character(L_CHAIN)  :: nomsub
 character(L_VAR)    :: signe_facteur
 character(L_VAR)    :: signe2
 Integer             :: lon_c_ili,lon_c_jcol,lon_ichain,lon_ichain2,lon_c_iorg
 Integer             :: size_sign2

 Open(19,FILE="calc_extinc.inc")
!
 write(19,*) '!-- Ce fichier est genere automatiquement lors de l execution du     --'
 write(19,*) '!-- programme a partir des informations du fichier config.ini        --'
 write(19,*) '!-- Il contient les expressions necessaires  au calcul du terme      --'
 write(19,*) '!-- exponentiel d extinction de la lumiere dans l eau                --'
        
 write(19,*) 
 idtemp = IRR_PAR%idproc
 write(nomsub,*) PROC_MOD(idtemp)%nomsub
 nomsub=trim(adjustl(nomsub))
 dim_param = PROC_MOD(idtemp)%nbpar
 
 write(19,'(A12)')'kextinc = &'
 write(19,*) nomsub(1:len_trim(nomsub)),'(& '
  do kk=1,dim_param
   write(19,*)'     IRR_PAR%valpar(',  kk  ,'),  & '
  enddo
  write(19,*) 'Chl_tot,dz) ' !derniers arguments de la fonction
!-- Fermeture du fichier "calc_extinc.inc"
 close(19)
         print * , " -- Sub_calc_extinc ! -- "
 
End Subroutine Sub_calc_extinc
!-------------------------------------------------------------------------------
!
                   SUBROUTINE Sub_save_flux     
!	    
!
!-- Routine de création du fichier "save_flux.inc"
!
!   Création du fichier "save_flux.inc" à partir de la matrice FLUX_PAR.
!   Ce fichier contient les instructions pour pouvoir sauvegarder les 
!   differentes composantes d''un flux donne entre deux variables d''etat
!   du modele FLUX_VAL(i,j), ou entre une variable d''etat et une variable 
!   non representee: SELF_VAL(i)
!
! Dernière modification: 27/06/07
!
!-------------------------------------------------------------------------------
 Use DEF_TYPE
 Use MOD_FCHAIN
 Use COUPLEUR_PHY_BIO
 Use VAR_GLOBAL
 Use VAR_USER
 
 Implicit None
 
!--Variables locales
 integer            :: errlec,idtemp 
 integer            :: nb_fonc_cell
 integer            :: ili,jcol,ii,jj,iii,iorg,j,iitest
 integer(2)         :: kk,dim_param
 character(L_CHAIN) :: c_ili,c_jcol,ichain,ichain2,c_iorg
 character(L_CHAIN) :: c_ist,c_jst,c_kst
 integer            :: i_st,j_st,k_st   
 integer            :: lon_cist, lon_cjst, lon_ckst
 integer            :: lon_c_ili,lon_c_jcol,lon_ichain,lon_ichain2,lon_c_iorg
 integer            :: size_sign2
 integer            :: nb_fonc_max
 integer            :: numfich
  
 character          :: c_jcol_fict
 character(L_CHAIN) :: nomsub
 character(L_VAR)   :: signe_facteur
 character(L_VAR)   :: signe2
 character (20)     :: fich
! 

!-- 1/ recherche du flux pour lequel il y a un max de fonctions differentes
!   afin de fixer le nb de colonnes:

      nb_fonc_max = 0
      do ili=1,nbvar
       do jcol=1,nbvar
          if(ASSOCIATED(FLUX_VAL(ili,jcol)%idproc)) then
            if (size(FLUX_VAL(ili,jcol)%idproc) > nb_fonc_max) nb_fonc_max = size(FLUX_VAL(ili,jcol)%idproc)
	  endif
	enddo
       enddo
   
       do ili=1,nbvar
          if(ASSOCIATED(SELF_VAL(ili)%idproc)) then
            if (size(SELF_VAL(ili)%idproc) > nb_fonc_max) nb_fonc_max = size(SELF_VAL(ili)%idproc)
	  endif
	enddo
   
!-- 2/ Ouverture du fichier dans lequel figureront les instructions pour la sauvegarde des flux
 open(19, File = './BIO/call_save_flux.inc')
 write(19,*) '!-- Ce fichier est genere automatiquement lors de l execution du programme'
 write(19,*) '!-- a partir des informations du fichier config.ini   --'
 write(19,*) '!-- Il contient les instructions necessaires a la sauvegarde --'
 write(19,*) '!-- des differentes composantes d''un flux entre deux variables --'
 write(19,*) ' '
 write(19,*) ' '


!----------       Boucle sur la matrice FLUX_PAR pour les appels des fonctions           ----------------
!----------   permettant de calculer les éléments de la matrice FLUX_VAL(ili,jcol)       ----------------

ii=1 ! initialisation

do while (ii <= size(FLUX_PAR))
!-- Transformation en chaines de caractères des entiers impliques 
!-- dans les appels des processus (ie les entiers ii, ili, jcol)
! -- entier ii <--> chaine ichain  de longueur lon_ichain
    ichain=f_Int2chain(ii,lon_ichain)

!-- entier ili <--> chaine c_ili de longueur utile lon_c_ili
    ili=FLUX_PAR(ii)%ipos(1)
    c_ili = f_Int2chain(ili,lon_c_ili)

!-- entier jcol <--> chaine c_jcol de longueur utile lon_c_jcol
    jcol=FLUX_PAR(ii)%ipos(2)
    c_jcol=f_Int2chain(jcol,lon_c_jcol)
                    
!  Allocation de ss_flux:
write(19,*) 'dimfonc = size(FLUX_VAL(', c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol) ,')%idproc)'
write(19,*) 'Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))'
write(19,*) ' '
write(19,*) ' '

!----------------------------------------------------------------------
! Debut de la boucle sur les fonctions impliquees dans un flux donne
!----------------------------------------------------------------------
  do nb_fonc_cell = 0 , size(FLUX_VAL(ili,jcol)%idproc)-1

!-- Recuperation du nom de la fction appelée par le processus n°(ii+nb_fonc_cell)
  write(nomsub,*)PROC_MOD(FLUX_PAR(ii+nb_fonc_cell)%idproc)%nomsub
  nomsub=trim(adjustl(nomsub))
 
!             -----------------------------------------------
!             ---- cas d''une fonction avec parametres ------
!             -----------------------------------------------
  if (size(FLUX_PAR(ii+nb_fonc_cell)%valpar )/=0 ) then
    
! Dimension du vecteur des parametres :      
!------------------------------------
      dim_param     = size(FLUX_PAR(ii+nb_fonc_cell)%valpar)
      signe_facteur = TRIM(ADJUSTL(FLUX_PAR(ii+nb_fonc_cell)%signe))
!
      if (signe_facteur(1:1) /= 'x') then  
         if (len_trim(signe_facteur)>1) then   
            write(19,*) ' '
            write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =  &'
            write(19,'(a1,a8,a5)') '(',signe_facteur, ') * &'
         else
            write(19,*) ' '
            write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =   &'
            write(19,'(a,a1,a6)') '(',signe_facteur(1:1), '1) * &'
         endif 
      endif
      write(19,*) '',nomsub(1:len_trim(nomsub)),'(', &
                    c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ', &'
         
!-- Suite appel : un argument par ligne, à l''aide de l''opérateur &
!-- chaine de caracteres pour l''entier ii + nb_fonc_cell:

         ichain2=f_Int2chain(ii+nb_fonc_cell,lon_ichain2)
         do kk=1,dim_param-1
            write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
           ')%valpar(',  kk  ,')', ', & '
         enddo
	 
!-- Test: y-a-t-il d''autres fonctions  dans la meme cellule?
!----------------------------------------------------------
!
     if (nb_fonc_cell < size(FLUX_VAL(ili,jcol)%idproc)-1) then
!
            signe2 = TRIM(ADJUSTL(FLUX_PAR(ii + nb_fonc_cell + 1 )%signe)) 
            size_sign2 = LEN_TRIM(signe2)
            if (signe2(1:1) =='x') then 
               write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
               ')%valpar(',  dim_param  ,') ) '	
               if (size_sign2 > 1) then	
                 write(19,*) ' '
                 write(19,*) 'ss_flux(',nb_fonc_cell+2,',:,:,:) =  &'
                 write(19,'(a1,a8,a5)') ',',signe2(2:size_sign2),' * &'	
               else
                 write(19,*) ' '
                 write(19,*) 'ss_flux(',nb_fonc_cell+2,',:,:,:) =  &'
               endif		
            else !-- cas où on va commencer un nouveau processus a la prochaine iteration	     
               write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
              ')%valpar(',  dim_param  ,') )'
            endif	    
     elseif (nb_fonc_cell == size(FLUX_VAL(ili,jcol)%idproc)-1) then
            write(19,*) '     FLUX_PAR(',  ichain2(:lon_ichain2)  ,&
          ')%valpar(',  dim_param  ,') )'
     endif
!              -----------------------------------------------------------------
  else !        ---- cas d''une fonction sans arguments autres que ili,jcol ------
!              -----------------------------------------------------------------
         dim_param=0
         signe_facteur = TRIM(ADJUSTL(FLUX_PAR(ii+nb_fonc_cell)%signe))
         if (signe_facteur(1:1)/='x') then 
            if (len_trim(signe_facteur)>1) then
              write(19,*) ' '
              write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =  &'
              write(19,'(a1,a8,a5)') '(',signe_facteur, ') * &'
            else
              write(19,*) ' '
              write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =  &'
              write(19,'(a,a1,a6)') '(',FLUX_PAR(ii+nb_fonc_cell)%signe, '1) * &'
            endif
         endif
!-- Test: y-a-t-il d''autres fonctions  dans la meme cellule?
!----------------------------------------------------------
     if (nb_fonc_cell < size(FLUX_VAL(ili,jcol)%idproc)-1) then
         !signe_facteur = FLUX_PAR(ii+1+nb_fonc_cell)%signe
         signe2 = TRIM(ADJUSTL(FLUX_PAR(ii + nb_fonc_cell + 1 )%signe)) 
         size_sign2 = LEN_TRIM(signe2)
         if (signe2(1:1) =='x') then 
             write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(', &
              c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ')'
                if (size_sign2 > 1) then	
                 write(19,'(A1,A8,A3)') '*',signe_facteur(2:len_trim(signe_facteur)),'* &'
               else
                 write(19,*) ' '
                 write(19,*) 'ss_flux(',nb_fonc_cell+2,',:,:,:) =  &'
               endif
         else!--cas d''un nouveau processus à la prochaine iteration
             write(19,*)' ',nomsub(1:len_trim(nomsub)),'(', &
              c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ') '
            write(19,'(A3)') ''
         endif
     elseif (nb_fonc_cell == size(FLUX_VAL(ili,jcol)%idproc)-1 ) then
            write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(', &
           c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol), ') '
     endif
  endif
     
  enddo !fin de la boucle sur les fonctions de processus

!-----------------
!-- Ecriture dans les fichiers de resultats :
!-----------------
do iii=1,size(noms_fich_flux)
 
!-- coordonnees des points où sauver les valeurs des flux :
  i_st = coord_flux(iii,1)
  j_st = coord_flux(iii,2)
  k_st = coord_flux(iii,3)
   
!-- transformation des coordonnees du point de stockage en
!   chaines de caractère:
  c_ist = f_Int2chain(i_st,lon_cist)
  c_jst = f_Int2chain(j_st,lon_cjst)
  c_kst = f_Int2chain(k_st,lon_ckst)
    
!-- numero du fichier (créé et ouvert dans sub_affiche) où vont s''afficher les
!   valeurs de flux au point de coordonnées i_st, j_st, k_st :
       numfich = 200+coord_flux(iii,1)+coord_flux(iii,2)+coord_flux(iii,3)
      write(19,*) ' '
      write(19,*) ' '
      write(19,*) 'write(',numfich,',''(E10.2,$)'') tps/3600. '
      write(19,*) 'write(',numfich,',''(I4,$)'')', c_ili(:lon_c_ili) 
      write(19,*) 'write(',numfich,',''(I4,$)'')', c_jcol(:lon_c_jcol) 
      do jj=1,size(FLUX_VAL(ili,jcol)%idproc)
           write(19,*) 'write(',numfich,',',char(34),'(E10.2,$)',char(34),') ss_flux(',jj,',',c_ist(:lon_cist),',',c_jst(:lon_cjst),',',c_kst(:lon_ckst),')'
      enddo
      if ( size(FLUX_VAL(ili,jcol)%idproc) < nb_fonc_max) then
         do jj=size(FLUX_VAL(ili,jcol)%idproc)+1,nb_fonc_max
             write(19,*) 'write(',numfich,', ''(A10,$)'' )      '' NaN ''    '
         enddo
      endif
      write(19,*) 'write(',numfich,',''(E12.3E3)'')   FLUX_VAL(',c_ili(:lon_c_ili),',',c_jcol(:lon_c_jcol),')%val(',&
	  c_ist(:lon_cist),',',c_jst(:lon_cjst),',',c_kst(:lon_ckst),')'
enddo
write(19,*) ''
write(19,*) 'deallocate (ss_flux)'
write(19,*) ''
ii=ii+ size(FLUX_VAL(ili,jcol)%idproc) 
enddo

!------------------------------------------------------------------------------
!-------- Boucle sur la matrice SELF_PAR pour les appels des fonctions --------
!-------- permettant de calculer les éléments de la matrice SELF_VAL(ii) ------
!------------------------------------------------------------------------------
Write(19,*)
write(19,*) '!--- 2/ Expressions des flux (non nuls) n''impliquant qu''une ' 
write(19,*) '!       seule variable d''état' 
write(19,*)
 

ii=1 ! initialisation

do while (ii <= size(SELF_PAR))

!-- Mise sous forme de chaine de caractères des entiers impliques dans les
!--  appels des processus(variable ii, ili, jcol)
  ichain=f_Int2chain(ii,lon_ichain)

!-- transformation de ili en chaine de caracteres c_ili de longueur 
!-- utile lon_c_ili  
  ili=SELF_PAR(ii)%ipos(1) 
  c_ili = f_Int2chain(ili,lon_c_ili)
!-- utilisation fictive de jcol
  c_jcol_fict = '0'
 
!  Allocation de ss_flux:
write(19,*)'dimfonc = size(SELF_VAL(', c_ili(:lon_c_ili),')%idproc)'
write(19,*) 'Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))'
               
!----------------------------------------------------------------------
! Debut de la boucle sur les fonctions impliquees dans le flux de ili vers  jcol:
!----------------------------------------------------------------------

do nb_fonc_cell = 0,size(SELF_VAL(ili)%idproc)-1
 
!-- récuperation du nom de la fction appelée par le processus n°(i+nb_fonc_cell)
    write(nomsub,*)PROC_MOD(SELF_PAR(ii+nb_fonc_cell)%idproc)%nomsub
    nomsub=trim(adjustl(nomsub))
!
!       -----------------------------------------------------------
!       ---- cas d''une fonction de processus avec parametres ------
!       -----------------------------------------------------------
!
   if( size(SELF_PAR(ii+nb_fonc_cell)%valpar )/=0 ) then
    
! dimension du vecteur des parametres	
       dim_param=size(SELF_PAR(ii+nb_fonc_cell)%valpar) 
       signe_facteur = TRIM(ADJUSTL(SELF_PAR(ii+nb_fonc_cell)%signe))
        if (signe_facteur(1:1) /= 'x') then 
            if (LEN_TRIM(signe_facteur)>1) then
                write(19,*) ' '
                write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =  &'
                write(19,'(a,a8,a5)') '(',signe_facteur, ') * &'
            else
                write(19,*) ' '
                write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =  &'
                write(19,'(a,a1,a6)') '(',SELF_PAR(ii+nb_fonc_cell)%signe, '1) * &'
            endif
        endif  
        write(19,*) '  ',nomsub(1:len_trim(nomsub)),'(',c_ili(:lon_c_ili),',',c_jcol_fict,', &'
!
!suite appel : un argument par ligne, à l''aide de l''opérateur &
! entier ii + nb_fonc_cell <--> ichain2 de long. utile lon_ichain2
       ichain2=f_Int2chain(ii+nb_fonc_cell,lon_ichain2)

       do kk=1,dim_param-1
          write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
         ')%valpar(',  kk  ,')', ', & '
       enddo
!
! Autres processus ds la cellule?
!--------------------------------
     if(nb_fonc_cell < size(SELF_VAL(ili)%idproc)-1) then 
!
          signe2 = TRIM(ADJUSTL(SELF_PAR(ii + nb_fonc_cell + 1 )%signe))
          size_sign2 = LEN_TRIM(signe2)
          if (signe2(1:1)== 'x') then
             write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
             ')%valpar(',  dim_param  ,') )'	
             if (size_sign2 > 1 ) then
                 write(19,*) ' '
                 write(19,*) 'ss_flux(',nb_fonc_cell+2,',:,:,:) =  &'
                 write(19,'(A1,A8,A3)') '*',signe2(2:size_sign2),'* &'	
             else
                 write(19,*) ' '
                 write(19,*) 'ss_flux(',nb_fonc_cell+2,',:,:,:) =  &'
             endif
          else
              write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
            ')%valpar(',  dim_param  ,') )  '
              write(19,'(A3)') ' '
          endif
     elseif (nb_fonc_cell == size(SELF_VAL(ili)%idproc)-1) then 
          write(19,*) '     SELF_PAR(',  ichain2(:lon_ichain2)  ,&
          ')%valpar(',  dim_param  ,') ) '
     endif
!                -----------------------------------------------
  else !          ---- cas d''une fonction sans parametres ------
!                -----------------------------------------------
         dim_param=0
         signe_facteur = TRIM(ADJUSTL(SELF_PAR(ii+nb_fonc_cell)%signe))
         if (signe_facteur(1:1)/='x') then 
            if (len_trim(signe_facteur)>1) then
              write(19,*) ' '
              write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =  &'
              write(19,'(a1,a8,a5)') '(',signe_facteur, ') * &'
            else
              write(19,*) ' '
              write(19,*) 'ss_flux(',nb_fonc_cell+1,',:,:,:) =  &'
              write(19,'(a,a1,a6)') '(',SELF_PAR(ii+nb_fonc_cell)%signe, '1) * &'
            endif
         endif
!-- Test: y-a-t-il d''autres fonctions  dans la meme cellule?
!----------------------------------------------------------
     if (nb_fonc_cell < size(SELF_VAL(ili)%idproc)-1) then
         signe2 = TRIM(ADJUSTL(SELF_PAR(ii + nb_fonc_cell + 1 )%signe)) 
         size_sign2 = LEN_TRIM(signe2)
         if (signe2(1:1) =='x') then 
             write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(', &
              c_ili(:lon_c_ili),',',c_jcol_fict, ')'
                if (size_sign2 > 1) then	
                 write(19,'(A1,A8,A3)') '*',signe_facteur(2:len_trim(signe_facteur)),'* &'
               else
                 write(19,*) ' '
                 write(19,*) 'ss_flux(',nb_fonc_cell+2,',:,:,:) =  &'
               endif
         else!--cas d''un nouveau processus à la prochaine iteration
             write(19,*)' ',nomsub(1:len_trim(nomsub)),'(', &
              c_ili(:lon_c_ili),',',c_jcol_fict, ') '
            write(19,'(A3)') ''
         endif
     elseif (nb_fonc_cell == size(SELF_VAL(ili)%idproc)-1 ) then
            write(19,*)'  ',nomsub(1:len_trim(nomsub)),'(', &
           c_ili(:lon_c_ili),',',c_jcol_fict, ') '
     endif
  endif

 enddo !-- fin de boucle sur les fonctions

!-----------------
!-- Ecriture dans les fichiers de resultats :
!-----------------
do iii=1,size(noms_fich_flux)
 
!-- coordonnees des points où sauver les valeurs des flux :
  i_st = coord_flux(iii,1)
  j_st = coord_flux(iii,2)
  k_st = coord_flux(iii,3)
   
!-- transformation des coordonnees du point de stockage en
!   chaines de caractère:
  c_ist = f_Int2chain(i_st,lon_cist)
  c_jst = f_Int2chain(j_st,lon_cjst)
  c_kst = f_Int2chain(k_st,lon_ckst)
    
!-- numero du fichier (créé et ouvert dans sub_affiche) où vont s''afficher les
!   valeurs de flux au point de coordonnées i_st, j_st, k_st :
       numfich = 200+coord_flux(iii,1)+coord_flux(iii,2)+coord_flux(iii,3)
      write(19,*) ' '
      write(19,*) ' '
      write(19,*) 'write(',numfich,',''(E10.2,$)'') tps/3600. '
      write(19,*) 'write(',numfich,',''(I4,$)'')', c_ili(:lon_c_ili) 
      write(19,*) 'write(',numfich,',''(I4,$)'')', c_ili(:lon_c_ili) 
      do jj=1,size(SELF_VAL(ili)%idproc)
           write(19,*) 'write(',numfich,',',char(34),'(E10.2,$)',char(34),') ss_flux(',jj,',',c_ist(:lon_cist),',',c_jst(:lon_cjst),',',c_kst(:lon_ckst),')'
      enddo
      if ( size(SELF_VAL(ili)%idproc) < nb_fonc_max) then
         do jj=size(SELF_VAL(ili)%idproc)+1,nb_fonc_max
             write(19,*) 'write(',numfich,', ''(A10,$)'' )      '' NaN ''    '
         enddo
      endif
      write(19,*) 'write(',numfich,',''(E12.3E3)'')   SELF_VAL(',c_ili(:lon_c_ili),')%val(',&
	  c_ist(:lon_cist),',',c_jst(:lon_cjst),',',c_kst(:lon_ckst),')'

enddo
write(19,*) ''
write(19,*) 'ss_flux=0.d0'
write(19,*) 'deallocate(ss_flux)'
write(19,*) ''

!------ nouvel element de SELF_VAL
  ii=ii+ size(SELF_VAL(ili)%idproc) 
enddo
 
 !------------------------------------------------------------------------------
 !-- Fermeture du fichier de sauvegarde des flux
 close(19)
 !------------------------------------------------------------------------------



End Subroutine Sub_save_flux

!-------------------------------------------------------------------------------
