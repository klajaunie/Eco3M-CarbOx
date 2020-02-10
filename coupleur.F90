!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
! m.baklouti@univmed.fr; vincent.faure@univmed.fr
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
!---------------------------------------------------------------------------
!
!                      PROGRAMME PRINCIPAL
!
!  * Ce programme sert d interface de couplage entre Eco3M et un code physique
!    dans le cas où la clé de précompilation COUPL est utilisée
!  * Ce programme est aussi le programme principal en mode non couplé (clé
!     de précompilation NCOUPL), c'est à dire dans le cas d'une utilisation 
!     d'Eco3M sans intervention de l'hydrodynamique.
!  * Le sous-programme MAIN_BIO (fichier main_bio.F90) contient le programme
!    principal du code Eco3M destiné à la modélisation biogéochimique.
!  * Cette version du coupleur est amenée à évoluer de façon à être la plus
!    standart possible pour le couplage à différents codes physique
! 
!---------------------------------------------------------------------------
Use DEF_TYPE
Use VAR_USER
Use VAR_GLOBAL
Use COUPLEUR_PHY_BIO
Use F_PROCESS

Implicit None
#ifdef COUPL
 include './PHY/comms.f90'
#endif

!-- variables locales :
Integer :: imat,jmat,istat
Integer :: nbcall
Integer :: l,n,m,i,j,k

!---------------------------------------------------------------------------------
!--- 1ere execution de Eco3M pour introduire la structure du modele ds le code

#ifdef INI
  #ifdef COUPL  
! ce sous-programme permet de recuperer la dimension du maillage et 
! est a adapter a chaque code physique. 
    call Recup_dim_maill_Mobee ! specifique au couplage avec Mobeehdycs
  #endif
  call MAIN_BIO
  call sub_fin_INI

#endif
!---------------------------------------------------------------------------------

#ifdef CALC

#ifdef MODTEST
  write(11,*) ' '
  write(11,*) ' debut du mode CALC'
  write(11,*) ' '
#endif
!nombre d appel au module bio
      nbcallbio = -1 !-- debute a -1 car 2 etapes d initialisation sont necessaires
!
  write(*,*)'         ------------------------------------'
  write(*,*)'         debut de l  execution en mode calcul'
  write(*,*)'         ------------------------------------'
  call MAIN_BIO
  
#ifdef NCOUPL
  #ifdef M3D-NCOUPL
!-- definition le cas echeant des immersions et profondeurs de mailles
!  lorsqu''aucun code physique n''est couple. A ameliorer pour permettre
!  cette lecture dans un fichier de bathymetrie pour plus de generalite
!
Allocate(prof(nx_min:nx_max,ny_min:ny_max))
Allocate(alti_z(nx_min:nx_max,ny_min:ny_max,nz_max))

do i=nx_min,nx_max
  do j=ny_min,ny_max
        prof(i,j) = 2.0d0 !prof de 2 metres
!DEBUG Christel
write(*,*) i,j, 'prof=',prof(i,j)
    if (pos_nzmax == 'SURF') then 
       alti_z(i,j,nz_max) = prof(i,j) - 1.d0/10.d0
       do k=nz_max-1,1,-1
        alti_z(i,j,k)= alti_z(i,j,k+1) - 1.d0/5.d0 
	enddo
write(*,*) 'alti_z=',alti_z(i,j,k)
    elseif (pos_nzmax == 'FOND') then 
       alti_z(i,j,1) = prof(i,j) - 1.d0/10.d0
       do k=2,nz_max
         alti_z(i,j,k)= alti_z(i,j,k-1) - 1.d0/5.d0 
       enddo
write(*,*)  'alti_z=',alti_z(i,j,k)
    endif
   enddo
 enddo
  #endif
#endif
! 
#ifdef COUPL  
! renvoyer de cet appel de MAIN_BIO les infos au code physique
! concernant le nb de variables bio et le nb d elements ds
! lesquels elles sont exprimees.
!-- lancement de l initialisation de Mecca et en rapporter:
!-- la matrice temperature, la matrice salinite, la taille du
!   maillage, le pas de temps, la duree de simulation, et les
!   concentrations biogeo initiales
!
!-- 1er appel du code physique pour en rapporter (cf aMECCA.for):
!  - les dimensions du maillage (nx_max,ny_max,nz_max) par module
!  - hr1,nstmax par argument
!  - le pas de tps dt_phy par module
!  - les champs de salinite et temperature 
!  - l immersion des mailles
!
!-- initialisation des variables imposees par le code physique
 !-- variables transmises par argument:
!   hr1=0.0    ! heure de la simul si l heure initiale (hr0) n est pas = a zero (h)
    nstmax =0  ! nb d iterations max de la simul
!-- initialisation des variables imposees par l interface de couplage
    tps_hr=0.0 ! nb d heures depuis le debut de la simulation (h)
               ! tps_hr est une variable globale (module COUPL_PHY_BIO)
    nst=0
    nbcall = 0 ! nb d''appels du code physique
!-----------------------------------------------------------------------
!        1er appel : nbcall = 0 
!                 execution partielle pour recuperer des donnees 
!                 (ex. liees au domaine geographique et aux parametres temporels
!                 de la simulation) utiles au modèle bio
       
    call MAIN_PHY(nbcall)
!-----------------------------------------------------------------------
!        2nd appel : nbcall =1  du code physique pour en rapporter les
!                    concentrations biogeochimiques et les temperatures et salinites 
!                    initiales afin de calculer les tendances a l''instant 0

    call ALLOC_VAR_PHY 
    
    call MAIN_PHY(nbcall)
#ifdef MODTEST
    write(*,*) 'retour de MAIN_PHY, nbcall =',nbcall
#endif
!-----------------------------------------------------------------------

 !-- a l instant t=0, les conc prennent leur valeur initiale
   conc_phy = conc_phy0
   conc_phy_old = conc_phy0
#endif

!-- Allocation dynamique des tableaux utilises ds Eco3M (autres 
!  que salinite et temperature) :
 call ALLOC_VAR_Eco3M(nbcall)
 
!-- sauvegarde concentrations ds Mobeehdycs
#ifdef COUPL
   call  TRANSFERT_PHY_Eco3M
#ifdef MODTEST
    write(*,*) 'retour de TRANSFERT_PHY_Eco3M'
#endif
   call DEALLOC_VAR_PHY
#endif    


!-- Initialisation du temps et des parametres associes
!-----------------------------------------------------
tps=0.d0 ! tps (en seconde ecoule depuis le debut de la simul)
#ifdef COUPL
!specifique a mohbeedycs:
tps_hr=0.0
nst=0
#endif
! fin specifique a mohbeedycs
!------------------------------------------------------------------------------------------------
!-----  Début de la BOUCLE DE CALCUL  (nbcall = 1 et nbcallbio = 0 à la 1ere iteration) ---------
!------------------------------------------------------------------------------------------------
do
!
!-- Calcul des tendances a l instant tps par Eco3M
 call MAIN_BIO

#ifdef MODTEST
    write(*,*) 'retour de MAIN_BIO, nbcallbio =',nbcallbio
#endif
    
 !------------------- MODE NON COUPLE    ----------------------------------
  #ifdef NCOUPL
 !-- Boucle de calcul (sur le temps en 0D) :
  tps = tps + dt_bio !-- calcul du tps de simulation en mode 0D
! modif marion pour simuler un upwelling  quand le system est stable 
 !   if ((tps>=17280000).and.(tps<=17286400))then
 !  VAR(15)%conc=VAR(15)%conc+0.03105
 !    TEMP_BIO=14.d0
 !     endif
 !      if (tps>=31104000)then
 !         TEMP_BIO=17.5d0
 !       endif
 !  write(*,*)'TEMP_BIO coupleur=',TEMP_BIO
!fin marion
  if (mod(tps,3600.) < 1d-5 )then
   write(*,*)
   write(*,*)'-------------------------------------------'
   write(*,*) '           tps(h) =',tps/3600
   write(*,*)'-------------------------------------------'
   write(*,*)
  endif
! NB: le test sur le temps pour stopper la boucle d iteration en mode 0D 
!     est realise dans le prog. MAIN_BIO qd tps est > a duree
 #endif
 !------------------- FIN  MODE NON COUPLE    -----------------------------

 !------------------- MODE  COUPLE    -------------------------------------
 #ifdef COUPL
!-- Transfert au code physique  des variables necessaires (concentrations et tendances i.e.
!   termes source - puits)
 call ALLOC_VAR_PHY
#ifdef MODTEST
  write(*,*) 'retour de ALLOC_VAR_PHY'
#endif

!-- transfert des calculs Eco3M au code physique (apres reallocation)
 call  TRANSFERT_Eco3M_PHY
#ifdef MODTEST
  write(*,*) 'retour de transfert eco3m --> phy'
#endif


!-- Desallocation partielle des concentrations et flux Eco3M :
  call DEALLOC_VAR_Eco3M(nbcall)
#ifdef MODTEST
  write(*,*) 'retour de dealloc  eco3m'
#endif

 
!---  Calcul physique en boucle tant que ------
!---  la bio ne doit pas etre calculee   ------
#ifdef MODTEST
  write(*,*) 'avant physique'
#endif
 do
   nst=nst+1
   tps_hr=float(nst)*dt_phy/3600.
   tps=float(nst)*dt_phy
   tps_hr1=tps_hr + tps_hr0
   
!-- Appels du code physique:

  call MAIN_PHY(nbcall)
!
! le transfert depuis la physique des concentrations des variables d etat bio,
! de la salinite, la temperature et la  conc de MES se font par module
   if (mod(nst,500)==0) then
      write(*,*) 'tps(h)',tps/3600
   endif
   if (mod(tps,dt_bio)==0 .or. nst==nstmax) then
      exit
   endif
 enddo
#ifdef MODTEST
  write(*,*) 'apres physique, nbcall =',nbcall
#endif
!--------------- Fin du calcul physique en boucle ----------------
   
 if (nst == nstmax) then
      write(*,*) 'temps final', tps/3600/24.,'(j) atteint'
      call sub_fin !desallocation des tableaux de MAIN_BIO
!      call MAIN_PHY(nst,nstmax,nbcall)! dernier appel de MAIN_PHY
      call DEALLOC_VAR_PHY
      write(*,*) 'FIN DU PROGRAMME'
      exit
 endif

!-- Allocation dynamique des tableaux lies au maillage :
 call ALLOC_VAR_Eco3M(nbcall)
#ifdef MODTEST
  write(*,*) 'retour de ALLOC_VAR_Eco3M'
#endif

!-- sauvegarde concentrations ds Mobeehdycs
 call  TRANSFERT_PHY_Eco3M
#ifdef MODTEST
  write(*,*) 'retour de TRANSFERT_PHY_Eco3M'
#endif
     
!-- deallocation des variables globales du code physique     
 call DEALLOC_VAR_PHY
#ifdef MODTEST
  write(*,*) 'retour de DEALLOC_VAR_PHY'
#endif

#endif 
!-- fin du ifdef COUPL 
      
enddo !-- fin de l iteration de calcul
#endif
 !---------------------   FIN  MODE NON COUPLE    ------------------

END 
!-------------------------------------------------------------------
!----------------- FIN DU PROGRAMME PRINCIPAL ----------------------
!-------------------------------------------------------------------
#ifdef CALC
!-------------------------------------------------------------------
!
              SUBROUTINE ALLOC_VAR_Eco3M(nbcall)   
!     
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-- subroutine permettant d allouer dynamiquement les concentrations 
!-- et les flux et de les initialiser  en mode couplé

 Use DEF_TYPE
 Use VAR_USER
 Use VAR_GLOBAL
 Use COUPLEUR_PHY_BIO
 
Implicit None

!-- variables locales:
Integer :: i,j,imat,jmat,istat,nn
Integer:: nbcall

#ifdef MODTEST
  write(*,*)'nbcall ds ALLOC_VAR_Eco3M=',nbcall
  write(*,*) 'nx_min,nx_max,ny_min,ny_max,nz_max=',nx_min,nx_max,ny_min,ny_max,nz_max
#endif

!-----------------------------------------------------------------
!  Allocations  de la salinite et de la temperature du modele bio 
!-----------------------------------------------------------------
!-- pour le programme Moohbeedycs, cette allocation est faite dans
!   le programme physique 
!Temperature
  Allocate(TEMP_BIO(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)  
   if (istat/=0) write(*,*) 'pb d allocation de TEMP_BIO'
            
!-- Salinite
 Allocate(SAL_BIO(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat) 
  if (istat/=0) write(*,*) 'pb d allocation de SAL_BIO'

!-- Vent --! Ajout Katixa Lajaunie-Salla 25102018
 Allocate(WIND(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat) 
  if (istat/=0) write(*,*) 'pb d allocation de WIND'

!-- CO2 atm --! Ajout Katixa Lajaunie-Salla 25102018
 Allocate(CO2atm(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat) 
  if (istat/=0) write(*,*) 'pb d allocation de CO2atm'
 Allocate(NH4(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat) 
  if (istat/=0) write(*,*) 'pb d allocation de NH4'
 Allocate(NO3(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat) 
  if (istat/=0) write(*,*) 'pb d allocation de NO3'
 Allocate(PO4(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat) 
  if (istat/=0) write(*,*) 'pb d allocation de PO4'
!-- Lors du 1er appel de cette subroutine, seules la temperature et
!   la salinité sont allouées dynamiquement
!  
!if (nbcall == 1) return  


!-- allocation du champ conc au cas ou on n''utilise pas le fichier declar.inc
!-----------------------------------------------------------------------------
! Allocate(conc_phy(nbvar),STAT=istat)
! if (istat/=0) write(*,*) 'pb d allocation de conc_phy'
! Allocate(conc_phy0(nbvar),STAT=istat)
! if (istat/=0) write(*,*) 'pb d allocation de conc_phy0'
! Allocate(conc_phy_old(nbvar),STAT=istat)
! if (istat/=0) write(*,*) 'pb d allocation de conc_phy_old'

    
! do i=1,nbvar
!   NULLIFY (conc_phy(i)%conc)
!   NULLIFY (conc_phy0(i)%conc)
!   NULLIFY (conc_phy_old(i)%conc)
!   
!   Allocate(conc_phy(i)%conc(nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
!   if (istat/=0) write(*,*) 'pb d allocation de conc_phy%conc'
!   Allocate(conc_phy0(i)%conc(nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
!   if (istat/=0) write(*,*) 'pb d allocation de conc_phy0%conc'
!   Allocate(conc_phy_old(i)%conc(nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
!   if (istat/=0) write(*,*) 'pb d allocation de conc_phy_old%conc'
! enddo

! Allocate(flux_phy(nbvar,nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
! if (istat/=0) write(*,*) 'pb d allocation de flux_phy'

!-----------------------------------------------------------------
!-- Allocations des variables du systeme des carbonates 
!-----------------------------------------------------------------
 Allocate(K0(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de K0' 
#ifdef MODTEST
  write(11,*) 'Alloc de K0'
#endif 
 Allocate(K1(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de K1' 
#ifdef MODTEST
  write(11,*) 'Alloc de K1'
#endif 
 Allocate(K2(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de K2' 
#ifdef MODTEST
  write(11,*) 'Alloc de K2'
#endif  
 Allocate(Ke(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de Ke' 
#ifdef MODTEST
  write(11,*) 'Alloc de Ke'
#endif 
 Allocate(KB(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de KB' 
#ifdef MODTEST
  write(11,*) 'Alloc de KB'
#endif 
 Allocate(Ksp(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de Ksp' 
#ifdef MODTEST
  write(11,*) 'Alloc de Ksp'
#endif 
 Allocate(Ca2(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de Ca2' 
#ifdef MODTEST
  write(11,*) 'Alloc de Ca2'
#endif 
 Allocate(Omega(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de Omega' 
#ifdef MODTEST
  write(11,*) 'Alloc de Omega'
#endif 
 Allocate(AC(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de AC' 
#ifdef MODTEST
  write(11,*) 'Alloc de AC'
#endif 
 Allocate(AC1(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de AC1' 
#ifdef MODTEST
  write(11,*) 'Alloc de AC1'
#endif 
 Allocate(X(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de X' 
#ifdef MODTEST
  write(11,*) 'Alloc de X'
#endif 
 Allocate(Y(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de Y' 
#ifdef MODTEST
  write(11,*) 'Alloc de Y'
#endif 
 Allocate(H(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de H' 
#ifdef MODTEST
  write(11,*) 'Alloc de H'
#endif
 Allocate(OH(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de OH' 
#ifdef MODTEST
  write(11,*) 'Alloc de OH'
#endif
 Allocate(BOH4(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de BOH4' 
#ifdef MODTEST
  write(11,*) 'Alloc de BOH4'
#endif
 Allocate(CO2(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de CO2' 
#ifdef MODTEST
  write(11,*) 'Alloc de CO2'
#endif
 Allocate(CO3(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de CO3' 
#ifdef MODTEST
  write(11,*) 'Alloc de CO3'
#endif
 Allocate(HCO3(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de HCO3' 
#ifdef MODTEST
  write(11,*) 'Alloc de HCO3'
#endif
 Allocate(AT_test(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 Allocate(rho_sw(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 Allocate(TEMP(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de TEMP' 
#ifdef MODTEST
  write(11,*) 'Alloc de TEMP'
#endif
 !-- hauteur d''eau:
 Allocate(prof(nx_min:nx_max,ny_min:ny_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de prof'
#ifdef MODTEST
  write(11,*) 'Alloc de prof'
#endif 
!-- altitude de la maille (par rapport au fond):
 Allocate(alti_z(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de alti_z'
#ifdef MODTEST
  write(11,*) 'Alloc de alti_z'
#endif

!-------------------------------------------------------------------------------
!---               Allocation des variables liées à l''irradiance  
!-------------------------------------------------------------------------------
!---1/  Cas ou l''irradiance varie en fonction de la position horizontale dans l''espace 
!
!-- tableau d''irradiance E(0+) 
! Allocate (irrad(nx_min:nx_max,ny_min:ny_max),STAT=istat)
!  if (istat/=0) write(*,*) 'pb d allocation de irrad'
 
!-- tableau de PAR EPAR(0+) 
! Allocate (E_PAR(nx_min:nx_max,ny_min:ny_max),STAT=istat)
!  if (istat/=0) write(*,*) 'pb d allocation de E_PAR'
 
 
!---2/  Cas ou l''irradiance est la même quelque soient x et y  
!
!-- tableau d''irradiance E(0+) 
 Allocate (irrad(1,1),STAT=istat)
  if (istat/=0) write(*,*) 'pb d allocation de irrad'
#ifdef MODTEST
 write(11,*) 'Alloc de irrad'
#endif 
!-- tableau de PAR EPAR(0+) 
 Allocate (E_PAR(1,1),STAT=istat)
  if (istat/=0) write(*,*) 'pb d allocation de E_PAR'
#ifdef MODTEST
 write(11,*) 'Alloc de E_PAR'
#endif 
 
!-- tableau de PAR selon la verticale
 Allocate (E_PARZ(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
  if (istat/=0) write(*,*) 'pb d allocation de E_PARZ'
#ifdef MODTEST
 write(11,*) 'Alloc de E_PARZ'
#endif 

!-------------------------------------------------------------------------------
!---                Allocation des variables d etat du modele
!-------------------------------------------------------------------------------
 do i=1,nbvar  
  Allocate(VAR(i)%conc(nx_min:nx_max,ny_min:ny_max,nz_max), STAT=istat)
  if (istat/=0) write(*,*) 'pb d allocation de VAR (',i,')'
!  L''initialisation des concentrations provient du modele BIO en mode non couple (subroutine
!  sub_var_ini) et diu modele physique en mode couple. 
!  En attendant, on initialise a zero
   VAR(i)%conc = 0.d0
enddo
#ifdef MODTEST
 write(11,*) 'Alloc de VAR(i)%conc'
#endif 
!-------------------------------------------------------------------------------
!---                Allocation des tableaux de flux
!-------------------------------------------------------------------------------
do i=1,nbvar  
     if (Associated (SELF_VAL(i)%idproc )) then
       Allocate(SELF_VAL(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max), STAT = istat)
       if (istat/=0) write(*,*) 'pb d allocation de FLUX_VAL'
!       SELF_VAL(i)%val = 0.0d0 !--mise a zero des flux biogeochimiques
#ifdef MODTEST
        write(11,*) 'Alloc de SELF_VAL(',i,')%val'
#endif 
     endif
  do j=1, nbvar
     if (Associated (FLUX_VAL(i,j)%idproc )) then
         Allocate(FLUX_VAL(i,j)%val(nx_min:nx_max,ny_min:ny_max,nz_max), STAT=istat)
         if (istat/=0) write(*,*)    'pb d allocation de SELF_VAL(',i,j,')%val'
         FLUX_VAL(i,j)%val = 0.0d0 !--mise a zero des flux biogeochimiques
#ifdef MODTEST
         write(11,*) 'Alloc de FLUX_VAL(',i,j,')%val'
#endif 
     endif
  enddo
enddo

!----------------------------------------------------------------------------------
!-- Allocation et initialisation des termes de tendance :
!----------------------------------------------------------------------------------
Allocate (TEND(nbvar,nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat/=0) write(*,*) 'pb d allocation de TEND'
#ifdef MODTEST
 write(11,*) 'Alloc de TEND'
#endif
TEND = 0.0d0
Allocate (capp(nbvar,nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat/=0) write(*,*) 'pb d allocation de capp'
#ifdef MODTEST
 write(11,*) 'Alloc de capp'
#endif
capp = 0.0d0

!----------------------------------------------------------------------------------
!-- Allocation dyn. de variables globales supplementaires (à completer si besoin):
!----------------------------------------------------------------------------------
!If (CHL_C_BOOL ) then
!   Allocate(CHL_C(nscp_phy))
!     write(11,*) 'Alloc de CHL_C'
!   do i=1,nscp_phy
!!     NULLIFY(Chl_C(i)%val)
!     Allocate(Chl_C(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max), STAT = istat) !--cas CHL:C variable
!     if (istat/=0) write(*,*) 'pb d allocation de CHL_C'
!     write(11,*) 'Alloc de CHL_C(',i,')%val'
!     CHL_C(i)%val(:,:,:)=0.d0
!     Chl_C(i)%idorg=iscp_phy(i)
!   enddo
!Endif

!-- vitesse specifique de production primaire 
if (nscp_phy /=0) Allocate(mu_PPB(nscp_phy),STAT=istat)
  if (istat/=0) write(*,*) 'pb d allocation de mu_PPB, istat=',istat
#ifdef MODTEST
     write(11,*) 'Alloc de mu_PPB'
#endif
!-- vitesse specifique de production primaire (NR = nutrient replete)
if (nscp_phy /=0) Allocate(mu_PPB_NR(nscp_phy),STAT=istat)
  if (istat/=0) write(*,*) 'pb d allocation de mu_PPB_NR'
#ifdef MODTEST
     write(11,*) 'Alloc de mu_PPB_NR'
#endif
!-- vitesse specifique de croissance phytoplanctonique
Allocate(mu_phy(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de mu_phy'
#endif
!!!!!!!!!!!!!!!!!!!!!RAJOUT CHRISTEL 30/10/2007
!-- production primaire phytoplanctonique
Allocate(PPB(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de PPB'
#endif
!-- vitesse d uptake de nh
Allocate(upnh(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de upnh'
#endif
!-- vitesse d uptake de no
Allocate(upno(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de upno'
#endif
!-- fonction quota uptake de nh (ajout Katixa 23052018)
Allocate(hQnh(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de hQnh'
#endif
!-- fonction quota uptake de no (ajout Katixa 23052018)
Allocate(hQno(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de hQno'
#endif
!-- fonction quota uptake de po (ajout Katixa 23052018)
Allocate(hQpo(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de hQpo'
#endif
Allocate(fupnh(nscp_phy),STAT=istat)
Allocate(fupno(nscp_phy),STAT=istat)
Allocate(fuppo(nscp_phy),STAT=istat)
Allocate(Qnh(nscp_phy),STAT=istat)
Allocate(Qpo(nscp_phy),STAT=istat)
!-- respiration du phyto en s-1
Allocate(pr(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de pr'
#endif
!-- production bactérienne en carbone en µmol C l-1 s-1
Allocate(bp(nscp_bac),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de bp'
#endif
! respiration bactérienne en carbone en µmol C l-1 s-1
Allocate(br(nscp_bac),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de br'
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!AJOUT MARION 22/04/10
!-- vitesse d uptake de po
Allocate(uppo(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de uppo'
#endif
!-- limitation phyto N
Allocate(limn(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de limn'
#endif
!-- limitation phyto P
Allocate(limp(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de limp'
#endif
!-- limitation phyto finale
Allocate(lim_nut(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de lim_nut'
#endif
!-- limitation bac C
Allocate(limbaC(nscp_bac),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de limbaC'
#endif
!-- limitation bac N
Allocate(limbaN(nscp_bac),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de limbaN'
#endif
!-- limitation bac P
Allocate(limbaP(nscp_bac),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de limbaP'
#endif
!-- limitation bac finale
Allocate(lim_ba(nscp_bac),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de lim_ba'
#endif
!-- mineralisation
Allocate(remnh(nscp_bac),STAT=istat)
Allocate(rempo(nscp_bac),STAT=istat)
!-- limitation temperature ppb
Allocate(limT_ppb(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de limT_ppb'
#endif
!-- limitation temperature graz 
Allocate(limT_graz(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de limT_graz'
#endif
!-- limitation lumiere
Allocate(lim_lum(nscp_phy),STAT=istat)
#ifdef MODTEST
write(11,*) 'Alloc de lim_lum'
#endif
!!!!!!!!!!!!!!
if (allocated (mu_PPB)) then
  do i=1,nscp_phy
!
     NULLIFY(mu_PPB(i)%val)
#ifdef MODTEST
     write(11,*) 'nullify de mu_PPB(',i,')%val'
#endif
 Allocate(mu_PPB(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
#ifdef MODTEST
     write(11,*) 'Alloc de mu_PPB(',i,')%val'
#endif
     mu_PPB(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de mu_PPB(',i,')%idorg'
#endif
     mu_PPB(i)%val(:,:,:)=0.d0
 enddo
endif
!
if (allocated (mu_PPB_NR)) then
  do i=1,nscp_phy
     NULLIFY(mu_PPB_NR(i)%val)
#ifdef MODTEST
     write(11,*) 'nullify de mu_PPB_NR(',i,')%val'
#endif
     Allocate(mu_PPB_NR(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
#ifdef MODTEST
     write(11,*) 'Alloc de mu_PPB_NR(',i,')%val'
#endif
     mu_PPB_NR(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de mu_PPB_NR(',i,')%idorg'
#endif
     mu_PPB_NR(i)%val(:,:,:)=0.d0
 enddo
endif
!
if (allocated (mu_phy)) then
  do i=1,nscp_phy
     NULLIFY(mu_phy(i)%val)
     Allocate(mu_phy(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
     mu_phy(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de mu_phy(',i,')%idorg'
	write(11,*) 'Nullify de mu_phy(',i,')%idorg'
#endif

     mu_phy(i)%val(:,:,:) = 0.d0
 enddo
endif
!!!!!!!!!!!!!!!!!!!!! RAJOUT CHRISTEL 30/10/2007
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if (allocated (PPB)) then
  do i=1,nscp_phy
     NULLIFY(PPB(i)%val)
     Allocate(PPB(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     PPB(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de PPB(',i,')%idorg'
        write(11,*) 'Nullify de PPB(',i,')%idorg'
#endif

     PPB(i)%val(:,:,:) = 0.d0
 enddo
endif
!
if (allocated (upnh)) then
  do i=1,nscp_phy
     NULLIFY(upnh(i)%val)
     Allocate(upnh(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
     upnh(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de upnh(',i,')%idorg'
	write(11,*) 'Nullify de upnh(',i,')%idorg'
#endif

     upnh(i)%val(:,:,:) = 0.d0
 enddo
endif
!
if (allocated (upno)) then
  do i=1,nscp_phy
     NULLIFY(upno(i)%val)
     Allocate(upno(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
     upno(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de upno(',i,')%idorg'
	write(11,*) 'Nullify de upno(',i,')%idorg'
#endif

     upno(i)%val(:,:,:) = 0.d0
 enddo
endif
!
!!!!!!!!AJOUT MARION 22/04/10
if (allocated (uppo)) then
  do i=1,nscp_phy
     NULLIFY(uppo(i)%val)
     Allocate(uppo(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     uppo(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de uppo(',i,')%idorg'
        write(11,*) 'Nullify de uppo(',i,')%idorg'
#endif

     uppo(i)%val(:,:,:) = 0.d0
 enddo
endif
! -- Ajout Katixa 23/05/2018
if (allocated (hQnh)) then
  do i=1,nscp_phy
     NULLIFY(hQnh(i)%val)
     Allocate(hQnh(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     hQnh(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de hQnh(',i,')%idorg'
        write(11,*) 'Nullify de hQnh(',i,')%idorg'
#endif

     hQnh(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (hQno)) then
  do i=1,nscp_phy
     NULLIFY(hQno(i)%val)
     Allocate(hQno(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     hQno(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de hQno(',i,')%idorg'
        write(11,*) 'Nullify de hQno(',i,')%idorg'
#endif

     hQno(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (hQpo)) then
  do i=1,nscp_phy
     NULLIFY(hQpo(i)%val)
     Allocate(hQpo(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     hQpo(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de hQpo(',i,')%idorg'
        write(11,*) 'Nullify de hQpo(',i,')%idorg'
#endif

     hQpo(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (fupnh)) then
  do i=1,nscp_phy
     NULLIFY(fupnh(i)%val)
     Allocate(fupnh(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     fupnh(i)%idorg=iscp_phy(i)
     fupnh(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (fupno)) then
  do i=1,nscp_phy
     NULLIFY(fupno(i)%val)
     Allocate(fupno(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     fupno(i)%idorg=iscp_phy(i)
     fupno(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (fuppo)) then
  do i=1,nscp_phy
     NULLIFY(fuppo(i)%val)
     Allocate(fuppo(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     fuppo(i)%idorg=iscp_phy(i)
     fuppo(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (Qnh)) then
  do i=1,nscp_phy
     NULLIFY(Qnh(i)%val)
     Allocate(Qnh(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     Qnh(i)%idorg=iscp_phy(i)
     Qnh(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (Qpo)) then
  do i=1,nscp_phy
     NULLIFY(Qpo(i)%val)
     Allocate(Qpo(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     Qpo(i)%idorg=iscp_phy(i)
     Qpo(i)%val(:,:,:) = 0.d0
 enddo
endif
!
if (allocated (limn)) then
  do i=1,nscp_phy
     NULLIFY(limn(i)%val)
     Allocate(limn(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     limn(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de limn(',i,')%idorg'
        write(11,*) 'Nullify de limn(',i,')%idorg'
#endif

     limn(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (limp)) then
  do i=1,nscp_phy
     NULLIFY(limp(i)%val)
     Allocate(limp(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     limp(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de limp(',i,')%idorg'
        write(11,*) 'Nullify de limp(',i,')%idorg'
#endif

     limp(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (lim_nut)) then
  do i=1,nscp_phy
     NULLIFY(lim_nut(i)%val)
     Allocate(lim_nut(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     lim_nut(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de lim_nut(',i,')%idorg'
        write(11,*) 'Nullify de lim_nut(',i,')%idorg'
#endif

     lim_nut(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (limbaC)) then
  do i=1,nscp_bac
     NULLIFY(limbaC(i)%val)
     Allocate(limbaC(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     limbaC(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de limbaC(',i,')%idorg'
        write(11,*) 'Nullify de limbaC(',i,')%idorg'
#endif

     limbaC(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (limbaN)) then
  do i=1,nscp_bac
     NULLIFY(limbaN(i)%val)
     Allocate(limbaN(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     limbaN(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de limbaN(',i,')%idorg'
        write(11,*) 'Nullify de limbaN(',i,')%idorg'
#endif

     limbaN(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (limbaP)) then
  do i=1,nscp_bac
     NULLIFY(limbaP(i)%val)
     Allocate(limbaP(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     limbaP(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de limbaP(',i,')%idorg'
        write(11,*) 'Nullify de limbaP(',i,')%idorg'
#endif

     limbaP(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (lim_ba)) then
  do i=1,nscp_bac
     NULLIFY(lim_ba(i)%val)
     Allocate(lim_ba(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     lim_ba(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de lim_ba(',i,')%idorg'
        write(11,*) 'Nullify de lim_ba(',i,')%idorg'
#endif
     lim_ba(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (remnh)) then
  do i=1,nscp_bac
     NULLIFY(remnh(i)%val)
     Allocate(remnh(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     remnh(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de remnh(',i,')%idorg'
        write(11,*) 'Nullify de remnh(',i,')%idorg'
#endif
     remnh(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (rempo)) then
  do i=1,nscp_bac
     NULLIFY(rempo(i)%val)
     Allocate(rempo(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     rempo(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de rempo(',i,')%idorg'
        write(11,*) 'Nullify de rempo(',i,')%idorg'
#endif
     rempo(i)%val(:,:,:) = 0.d0
 enddo
endif
if (allocated (limT_ppb)) then
  do i=1,nscp_phy
     NULLIFY(limT_ppb(i)%val)
     Allocate(limT_ppb(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     limT_ppb(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de limT_ppb(',i,')%idorg'
        write(11,*) 'Nullify de limT_ppb(',i,')%idorg'
#endif

     limT_ppb(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (limT_graz)) then
  do i=1,nscp_phy
     NULLIFY(limT_graz(i)%val)
     Allocate(limT_graz(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     limT_graz(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de limT_graz(',i,')%idorg'
        write(11,*) 'Nullify de limT_graz(',i,')%idorg'
#endif

     limT_graz(i)%val(:,:,:) = 0.d0
 enddo
endif

if (allocated (lim_lum)) then
  do i=1,nscp_phy
     NULLIFY(lim_lum(i)%val)
     Allocate(lim_lum(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))
     lim_lum(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de lim_lum(',i,')%idorg'
        write(11,*) 'Nullify de lim_lum(',i,')%idorg'
#endif

     lim_lum(i)%val(:,:,:) = 0.d0
 enddo
endif
!!!!!!!!!!!!! fin marion
if (allocated (pr)) then
  do i=1,nscp_phy
     NULLIFY(pr(i)%val)
     Allocate(pr(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
     pr(i)%idorg=iscp_phy(i)
#ifdef MODTEST
     write(11,*) 'Alloc de pr(',i,')%idorg'
	write(11,*) 'Nullify de pr(',i,')%idorg'
#endif

     pr(i)%val(:,:,:) = 0.d0
 enddo
endif
!
if (allocated (bp)) then
  do i=1,nscp_bac
     NULLIFY(bp(i)%val)
     Allocate(bp(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
     bp(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de bp(',i,')%idorg'
	write(11,*) 'Nullify de bp(',i,')%idorg'
#endif

     bp(i)%val(:,:,:) = 0.d0
 enddo
endif
!
if (allocated (br)) then
  do i=1,nscp_bac
     NULLIFY(br(i)%val)
     Allocate(br(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max)) 
     br(i)%idorg=iscp_bac(i)
#ifdef MODTEST
     write(11,*) 'Alloc de br(',i,')%idorg'
	write(11,*) 'Nullify de br(',i,')%idorg'
#endif

     br(i)%val(:,:,:) = 0.d0
 enddo
endif
!-- Init du nombre de bactéries
    allocate(nb_bact(nx_min:nx_max,ny_min:ny_max))

! lecture d un fichier init, genere par matlab, rangé par ligne (171),
! puis colonne (91). Pas de vertical !!
!   ficbac = trim(repconcinit) // 'nb_bact.DAT'
!   inquire(file=ficbac,EXIST=present)
!   if(.NOT.present) then
!     write(*,*) 'fichier nb_bact non trouvee'
!     nb_bact(:,:) = 1.
!   else
!     open(135,FILE=ficbac)
!     do col=ny_min,ny_max
!       Read(135,'(171(2X,e14.7))',iostat=errlec) (nb_bact(ili,col),ili=nx_min,nx_max)
!       where (isnan(nb_bact)) nb_bact=0.
!       if(errlec==-1) then
!         print *, 'fin lecture fichier init nb_bact',ficbac
!         exit
!       endif
!       if (errlec /=0) then
!         print *, 'pb de lecture du fichier de concentration initiale :', ficbac
!         stop
!       endif
!     enddo
!     write(*,*) 'nb_bact init a partir de', ficbac
!   endif
nb_bact(:,:) = 0.2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!--Allocation de la matrice de vitesse de grazing

nn=nscp_phy+nscp_zoo+nscp_bac
if (nscp_zoo /=0) allocate (mu_graz(nscp_zoo,nn),STAT=istat)
if (istat/=0) write(*,*) 'pb d allocation de mu_graz'
#ifdef MODTEST
write(11,*) 'Alloc de mu_graz'
#endif
  
if (nscp_zoo /= 0 .and. allocated(mu_graz)) then
do i=1,nscp_zoo
  do j=1,nn
   NULLIFY(mu_graz(i,j)%val)
#ifdef MODTEST
   write(11,*) 'nullify de mu_graz(',i,j,')%val'
#endif
    Allocate(mu_graz(i,j)%val(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)  !--matrice de vitesse de grazing
#ifdef MODTEST
    write(11,*)'alloc de mu_graz(',i,j,')%val'
#endif
    if (istat/=0) write(*,*) 'pb d allocation de mu_graz%val' 
    mu_graz(i,j)%val(:,:,:)=0.d0
  enddo
 enddo
endif

!--Allocation de la matrice de vitesse de production bact.
!Allocate(mu_BP(nscp_bac))           
!  do i=1,nscp_bac
!      Allocate(mu_BP(i)%val(nx_min:nx_max,ny_min:ny_max,nz_max))           
!      mu_BP(i)%idorg = iscp_bac(i)
!      mu_BP(i)%val = 0.d0
!enddo      
!

#ifdef MODTEST
     write(*,*) 'Fin du ss-programme ALLOC_VAR_Eco3M, nbcall=',nbcall
     write(11,*) 'Fin du ss-programme ALLOC_VAR_Eco3M, nbcall=',nbcall
#endif

End SUBROUTINE ALLOC_VAR_Eco3M
#endif
!-----------------------------------------------------------
!
               SUBROUTINE DEALLOC_VAR_Eco3M(nbcall)   
!
!------------------------------------------------------------
Use DEF_TYPE
Use VAR_USER
Use COUPLEUR_PHY_BIO
Use VAR_GLOBAL
Implicit None

!-- arguments:
integer:: nbcall

! Variables locales 
 Integer::i,j,k,istat,nn
 LOGICAL :: test

!--initialisation
!istat = 0

#ifdef MODTEST
write(*,*)'Dans DEALLOC_VAR_Eco3M, nbcall =',nbcall
#endif

!
if (nbcall==2) then 
  deallocate (tab_noms_fich, STAT=istat)
  if (istat /=0) write(*,*) 'pb de desallocation de tab_noms_fich'
#ifdef MODTEST
  write(11,*) 'dealloc de tab_noms'
#endif
endif

 do  i=1,nbvar 
!-------------------------------------------------------------------------------
!---               Desallocation des variables d etat du modele
!-------------------------------------------------------------------------------
       if(associated (VAR(i)%conc))  deallocate(VAR(i)%conc)
!       NULLIFY(VAR(i)%conc)
#ifdef MODTEST
       write(11,*)'dealloc de VAR(',i,')%conc'
#endif
!-------------------------------------------------------------------------------
!---               Desallocation des tableaux de flux
!-------------------------------------------------------------------------------
!      if(associated (SELF_VAL(i)%idproc))    NULLIFY(SELF_VAL(i)%val)
       if (associated (SELF_VAL(i)%val)) deallocate(SELF_VAL(i)%val)
#ifdef MODTEST
       write(11,*)'dealloc de SELF_VAL(',i,')%val'
#endif
       do j=1,nbvar 
         if(associated (FLUX_VAL(i,j)%idproc))  then 
           deallocate (FLUX_VAL(i,j)%val)
#ifdef MODTEST
           write(11,*)'dealloc de FLUX_VAL(',i,j,')%val'
#endif
!           deallocate(FLUX_VAL(i,j)%val, STAT=istat)
!           if (istat /= 0)  stop 'PB Dans DEALLOC_VAR_Eco3M'
         endif
       enddo
 enddo

!----------------------------------------------------------------------------------
!-- Desallocation de variables globales supplementaires (à completer si besoin):
!----------------------------------------------------------------------------------
!
!-- deallocation des variables liees a la production primaire 
if (nscp_phy /= 0) then
 
 
! if (allocated (mu_phy)) then 
!  do i=1,nscp_phy
!       if(associated(mu_phy(i)%val))  deallocate (mu_phy(i)%val,STAT=istat)
!       if (istat /=0) write(*,*) 'pb de desallocation de mu_phy%val'
!  enddo
!  Deallocate(mu_phy,STAT=istat) !--matrice de vitesse de prod. primaire
!  if (istat /=0) write(*,*) 'pb de desallocation de mu_phy'
! endif  

 if (allocated (mu_PPB)) then
  do i=1,nscp_phy
      if(associated(mu_PPB(i)%val))  then
        deallocate (mu_PPB(i)%val)
#ifdef MODTEST
        write(11,*)'dealloc de mu_PPB(',i,')%val'
#endif
     endif
  enddo
  Deallocate(mu_PPB,STAT=istat) !--vitesse specifique de production primaire
  if (istat /=0) write(*,*) 'pb de desallocation de mu_PPB'
#ifdef MODTEST
  write(11,*)'desallocation de  mu_PPB'
#endif
 endif  

 if (allocated (mu_PPB_NR)) then
  do i=1,nscp_phy
      if(associated(mu_PPB_NR(i)%val))  then
        deallocate (mu_PPB_NR(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de mu_PPB_NR%val'
#ifdef MODTEST
        write(11,*)'dealloc de mu_PPB_NR(',i,')%val'
#endif
      endif
  enddo
  Deallocate(mu_PPB_NR,STAT=istat) !--vitesse specifique de prod. primaire NR
  if (istat /=0) write(*,*) 'pb de desallocation de mu_PPB_NR'
#ifdef MODTEST
  write(11,*)'desallocation de  mu_PPB_NR'
#endif
 endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!RAJOUT CHRISTEL 30/10/2007
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (allocated (mu_phy)) then
  do i=1,nscp_phy
      if(associated(mu_phy(i)%val))  then
        deallocate (mu_phy(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de mu_phy%val'
#ifdef MODTEST
        write(11,*)'dealloc de mu_phy(',i,')%val'
#endif
      endif
  enddo
  Deallocate(mu_phy,STAT=istat) !--vitesse specifique de croissance du phyto
  if (istat /=0) write(*,*) 'pb de desallocation de mu_phy'
#ifdef MODTEST
  write(11,*)'desallocation de  mu_phy'
#endif
 endif  
if (allocated (PPB)) then
  do i=1,nscp_phy
      if(associated(PPB(i)%val))  then
        deallocate (PPB(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de PPB%val'
#ifdef MODTEST
        write(11,*)'dealloc de PPB(',i,')%val'
#endif
      endif
  enddo
  Deallocate(PPB,STAT=istat) !--production primaire du phyto
  if (istat /=0) write(*,*) 'pb de desallocation de PPB'
#ifdef MODTEST
  write(11,*)'desallocation de  PPB'
#endif
 endif
if (allocated (upnh)) then
  do i=1,nscp_phy
      if(associated(upnh(i)%val))  then
        deallocate (upnh(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de upnh%val'
#ifdef MODTEST
        write(11,*)'dealloc de upnh(',i,')%val'
#endif
      endif
  enddo
  Deallocate(upnh,STAT=istat) !--vitesse d uptake nh du phyto
  if (istat /=0) write(*,*) 'pb de desallocation de nh'
#ifdef MODTEST
  write(11,*)'desallocation de  upnh'
#endif
 endif  
if (allocated (upno)) then
  do i=1,nscp_phy
      if(associated(upno(i)%val))  then
        deallocate (upno(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de upno%val'
#ifdef MODTEST
        write(11,*)'dealloc de upno(',i,')%val'
#endif
      endif
  enddo
  Deallocate(upno,STAT=istat) !--vitesse d uptake no du phyto
  if (istat /=0) write(*,*) 'pb de desallocation de no'
 endif 
if (allocated (fupnh)) then
  do i=1,nscp_phy
      if(associated(fupnh(i)%val))  then
        deallocate (fupnh(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de fupnh%val'
      endif
  enddo
endif
if (allocated (fupno)) then
  do i=1,nscp_phy
      if(associated(fupno(i)%val))  then
        deallocate (fupno(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de fupno%val'
      endif
  enddo
endif
if (allocated (fuppo)) then
  do i=1,nscp_phy
      if(associated(fuppo(i)%val))  then
        deallocate (fuppo(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de fuppo%val'
      endif
  enddo
endif
if (allocated (Qnh)) then
  do i=1,nscp_phy
      if(associated(Qnh(i)%val))  then
        deallocate (Qnh(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de Qnh%val'
      endif
  enddo
endif
if (allocated (Qpo)) then
  do i=1,nscp_phy
      if(associated(Qpo(i)%val))  then
        deallocate (Qpo(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de Qpo%val'
      endif
  enddo
endif
!!!!!!!!!!AJOUT MARION 22/04/10
if (allocated (uppo)) then
  do i=1,nscp_phy
      if(associated(uppo(i)%val))  then
        deallocate (uppo(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de uppo%val'
#ifdef MODTEST
        write(11,*)'dealloc de uppo(',i,')%val'
#endif
      endif
  enddo
  Deallocate(uppo,STAT=istat) !--vitesse d uptake po du phyto
  if (istat /=0) write(*,*) 'pb de desallocation de po'
#ifdef MODTEST
  write(11,*)'desallocation de  uppo'
#endif
 endif
! ajout Katixa 01062018
!--fonction quota d uptake no du phyto
if (allocated (hQno)) then
  do i=1,nscp_phy
      if(associated(hQno(i)%val))  then
        deallocate (hQno(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de hQno%val'
#ifdef MODTEST
        write(11,*)'dealloc de hQno(',i,')%val'
#endif
      endif
  enddo
  Deallocate(hQno,STAT=istat)
  if (istat /=0) write(*,*) 'pb de desallocation de hQno'
#ifdef MODTEST
  write(11,*)'desallocation de  hQno'
#endif
 endif

!--fonction quota d uptake nh du phyto
if (allocated (hQnh)) then
  do i=1,nscp_phy
      if(associated(hQnh(i)%val))  then
        deallocate (hQnh(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de hQnh%val'
#ifdef MODTEST
        write(11,*)'dealloc de hQnh(',i,')%val'
#endif
      endif
  enddo
  Deallocate(hQnh,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de hQnh'
#ifdef MODTEST
  write(11,*)'desallocation de  hQnh'
#endif

 endif
!--fonction quota uptake po du phyto
if (allocated (hQpo)) then
  do i=1,nscp_phy
      if(associated(hQpo(i)%val))  then
        deallocate (hQpo(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de hQpo%val'
#ifdef MODTEST
        write(11,*)'dealloc de hQpo(',i,')%val'
#endif
      endif
  enddo
  Deallocate(hQpo,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de hQpo'
#ifdef MODTEST
  write(11,*)'desallocation de  hQpo'
#endif
 endif
!--limitation N d uptake  du phyto 
 if (allocated (limn)) then
  do i=1,nscp_phy
      if(associated(limn(i)%val))  then
        deallocate (limn(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de limn%val'
#ifdef MODTEST
        write(11,*)'dealloc de limn(',i,')%val'
#endif
      endif
  enddo
  Deallocate(limn,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de limn'
#ifdef MODTEST
  write(11,*)'desallocation de limn'
#endif
 endif

!--limitation P d uptake  du phyto 
 if (allocated (limp)) then
  do i=1,nscp_phy
      if(associated(limp(i)%val))  then
        deallocate (limp(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de limp%val'
#ifdef MODTEST
        write(11,*)'dealloc de limp(',i,')%val'
#endif
      endif
  enddo
  Deallocate(limp,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de limp'
#ifdef MODTEST
  write(11,*)'desallocation de limp'
#endif
 endif

!--limitation finale d uptake  du phyto 
  if (allocated (lim_nut)) then
  do i=1,nscp_phy
      if(associated(lim_nut(i)%val))  then
        deallocate (lim_nut(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de lim_nut%val'
#ifdef MODTEST
        write(11,*)'dealloc de lim_nut(',i,')%val'
#endif
      endif
  enddo
  Deallocate(lim_nut,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de lim_nut'
#ifdef MODTEST
  write(11,*)'desallocation de lim_nut'
#endif
 endif

!--limitation C d uptake  des bacteries 
  if (allocated (limbaC)) then
  do i=1,nscp_bac
      if(associated(limbaC(i)%val))  then
        deallocate (limbaC(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de limbaC%val'
#ifdef MODTEST
        write(11,*)'dealloc de limbaC(',i,')%val'
#endif
      endif
  enddo
  Deallocate(limbaC,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de limbaC'
#ifdef MODTEST
  write(11,*)'desallocation de limbaC'
#endif
 endif

!--limitation P d uptake  des bacteries
   if (allocated (limbaP)) then
  do i=1,nscp_bac
      if(associated(limbaP(i)%val))  then
        deallocate (limbaP(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de limbaP%val'
#ifdef MODTEST
        write(11,*)'dealloc de limbaP(',i,')%val'
#endif
      endif
  enddo
  Deallocate(limbaP,STAT=istat)
  if (istat /=0) write(*,*) 'pb de desallocation de limbaP'
#ifdef MODTEST
  write(11,*)'desallocation de limbaP'
#endif
 endif

!--limitation  d uptake  des bacteries 
    if (allocated (lim_ba)) then
  do i=1,nscp_bac
      if(associated(lim_ba(i)%val))  then
        deallocate (lim_ba(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de lim_ba%val'
#ifdef MODTEST
        write(11,*)'dealloc de lim_ba(',i,')%val'
#endif
      endif
  enddo
  Deallocate(lim_ba,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de lim_ba'
#ifdef MODTEST
  write(11,*)'desallocation de lim_ba'
#endif
 endif
!--mineralisation nh4
    if (allocated (remnh)) then
  do i=1,nscp_bac
      if(associated(remnh(i)%val))  then
        deallocate (remnh(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de remnh%val'
#ifdef MODTEST
        write(11,*)'dealloc de remnh(',i,')%val'
#endif
      endif
  enddo
  Deallocate(remnh,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de remnh'
#ifdef MODTEST
  write(11,*)'desallocation de remnh'
#endif
 endif
!--mineralisation po4
    if (allocated (rempo)) then
  do i=1,nscp_bac
      if(associated(rempo(i)%val))  then
        deallocate (rempo(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de rempo%val'
#ifdef MODTEST
        write(11,*)'dealloc de rempo(',i,')%val'
#endif
      endif
  enddo
  Deallocate(rempo,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de rempo'
#ifdef MODTEST
  write(11,*)'desallocation de rempo'
#endif
 endif
!--limitation  temperature ppb
    if (allocated (limT_ppb)) then
  do i=1,nscp_phy
      if(associated(limT_ppb(i)%val))  then
        deallocate (limT_ppb(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de limT_ppb%val'
#ifdef MODTEST
        write(11,*)'dealloc de limT_ppb(',i,')%val'
#endif
      endif
  enddo
  Deallocate(limT_ppb,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de limT_ppb'
#ifdef MODTEST
  write(11,*)'desallocation de limT_ppb'
#endif
 endif
!--limitation  temperature grazing
    if (allocated (limT_graz)) then
  do i=1,nscp_phy
      if(associated(limT_graz(i)%val))  then
        deallocate (limT_graz(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de limT_graz%val'
#ifdef MODTEST
        write(11,*)'dealloc de limT_graz(',i,')%val'
#endif
      endif
  enddo
  Deallocate(limT_graz,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de limT_graz'
#ifdef MODTEST
  write(11,*)'desallocation de limT_graz'
#endif
 endif
!--limitation lumiere
    if (allocated (lim_lum)) then
  do i=1,nscp_phy
      if(associated(lim_lum(i)%val))  then
        deallocate (lim_lum(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de lim_lum%val'
#ifdef MODTEST
        write(11,*)'dealloc de lim_lum(',i,')%val'
#endif
      endif
  enddo
  Deallocate(lim_lum,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de lim_lum'
#ifdef MODTEST
  write(11,*)'desallocation de lim_lum'
#endif
 endif
!--vitesse de respiration du phyto
if (allocated (pr)) then
  do i=1,nscp_phy
      if(associated(pr(i)%val))  then
        deallocate (pr(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de pr%val'
#ifdef MODTEST
        write(11,*)'dealloc de pr(',i,')%val'
#endif
      endif
  enddo
  Deallocate(pr,STAT=istat) 
  if (istat /=0) write(*,*) 'pb de desallocation de pr'
#ifdef MODTEST
  write(11,*)'desallocation de  pr'
#endif
 endif  
 !--vitesse de production bacterienne
if (allocated (bp)) then
  do i=1,nscp_bac
      if(associated(bp(i)%val))  then
        deallocate (bp(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de bp%val'
#ifdef MODTEST
        write(11,*)'dealloc de bp(',i,')%val'
#endif
      endif
  enddo
  Deallocate(bp,STAT=istat)
  if (istat /=0) write(*,*) 'pb de desallocation de bp'
#ifdef MODTEST
  write(11,*)'desallocation de  bp'
#endif
 endif  
!--vitesse de respiration des bacteries
if (allocated (br)) then
  do i=1,nscp_bac
      if(associated(br(i)%val))  then
        deallocate (br(i)%val,STAT=istat)
        if (istat /=0) write(*,*) 'pb de desallocation de br%val'
#ifdef MODTEST
        write(11,*)'dealloc de br(',i,')%val'
#endif
      endif
  enddo
  Deallocate(br,STAT=istat)
  if (istat /=0) write(*,*) 'pb de desallocation de br'
#ifdef MODTEST
  write(11,*)'desallocation de  br'
#endif
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif
!-- desallocation de la vitesse specifique de grazing
 if (nscp_zoo /= 0 .and. allocated(mu_graz)) then
  nn=nscp_phy+nscp_zoo+nscp_bac
  do i=1,nscp_zoo
   do j=1,nn
       if(associated(mu_graz(i,j)%val))  then
         deallocate(mu_graz(i,j)%val)
#ifdef MODTEST
         write(11,*)'dealloc de mu_graz(',i,')%val'
#endif
       endif
   enddo
  enddo
 
  Deallocate(mu_graz,STAT=istat)  !--matrice de vitesse de grazing
  if (istat /=0) write(*,*) 'pb de desallocation de mu_graz'
#ifdef MODTEST
  write(11,*)'dealloc mu_graz'
#endif
endif  

!if (allocated(bp)) then
!  Deallocate (bp,STAT=istat) !--matrice de vitesse de prod. bact.
!  if (istat /=0) write(*,*) 'pb de desallocation de bp'  
!endif
!
if (CHL_C_BOOL ) then
     if (nscp_phy /= 0 .and. allocated (CHL_C)) then  
        do i=1,nscp_phy
          if(associated (CHL_C(i)%val))  deallocate (CHL_C(i)%val)
#ifdef MODTEST
          write(11,*) 'dealloc  CHL_C(',i,')%val'
#endif
        enddo
     Deallocate (CHL_C)
#ifdef MODTEST
     write(11,*) 'dealloc CHL_C'
#endif
     endif
 endif
!----------------------------------------------------------------------------------
!-- Desallocation des termes de tendance :
!----------------------------------------------------------------------------------
if (allocated(TEND)) then
   deallocate (TEND,STAT=istat)
   if (istat /=0) write(*,*) 'pb de desallocation de TEND'
#ifdef MODTEST
   write(11,*) 'dealloc  TEND'
#endif
endif
if (allocated(capp)) then
   deallocate (capp,STAT=istat)
   if (istat /=0) write(*,*) 'pb de desallocation de capp'
#ifdef MODTEST
   write(11,*) 'dealloc  capp'
#endif
endif
 !-- hauteur d''eau:
if (allocated (prof)) then
   deallocate(prof,STAT=istat)
   if (istat/=0) write(*,*) 'pb de desallocation de prof'
#ifdef MODTEST
   write(11,*) 'dealloc  prof'
#endif
endif
 
!-- altitude de la maille (par rapport au fond):
if (allocated (alti_z)) then
    deallocate(alti_z,STAT=istat)
    if (istat/=0) write(*,*) 'pb de desallocation de alti_z'
#ifdef MODTEST
    write(11,*) 'dealloc  alti_z'
#endif
endif

 
!-------------------------------------------------------------------------------
!---         Desallocation des variables liées à l''irradiance  
!-------------------------------------------------------------------------------
!-- tableau d''irradiance E(0+) desalloue dans main_BIO
if (allocated (irrad)) then
    deallocate (irrad,STAT=istat)
    if (istat/=0) write(*,*) 'pb de desallocation de irrad'
#ifdef MODTEST
    write(11,*) 'dealloc irrad'
#endif
endif
  
!-- tableau de PAR EPAR(0+) desalloue dans main_BIO
if (allocated (E_PAR)) then
   deallocate (E_PAR,STAT=istat)
    if (istat/=0) write(*,*) 'pb de desallocation de E_PAR'
#ifdef MODTEST
    write(11,*) 'dealloc E_PAR'
#endif
endif
  
!-- tableau de PAR selon la verticale
if (allocated (E_PARZ)) then
   deallocate (E_PARZ,STAT=istat)
   if (istat/=0) write(*,*) 'pb de desallocation de E_PARZ'
#ifdef MODTEST
   write(11,*) 'dealloc E_PARZ'
#endif
endif
  
!-- Temperature
if (allocated (TEMP_BIO)) then
   deallocate(TEMP_BIO,STAT=istat)
   if (istat /=0) write(*,*) 'pb de desallocation de TEMP_BIO'
#ifdef MODTEST
   write(11,*) 'dealloc TEMP_BIO'
#endif
endif
             
!-- Salinite
if (allocated (SAL_BIO)) then
   deallocate(SAL_BIO, STAT=istat)
   if (istat /=0) write(*,*) 'pb de desallocation de SAL_BIO'
#ifdef MODTEST
   write(11,*) 'dealloc SAL_BIO'
#endif
endif

!-- Vent
if (allocated (WIND)) then
   deallocate(WIND, STAT=istat)
   if (istat /=0) write(*,*) 'pb de desallocation de WIND'
#ifdef MODTEST
   write(11,*) 'dealloc WIND'
#endif
endif

!-- CO2 atm
if (allocated (CO2atm)) then
   deallocate(CO2atm, STAT=istat)
endif
! densite rho_sw
if (allocated (rho_sw)) then
   deallocate(rho_sw, STAT=istat)
endif
!-- NH4 atm
if (allocated (NH4)) then 
   deallocate(NH4, STAT=istat)
endif
!-- NO3 atm
if (allocated (NO3)) then
   deallocate(NO3, STAT=istat)
endif
!-- PO4 atm
if (allocated (PO4)) then
   deallocate(PO4, STAT=istat)
endif
End SUBROUTINE DEALLOC_VAR_Eco3M 
!-----------------------------------------------------------

#ifdef COUPL
!---------------------------------------------------------
            Subroutine Recup_dim_maill_Mobee

! Cette subroutine est specifique au couplage avec le code physique Mobeehdycs
! mais peut etre adaptéé a chaque code physique. L objectif en est de recuperer la
! dimension du maillage, en allant la lire directement dans le fichier de donnees 
! du code physique
!------------------------------------------------------------
USE COUPLEUR_PHY_BIO
 Implicit none
 Integer:: lucon,toto
 real:: toto_r
 character(100)::fcon
 character(50)::fgeo,fich

 lucon=120
! fich = trim(adjustl(rep_physique))//'Data_Old/'//'penteb.con'
! fich = 'penteb.con'
fich = 'senlagune.con'
 open(lucon,File=fich)
 write(*,*) 'lecture du nb de couches ds le fichier ',fich
      read(lucon,*)fcon   
      read(lucon,*)fcon   
      read(lucon,*)fcon   
      read(lucon,*)fcon   
      read(lucon,*)fcon   
      read(lucon,*)fcon   
      read(lucon,*)fcon   
      read(lucon,*)toto,toto, nz_max
      write(*,*) 'nb de couches',nz_max
 close(lucon)
 
!fich = trim(adjustl(rep_physique))//'Data_Old/'//'penteb.geo'
!fich = 'penteb.geo'
fich = 'senlagune.geo'

open(lucon,File=fich)
write(*,*) 'Recuperation du maillage horizontal ds le fichier',fich
      read(lucon,*)fgeo   
      read(lucon,*)fgeo   
      read(lucon,*)fgeo   
      read(lucon,*)fgeo   
      read(lucon,*)ny_max,nx_max,toto_r
      write(*,*) 'taille du maillage en x et y',nx_max,ny_max
close(lucon)

! Specifique a Mobeehdycs :
    nx_min=1 !maille initiale suivant  x (N/S)
    ny_min=1 !maille initiale suivant y (E/W)

#ifdef INI
 open(100,File='declar.inc')
  write(100,*)'! Ce fichier est créé automatiquement par le code lors de son exécution en'
  write(100,*)'! mode couplé' 
  write(100,*) '' 
  write(100,*)'!-- Definition d un nouveau type (specifique au couplage avec la physique)'
  write(100,*)'TYPE VAR_COUPL'
  write(100,*)'Integer       :: idorg   ! identifiant unique pour chaque organisme'
  write(100,*)'Integer       :: ielmt    ! ielmt = 1 pour la conc de ref sur laquelle on'
  write(100,*)'                         ! calcule le transport '
  write(100,*)'! matrice contenant les concentrations ds l espace :'
  write(100,'(A16,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1)')'Real(8):: conc(',nz_max,',',ny_min,':',ny_max,',',nx_min,':',nx_max,')'
  write(100,*)'END TYPE VAR_COUPL'
  close(100)
  write(*,*) 'declar.inc cree'
#endif

end Subroutine Recup_dim_maill_Mobee
!---------------------------------------------------------
#ifdef CALC
!---------------------------------------------------------
      SUBROUTINE TRANSFERT_Eco3M_PHY
!---------------------------------------------------------
!-- subroutine permettant de transferer les variables 
!   biogeochimiques au code physique Mobeehdycs sous une forme 
!  adéquate pour différencier les termes sur lesquels
!  on doit calculer le transport hydrodynamique

Use COUPLEUR_PHY_BIO
Implicit None

!-- variables locales:
integer :: ii,num_elmt
integer :: l,n,m
 
!-- Initialisation
do ii = 1,nbvar
      conc_phy(ii)%conc  = 0.d0
      conc_phy(ii)%idorg = 0
      conc_phy(ii)%ielmt = 0
      flux_phy(ii,:,:,:) =0.d0
enddo
!
num_elmt = 1
do ii=1,nbvar
   do l=1,nz_max
    do n=ny_min,ny_max
     do m=nx_min,nx_max
       flux_phy(ii,l,n,m)= TEND(ii,m,n,l)
       conc_phy(ii)%conc(l,n,m)=VAR(ii)%conc(m,n,l)
     enddo
    enddo
   enddo
   conc_phy(ii)%idorg= VAR(ii)%idorg

!  Affectation du numero d element (1 pour la variable transportee
!  et autre numero pour les variables associees)
   conc_phy(ii)%ielmt= num_elmt
!-- incrementation de l indice correspondant aux elements
!   dans l ordre ou ils sont donnes dans config.ini
!  Autrement dit, si le phyto est exprime en 1er en carbone
!  dans config.ini, le transport sera calcule par rapport a la conc.
!  en carbone
   num_elmt = num_elmt +1
   if ( num_elmt > nb_elmt(VAR(ii)%idorg)) num_elmt=1 ! on remet a 1 pour un
                                                         ! nouvel organisme
   enddo

End SUBROUTINE TRANSFERT_Eco3M_PHY

!---------------------------------------------------------
      SUBROUTINE TRANSFERT_PHY_Eco3M
!---------------------------------------------------------
!-- subroutine permettant de transferer les variables 
!   biogeochimiques du code Mobeehdycs vers Eco3M
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Implicit None


!-- fonctions externes
Real(8) :: altitude_mil
Real(8) :: profondeur

!-- variables locales:
integer :: ii, l,n,m


do ii=1,nbvar
   do l=1,nz_max
    do n=ny_min,ny_max
     do m=nx_min,nx_max
      VAR(ii)%conc(m,n,l)=conc_phy(ii)%conc(l,n,m)
!-- specifique a Moohbeedycs: le transfert de Salinite et
!  temperature se fait dans le code physique
!      TEMP_BIO(m,n,l) = TEMP_PHY(l,n,m)
!      SAL_BIO(m,n,l)  = TEMP_PHY(l,n,m)
     enddo
    enddo
   enddo
enddo


 !-- profondeur (hauteur d''eau):
  do m = ny_min, ny_max
    do n = nx_min, nx_max
         prof(n,m) = profondeur(m,n)
    enddo
  enddo
 
 !-- Immersions des mailles :
  do l = 1,nz_max
    do m = ny_min, ny_max
      do n = nx_min, nx_max
        alti_z(n,m,l) = altitude_mil(l,m,n)
      enddo
    enddo
  enddo

End SUBROUTINE TRANSFERT_PHY_Eco3M

!-----------------------------------------------------------
     SUBROUTINE ALLOC_VAR_PHY 
!-----------------------------------------------------------
!-----------------------------------------------------------
!-- subroutine permettant d allouer dynamiquement les concentrations et les flux
!-- qui vont etre utilisés ds le code physique

! ATTENTION : 
! ds le code bio, les indices sont ds l ordre suivant les axes  x,y,z
! ds le code phy, les indices sont ds l ordre suivant les axes  z,y,x

Use COUPLEUR_PHY_BIO
Implicit None

!-- variables locales:
Integer :: i,istat

!-- Allocation de conc_phy  et flux_phy :
if (.not. allocated (conc_phy)) then
 Allocate(conc_phy(nbvar),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de conc_phy'
#ifdef MODTEST
   write(11,*)
   write(11,*)'allocation de conc_phy a ',nbvar
#endif
endif
if (.not. allocated (conc_phy0)) then
 Allocate(conc_phy0(nbvar),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de conc_phy0'
#ifdef MODTEST
  write(11,*)'allocation de conc_phy0 a ',nbvar
#endif
endif
if (.not. allocated (conc_phy_old)) then
 Allocate(conc_phy_old(nbvar),STAT=istat)
 if (istat/=0) write(*,*) 'pb d allocation de conc_phy_old'
#ifdef MODTEST
  write(11,*)'allocation de conc_phy_old a ',nbvar
#endif
endif

    
!-- allocation du champ conc au cas ou on n''utilise pas le fichier declar.inc
!-----------------------------------------------------------------------------
! do i=1,nbvar
!   NULLIFY (conc_phy(i)%conc)
!   NULLIFY (conc_phy0(i)%conc)
!   NULLIFY (conc_phy_old(i)%conc)
!   
!   Allocate(conc_phy(i)%conc(nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
!   if (istat/=0) write(*,*) 'pb d allocation de conc_phy%conc'
!   Allocate(conc_phy0(i)%conc(nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
!   if (istat/=0) write(*,*) 'pb d allocation de conc_phy0%conc'
!   Allocate(conc_phy_old(i)%conc(nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
!   if (istat/=0) write(*,*) 'pb d allocation de conc_phy_old%conc'
! enddo

if (.not. allocated (flux_phy)) then
  Allocate(flux_phy(nbvar,nz_max,ny_min:ny_max,nx_min:nx_max),STAT=istat)
  if (istat/=0) write(*,*) 'pb d allocation de flux_phy'
#ifdef MODTEST
  write(11,*)'allocation de flux_phy '
#endif
endif 


 End SUBROUTINE ALLOC_VAR_PHY   
!-----------------------------------------------------------
!----------------------------------------------------------
!
     SUBROUTINE DEALLOC_VAR_PHY   
!     
!-----------------------------------------------------------
!-- subroutine permettant de desallouer dynamiquement les concentrations et les flux
!-- qui vont etre utilisés ds le code physique

Use COUPLEUR_PHY_BIO
Implicit None

!-- variables locales:
Integer :: i,istat

!-- Desallocation de conc_phy  et flux_phy :
if (allocated (conc_phy))  then
      Deallocate(conc_phy,STAT=istat)
      if (istat/=0) write(*,*) 'pb de desallocation de conc_phy'
#ifdef MODTEST
      write(11,*)'desallocation de conc_phy '
#endif
else
      write(*,*) 'pb de desallocation de conc_phy, conc_phy non alloué'
endif
 
!-- Desallocation de conc_phy  et flux_phy :
if (allocated (conc_phy0)) then
     Deallocate(conc_phy0,STAT=istat)
     if (istat/=0) write(*,*) 'pb de desallocation de conc_phy0'
#ifdef MODTEST
     write(11,*)'desallocation de conc_phy0'
#endif
endif
 
!-- Desallocation de conc_phy  et flux_phy :
if (allocated (conc_phy_old)) then
     Deallocate(conc_phy_old,STAT=istat)
     if (istat/=0) write(*,*) 'pb de desallocation de conc_phy_old'
#ifdef MODTEST
     write(11,*)'desallocation de conc_phy_old'
#endif
endif
 
if (allocated (flux_phy)) then
     Deallocate(flux_phy,STAT=istat)
     if (istat/=0) write(*,*) 'pb de desallocation de flux_phy'
#ifdef MODTEST
     write(11,*)'desallocation de flux_phy'
#endif
endif
 
 End SUBROUTINE DEALLOC_VAR_PHY
!----------------------------------------------------------
#endif
#endif
