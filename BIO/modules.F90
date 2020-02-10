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
!
                            MODULE DEF_TYPE
!			
!-- Ce module contient les definitions des types derives, 
!   et de certaines variables globales spécifiques au modèle
!    biogéochimique (e.g. FLUX_VAL(:,:) ; SELF_VAL(:) ;...)
!
! Dernière modification: 27/06/07
!
!-------------------------------------------------------
USE LONG_CHAIN

Implicit None

!-- Definition de nouveaux types:

!-- type pour variables d etat 
TYPE VAR_ETAT
 Integer       :: idorg   ! identifiant unique pour chaque organisme
                          ! (permet de relier simplement les differentes
                          ! conc. d une meme espece)
 Character(3)  :: comp    ! nom du compartiment (exemple: phy)
 Character(10) :: scomp   ! nom du sous-compartiment (exemple: diatom)
 Character (5) :: elmt    ! exprime en mol de C ou de N, de P, de Si...
 Real(8),dimension(:,:,:),pointer:: conc => NULL() !matrice contenant les concentrations ds l espace
END TYPE VAR_ETAT

!-- type pour les parametres associes a chaque processus:
TYPE PARAM
 Integer      :: ipos(2)   ! couple d indices des variables concernees par le flux
 Integer      :: idproc    ! identifiant associe a un type de processus
 Character(L_VAR)      :: signe  ! + ou - ou *
 Real(8),dimension(:),pointer:: valpar => NULL()  ! valeur des paramètres associes au processus
                                        ! mis en jeu dans le flux de ipos(1) --> ipos(2)
END TYPE  PARAM

!-- Type pour les flux entre 2 variables:
TYPE FLUX
  Integer, dimension(:), pointer  :: idproc => NULL()
! nb de flux composant le flux global entre les variables   
  Integer                         :: nb_sflux
  Real(8),dimension(:,:,:),pointer:: val => NULL()
END TYPE FLUX

!-- Type pour les processus:
TYPE PROC
 Integer           :: idproc
 Character (14)    :: nomproc
 Character (L_VAR) :: nomsub 
 Integer(2)        :: nbpar
 Character(10),dimension(:),pointer::nompar => NULL()
END TYPE PROC

!-- Type pour les variables globales specifiques (definies dans le module VAR_GLOBAL)
 TYPE VAR_GLOB_USER
 Integer       :: idorg   ! identifiant unique pour chaque organisme
                          ! (permet de relier simplement les differentes
                          ! conc. d une meme espece)
 Real(8),dimension(:,:,:),pointer:: val => NULL() !matrice contenant les valeurs de la variable ds l espace
END TYPE VAR_GLOB_USER

!-- Type pour les variables globales specifiques (definies dans le module VAR_GLOBAL)
TYPE FLUX_GLOB_USER
 Real(8),dimension(:,:,:),pointer:: val => NULL()
END TYPE FLUX_GLOB_USER

!--- Declaration des tableaux:
!-- Tableau des modeles de processus pris en compte par le code compilé
 TYPE(PROC)      ,Allocatable    :: PROC_MOD(:)

!-- Matrice des flux entre variables d etat:
 TYPE(FLUX)     ,Allocatable     :: FLUX_VAL(:,:)

!-- Vecteur des flux n impliquant qu une seule des 
!   variables d etat
 TYPE(FLUX)     ,Allocatable     :: SELF_VAL(:)

!-- Matrice et vecteur des parametres associes aux flux
 TYPE(PARAM)     ,Allocatable    :: FLUX_PAR(:)
 TYPE(PARAM)     ,Allocatable    :: SELF_PAR(:)

END MODULE DEF_TYPE
!
!----------------------------------------------------------
!
                      Module VAR_USER
!
!-- Ce module contient des parametres de l utilisateur		
!----------------------------------------------------------
USE LONG_CHAIN

Implicit None

 !-- Paramètre defini par l utilisateur
Character(30),PARAMETER :: fichconf = './BIO/config.ini'
Character(25),PARAMETER :: fichconfmod = './BIO/modele.def'
 REAL(8) ,PARAMETER      :: cini        = 0.001d0
 
 !--
 Character(20),PARAMETER :: rep_sortie = './BIO/SORTIES/'
 Character(10)           :: suffbio 
 Character(80)           :: resu_bio
 
 !-- Fichier des donnees d irradiance :
 Character(15) :: fichirrad
 Character(40) :: fichirrad_long 

!-- Repertoire des fichiers de conc initiales :
 Character(40), PARAMETER :: repconcinit='./BIO/DATA/CI/'
!-- Repertoire des fichiers de conc initiales :
 Character(40), PARAMETER :: repconcapp='./BIO/DATA/apports/'

!-- tableau contenant les indices des points du
!   maillage ou l''on veut sauvegarder les flux 
 integer, Allocatable           :: coord_flux(:,:)
 character(20),Allocatable      :: noms_fich_flux(:) 
 
End Module VAR_USER
!
!----------------------------------------------------------
!
                   Module VAR_GLOBAL
!		
! Module contenant des variables du modele biogeochimique
! Ces variables sont globales et donc accessibles dans toute 
! unite de programmation moyennant l instruction USE VAR_GLOBAL
!----------------------------------------------------------
USE DEF_TYPE
Implicit None

!-- nb de ss-compartiments de phyto et de zoo:
integer, allocatable :: BIO_C(:)        !tableau contenant les indices des biomasses carbonées dans VAR
integer, allocatable :: BIO_N(:)          !tableau contenant les indices des biomasses azotées dans VAR
integer, allocatable :: BIO_P(:)          !tableau contenant les indices des biomasses en P dans VAR
real(8), allocatable :: CBIO_C(:)        !tableau contenant les indices des biomasses carbonées dans VAR
real(8), allocatable :: CBIO_N(:)          !tableau contenant les indices des biomasses azotées dans VAR
real(8), allocatable :: CBIO_P(:)          !tableau contenant les indices des biomasses en P dans VAR
real(8), allocatable :: nb_bact(:,:)

Integer                 :: nscp_phy,nscp_zoo,nscp_bac
Integer, Allocatable    :: iscp_phy(:),iscp_zoo(:),iscp_bac(:)
Integer                 :: nbproc_bib,nbproc_flux, nbproc_self ! nb de processus modelises

!-- Rapport CHL:C (variable logique) 
Logical                 :: CHL_C_BOOL !faux  si rapport cstant 
                                      !vrai  si rapport variable
TYPE(VAR_GLOB_USER) ,Allocatable    :: CHL_C(:)
REAL(8)                 :: CHL_C0

!-- Sens d''orientation de la verticale donnee par la position de la maille nz_max:
 Character(4)::pos_nzmax

!-- grandeurs calculées  par le modèle, mais potentiellement utiles a d autres processus
!-- vitesse de production primaire mol C m-3 s-1 :
TYPE(VAR_GLOB_USER) ,Allocatable    ::  PPB(:)
!-- vitesse de production primaire mol C m-3 s-1 (NR = nutrient replete):
TYPE(VAR_GLOB_USER) ,Allocatable    ::  PPB_NR(:)
! vitesse de croissance specifique, en s-1
TYPE(VAR_GLOB_USER) ,Allocatable    ::  mu_phy(:)    !-- taux de croissance des phyto
! vitesse specifique de production primaire, en s-1
TYPE(VAR_GLOB_USER) ,Allocatable    ::  mu_PPB(:)
! vitesse specifique de production primaire (nutrient_replete), en s-1
TYPE(VAR_GLOB_USER) ,Allocatable    ::  mu_PPB_NR(:)   
! respiration du phyto en s-1
TYPE(VAR_GLOB_USER) ,Allocatable    ::  pr(:)
! vitesse d uptake de nitrate en mol N (mol C)-1 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: upno(:)
! vitesse d uptake d  ammonium en mol N (mol C)-1 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: upnh(:)
! vitesse d uptake du po4 en mol P (mol C)-1 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: uppo(:)
! fonction quota d uptake de nitrate sans unité
TYPE(VAR_GLOB_USER) , Allocatable  :: hQno(:)
! fonction quota d uptake d'ammonium sans unité
TYPE(VAR_GLOB_USER) , Allocatable  :: hQnh(:)
! fonction quota d uptake du po4 sans unité
TYPE(VAR_GLOB_USER) , Allocatable  :: hQpo(:)
! flux nh4 uptake den mol N m-3 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: fupnh(:)
! flux nO3 uptake den mol N m-3 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: fupno(:)
! flux po4 uptake den mol N m-3 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: fuppo(:)
! quota de N dans phyto
TYPE(VAR_GLOB_USER) , Allocatable  :: Qnh(:)
! quota de P dans phyto
TYPE(VAR_GLOB_USER) , Allocatable  :: Qpo(:)
! production bactérienne en carbone en µmol C l-1 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: bp(:)
! respiration bactérienne en carbone en µmol C l-1 s-1
TYPE(VAR_GLOB_USER) , Allocatable  :: br(:)
! vitesse specifique de broutage, s-1
TYPE(FLUX_GLOB_USER) ,Allocatable   ::  mu_graz(:,:)   
!limitation phyto N
TYPE(VAR_GLOB_USER) , Allocatable  :: limn(:)
!limitation phyto P
TYPE(VAR_GLOB_USER) , Allocatable  :: limp(:)
!limitation phyto finale
TYPE(VAR_GLOB_USER) , Allocatable  :: lim_nut(:)
!limitation bac C
TYPE(VAR_GLOB_USER) , Allocatable  :: limbaC(:)
!limitation bac N
TYPE(VAR_GLOB_USER) , Allocatable  :: limbaN(:)
!limitation bac P
TYPE(VAR_GLOB_USER) , Allocatable  :: limbaP(:)
!limitation bac finale
TYPE(VAR_GLOB_USER) , Allocatable  :: lim_ba(:)
!limitation temperature ppb
TYPE(VAR_GLOB_USER) , Allocatable  :: limT_ppb(:)
!limitation temperature grazing
TYPE(VAR_GLOB_USER) , Allocatable  :: limT_graz(:)
!limitation lumiere
TYPE(VAR_GLOB_USER) , Allocatable  :: lim_lum(:)
!remin nh4
TYPE(VAR_GLOB_USER) , Allocatable  :: remnh(:)
!remin po4
TYPE(VAR_GLOB_USER) , Allocatable  :: rempo(:)
End Module VAR_GLOBAL

!----------------------------------------------------------------
!
                  Module COUPLEUR_PHY_BIO
!
! Module contenant toutes les déclarations des variables utiles
! au couplage physique biogéochimie
!
!----------------------------------------------------------------
USE DEF_TYPE
Implicit None

!-- Repretoire contenant le code physique
 Character(20),PARAMETER :: rep_physique = './PHY/'
 
!-- Parametres definis par l utilisateur dans son fichier de config

Integer            :: nx_min,nx_max,ny_min,ny_max,nz_max  !dimensions spatiales
Integer            :: dx          !dimension de la maille
Integer            :: nbcomp,nbscomp, nbvar !nb de compartiments, de
                                  !sous-compart. et de variables
Integer            :: pHscale ! determination de l'echelle de pH
real(8)            :: duree,dt_bio!duree (j)/pas de temps de calcul de la biogeochimie (s)
real(8)            :: dt_save_bio!pas de tps de sauvegarde de la biogeochimie (mn)
real(8)            :: tps,tps_hr,dt_phy !tps de simul en secondes, en heures /pas de tps phy
real(8)            :: tps_hr0,tps_hr1 !heure de lancement et heure reelle
!-- parametres relatifs a l''irradiance:
real(8)            :: dt_irrad  ! pas de temps de lecture de l irradiance
real(8)            :: irr2par   ! facteur de conversion de l irradiance en PAR. Lui donner la valeur 1 dans config.ini 
                                ! si les donnees de lumiere sont deja exprimees en PAR 
real(8)            :: albedo, irr_param  ! albedo, param necessaire au calcul de PAR(0-)
real   , PARAMETER :: journee= 24*60*60 !--nb de secondes dans une journee

integer            :: nbcallbio   ! nbre d appels du module bio 

!-- declaration du nombre d element par organisme (ou sous-compartiment):
INTEGER, Allocatable :: nb_elmt(:)

!---------------------------------------------------------------------------
!-- declaration de variables physiques utilisees dans Eco3M
!---------------------------------------------------------------------------
!-- tableau des hauteurs d''eau en fonction de la position horizontale
REAL(8) ,Allocatable :: prof(:,:)

!-- altitude des mailles (en volume finis = distance des centres des volumes par rapport au fond)
REAL(8) ,Allocatable :: alti_z(:,:,:)

!--irradiance incidente a la surface de l''eau E(0+)
REAL(8) ,Allocatable :: irrad(:,:)

!--irradiance max (cas ou l''irradiance est calculee par une fonction)
REAL(8)                         :: irrad_MAX
!-- Photosynthetic Active Radiation (PAR)
REAL(8) ,Allocatable :: E_PAR(:,:)
!-- Photosynthetic Active Radiation (PAR) en fonction de la profondeur
REAL(8) ,Allocatable :: E_PARZ(:,:,:)


!-- parametres associes au modele d''extinction lumineuse dans
!   la colonne d''eau
 TYPE(PARAM)   :: IRR_PAR

!-- temperature
REAL(8)         ,Allocatable    :: TEMP_BIO(:,:,:)
!-- salinite
REAL(8)         ,Allocatable    :: SAL_BIO(:,:,:)
!-- vent
REAL(8)        ,Allocatable    :: WIND(:,:,:)
!-- co2atm
REAL(8)        ,Allocatable    :: CO2atm(:,:,:)
!-- rho_sw
REAL(8)        ,Allocatable    :: rho_sw(:,:,:)
!-- nh4
REAL(8)        ,Allocatable    :: NH4(:,:,:)
!-- no3
REAL(8)        ,Allocatable    :: NO3(:,:,:)
!-- po4
REAL(8)        ,Allocatable    :: PO4(:,:,:)
!---- concentration en  MES 
REAL(8)         ,Allocatable    :: MES(:,:,:)
!-- Parametres pour systeme des carbonates Calcul DIC
Real(8),Allocatable::DIC_OLD(:,:,:)
Real(8),Allocatable::K0(:,:,:)
Real(8),Allocatable::Ke(:,:,:)
Real(8),Allocatable::K1(:,:,:)
Real(8),Allocatable::K2(:,:,:)
Real(8),Allocatable::KB(:,:,:)
Real(8),Allocatable::Ksp(:,:,:)
Real(8),Allocatable::Ca2(:,:,:)
Real(8),Allocatable::Omega(:,:,:)
Real(8),Allocatable::AC(:,:,:)
Real(8),Allocatable::AC1(:,:,:)
Real(8),Allocatable::X(:,:,:)
Real(8),Allocatable::Y(:,:,:)
Real(8),Allocatable::H(:,:,:)
Real(8),Allocatable::OH(:,:,:)
Real(8),Allocatable::BOH4(:,:,:)
Real(8),Allocatable::CO2(:,:,:)
Real(8),Allocatable::CO3(:,:,:)
Real(8),Allocatable::HCO3(:,:,:)
Real(8),Allocatable::AT_test(:,:,:)
Real(8),Allocatable::TEMP(:,:,:)
REAL(8),Allocatable :: capp(:,:,:,:)
!---------------------------------------------------------------------------
!-- declaration des variables biogeochimiques
!---------------------------------------------------------------------------
! vecteur contenant l ensemble des variables d etat du modele
! biogeochimique :
TYPE(VAR_ETAT) ,Allocatable :: VAR(:)

!-- matrice des termes de tendance:
REAL(8) ,Allocatable :: TEND(:,:,:,:)

!---------------------------------------------------------------------------
!-- declaration des variables circulant entre les modeles physique et bio
!---------------------------------------------------------------------------
 REAL(8),Allocatable :: flux_phy(:,:,:,:)
 Character(70),Allocatable :: tab_noms_fich(:)
!-- Definition d un nouveau type (specifique au couplage avec la physique)
#ifdef COUPL
  #ifdef CALC
!-- Definition d un nouveau type (specifique au couplage avec la physique)
!   pour lequel la taille de la matrice contenant les dimensions dans l espace
!   est ajustee à la taille du maillage. Ce type est defini dans le fichier
!   "declar.inc" lequel est cree automatiquement lors de la compilation (mode INI)
!   lorsque Eco3M est couplé à un code physique:

    include "../declar.inc"

!-- tableaux des concentrations utilisees dans le code physique
    TYPE(VAR_COUPL),Allocatable    :: conc_phy0(:)
    TYPE(VAR_COUPL),Allocatable    :: conc_phy(:)
    TYPE(VAR_COUPL),Allocatable    :: conc_phy_old(:)
    Integer:: ics_save
  #endif
#endif


End Module COUPLEUR_PHY_BIO
!-----------------------------------------------------------------
