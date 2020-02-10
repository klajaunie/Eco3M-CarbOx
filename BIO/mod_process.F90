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
!---------------------------------------------------------------------------
!
                     MODULE F_PROCESS
!
!---------------------------------------------------------------------------
!
! Ce module sert a integrer toutes les fonctions modélisant
! des processus que l''on désire intégrer dans le modèle
!
! Ce fichier est donc l''un des rares fichiers que l''utilisateur est
! amené à modifier
!
! ATTENTION: il doit y figurer obligatoirement le nom de la fonction
! calculant le coefficient d''extinction lumineuse
!
! Derniere modification: 06/07/07
!-----------------------------------------------------------------
 Use COUPLEUR_PHY_BIO
 Use MOD_FCHAIN

 CONTAINS

!-----------------------------------------------------------------
! --------- Modèles de Production primaire Brute (PPB) -----------
!-----------------------------------------------------------------
!
!-----------------------------------------------------------------------------------
!---- Fonctions de quota de limitation de l''assimilation de nutriments (uptake) ----
!-----------------------------------------------------------------------------------

include "../F_PROCESS/LIGHT/ATTENUATION/f_extinc_morel.inc"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! AJOUT DES PROCESS CHRISTEL 9/10/2007
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
include "../F_PROCESS/PPB/f_ppb_geider_withPlim.inc"
include "../F_PROCESS/RESP/f_resp.inc"
include "../F_PROCESS/RESP/f_resp_bac.inc"
include "../F_PROCESS/REMIN/f_remin_kls.inc"
include "../F_PROCESS/PERTES/f_out.inc"
include "../F_PROCESS/PERTES/f_out_det.inc"
include "../F_PROCESS/PERTES/f_out_mod.inc"
include "../F_PROCESS/PERTES/f_p_zoo.inc"
include "../F_PROCESS/UPT/f_upt_droop.inc"
include "../F_PROCESS/UPT/f_upt_bact_droop_kls.inc"
include "../F_PROCESS/GRAZ/f_gra_phy.inc"
include "../F_PROCESS/GRAZ/f_gra_bact.inc"
include "../F_PROCESS/GRAZ/f_gra_det.inc"
include "../F_PROCESS/PBACT/f_bp_c.inc"
include "../F_PROCESS/REMIN/f_nit.inc"
include "../F_PROCESS/REMIN/f_prodO2.inc"
!AJOUT Marion 
!fonction de limitation du grazing par la temperature
include "../F_PROCESS/GRAZ/f_lim_graz.inc"
include "../F_PROCESS/PPB/f_limtemp_ppb.inc"
!AJOUT Katixa Lajaunie-Salla 11/06/2018
include "../F_PROCESS/AERA/f_aera.inc"
include "../F_PROCESS/CARBO/f_precipcarbo.inc"
include "../F_PROCESS/CARBO/f_disscarbo.inc"
END MODULE F_PROCESS
