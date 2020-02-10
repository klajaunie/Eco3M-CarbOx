#  MAKEFILE ECO3M - Makefile pour la première compilation (mode INI)
#
# CLES DE PRECOMPILATION à mettre dans CPPFLAGS:
#------------------------------------------------
# * CALC     : clé utilisée obligatoirement dans les fichiers Makefile.calc; correspond au mode calcul du modele
# * COUPL    : utilisation en mode couplé avec un modele physique
# * INI      : utilisée obligatoirement dans le fichier Makefile, pour lancer le mode d'initialisation qui permet de générer le modèle via la création du fichier call.inc
# * M0D-NCOUPL : clé qui ne peut se mettre que si on est en mode NON couplé (NCOUPL) et en 0D
# * M3D-NCOUPL : clé qui ne peut se mettre que si on est en mode NON couplé (NCOUPL) mais avec une geometrie 3D (ou 1D)
# * NCOUPL   : utilisation en mode NON couplé
# * MODTEST : en phase de test => affichages intermediaires
# * RESBIO_A : sortie des resultats au format ASCII
# * RESBIO_N : sortie des resultats au format NETCDF
# * VAR_GLOB: pour l'utilisation de variables globales supplémentaires pour le modèle biogéochimique
#
#Compilateur
FC            = ifort

#linker
LINKER        = ifort

#preprocesseur
CPP           = cpp -nostdinc
#
# FLAGS de precompilation
#------------------------
# 1/ cas d'une utilisation en mode couplé
#CPPFLAGS      = -P -C -DINI -DRESBIO_A -DCOUPL  -DMODTEST -I./BIO/INCLUDE/
# 2/ cas d'une utilisation en mode NON couplé en 0D 
CPPFLAGS      = -P -C -DINI -DRESBIO_A -DNCOUPL  -DM0D-NCOUPL  -I./BIO/INCLUDE/
#sauvegarde de 2/marion avant supression modtest
#CPPFLAGS      = -P -C -DINI -DRESBIO_A -DNCOUPL  -DM0D-NCOUPL -DMODTEST -I./BIO/INCLUDE/
# 3/ cas d'une utilisation en mode NON couplé en multi-dimension sans physique)
#CPPFLAGS      = -P -C -DINI -DRESBIO_A -DNCOUPL  -DM3D-NCOUPL -DMODTEST -I./BIO/INCLUDE/

# Repertoires INCLUDE
INCDIR        = -I./BIO/INCLUDE/ -I./BIO/MODULES/

# FLAG de compilation
COMP_OPTS     = ${INCDIR} ${MODDIR}
#
PROGRAM	      = ./eco3M_ini.exe

#LIBRAIRIES
# pour cluster
LIBS          =  ./BIO/LIB/libionetcdfPC.a \
		./BIO/LIB/libnetcdf_351.a
 
#REPERTOIRES
SRCBIO	      = ./BIO/
SRC           = ./BIO/petitf/
OBJ           = ./BIO/OBJETS/

#MODULES
MODULS	      = ${SRCBIO}mod_longchain.f90 \
		${SRCBIO}modules.f90 \
		${SRCBIO}mod_fchain.f90 \
		${SRCBIO}mod_process.f90 
#
# pour ifort :
 MODDIR        = -module ./BIO/MODULES/

#
OBJS          = ${OBJ}mod_longchain.o \
		${OBJ}modules.o \
		${OBJ}mod_fchain.o \
		${OBJ}mod_process.o \
		coupleur.o \
		${OBJ}main_bio.o \
		${OBJ}sub_affiche.o \
		${OBJ}sub_call_tab.o \
		${OBJ}sub_calc_prelim.o \
		${OBJ}sub_fin.o \
		${OBJ}sub_init.o \
		${OBJ}sub_init_mod.o \
		${OBJ}sub_var_ini.o \
                ${OBJ}sub_resbio_ascii.o
#		${OBJ}sub_resbio_netcdf.o
#
#
all:		 ${PROGRAM} 
#
# CREATION PROGRAMME - EDITIONS DES LIENS
#
modules:
	${FC} ${COMP_OPTS} -c ${MODULS} ${MODDIR}
	@echo "fin creation modules"
#
${PROGRAM}:	${OBJS}
	@echo  "CREATION ${PROGRAM} ... "
	${LINKER} ${INCDIR} ${LDFLAG} ${OBJS} ${LIBS}  -o  ${PROGRAM} 
	@echo "COMPILATION - lIENS TERMINE"
	${PROGRAM} 
	@echo "Mode INI TERMINE"
	rm -f  ${OBJ}* coupleur.o
# cas ou Eco3M n est pas couplé
	make -f Makefile.calc.noncouple
# cas ou Eco3M est couplé à un code physique
#	make -f Makefile.calc.couple

#
clean:;		rm -f  ${SRCBIO}modules.f90 coupleur_ini.f90 coupleur_calc.f90 coupleur.o ${OBJ}* ${SRC}* ./MODULES/*  ${PROGRAM} *.mod *.inc *.exe
#
#
# COMPILATION

# 
coupleur.o: coupleur.F90 	
	${CPP}  ${CPPFLAGS} coupleur.F90 coupleur_ini.f90	
	${FC}	${COMP_OPTS} -c coupleur_ini.f90 -o coupleur.o

${OBJ}main_bio.o: ${SRCBIO}main_bio.F90 
	${CPP}  ${CPPFLAGS} ${SRCBIO}main_bio.F90 ${SRC}main_bio_ini.f90	
	${FC}	${COMP_OPTS} -c ${SRC}main_bio_ini.f90 -o ${OBJ}main_bio.o 

${OBJ}mod_longchain.o: ${SRCBIO}mod_longchain.f90 
	${FC}	${COMP_OPTS} -c ${SRCBIO}mod_longchain.f90 -o ${OBJ}mod_longchain.o
# 
${OBJ}mod_fchain.o: ${SRCBIO}mod_fchain.f90 	
	${FC}	${COMP_OPTS} -c ${SRCBIO}mod_fchain.f90 -o ${OBJ}mod_fchain.o
# 
${OBJ}mod_process.o: ${SRCBIO}mod_process.F90 	
	${CPP}  ${CPPFLAGS} ${SRCBIO}mod_process.F90 ${SRC}mod_process.f90	
	${FC}	${COMP_OPTS} -c ${SRC}mod_process.f90 -o ${OBJ}mod_process.o
# 
${OBJ}modules.o: ${SRCBIO}modules.F90 	
	${CPP}  ${CPPFLAGS} ${SRCBIO}modules.F90 ${SRCBIO}modules.f90	
	${FC}	${COMP_OPTS} -c ${SRCBIO}modules.f90 -o ${OBJ}modules.o
# 
${OBJ}sub_affiche.o: ${SRCBIO}sub_affiche.F90 	
	${CPP}  ${CPPFLAGS} ${SRCBIO}sub_affiche.F90 ${SRC}sub_affiche.f90	
	${FC}	${COMP_OPTS} -c ${SRC}sub_affiche.f90 -o ${OBJ}sub_affiche.o
# 
${OBJ}sub_calc_prelim.o: ${SRCBIO}sub_calc_prelim.F90 	
	${CPP}  ${CPPFLAGS} ${SRCBIO}sub_calc_prelim.F90 ${SRC}sub_calc_prelim.f90	
	${FC}	${COMP_OPTS} -c ${SRC}sub_calc_prelim.f90 -o ${OBJ}sub_calc_prelim.o
# 
 ${OBJ}sub_call_tab.o: ${SRCBIO}sub_call_tab.f90 	
	${FC}	${COMP_OPTS} -c ${SRCBIO}sub_call_tab.f90 -o ${OBJ}sub_call_tab.o
# 
${OBJ}sub_fin.o: ${SRCBIO}sub_fin.f90 	
	${FC}	${COMP_OPTS} -c ${SRCBIO}sub_fin.f90 -o ${OBJ}sub_fin.o
# 
${OBJ}sub_init.o: ${SRCBIO}sub_init.F90 	
	${CPP}  ${CPPFLAGS} ${SRCBIO}sub_init.F90 ${SRC}sub_init_ini.f90	
	${FC}	${COMP_OPTS} -c ${SRC}sub_init_ini.f90 -o ${OBJ}sub_init.o
# 
${OBJ}sub_init_mod.o: ${SRCBIO}sub_init_mod.f90 	
	${FC}	${COMP_OPTS} -c ${SRCBIO}sub_init_mod.f90 -o ${OBJ}sub_init_mod.o
# 
#${OBJ}sub_phys_ini.o: ${SRCBIO}sub_phys_ini.f90 	
#	${FC}	${COMP_OPTS} -c ${SRCBIO}sub_phys_ini.f90 -o ${OBJ}sub_phys_ini.o
# 
${OBJ}sub_var_ini.o: ${SRCBIO}sub_var_ini.F90 	
	${CPP}  ${CPPFLAGS} ${SRCBIO}sub_var_ini.F90 ${SRC}sub_var_ini.f90	
	${FC}	${COMP_OPTS} -c ${SRC}sub_var_ini.f90 -o ${OBJ}sub_var_ini.o
#	
${OBJ}sub_resbio_ascii.o: ${SRCBIO}sub_resbio_ascii.F90 
	${CPP}  ${CPPFLAGS} ${SRCBIO}sub_resbio_ascii.F90 ${SRC}sub_resbio_ascii.f90	
	${FC}   ${COMP_OPTS} -c ${SRC}sub_resbio_ascii.f90 -o ${OBJ}sub_resbio_ascii.o
# 
#${OBJ}sub_resbio_netcdf.o: ${SRCBIO}sub_resbio_netcdf.f90 
#	${FC}   ${COMP_OPTS} -c ${SRCBIO}sub_resbio_netcdf.f90 -o ${OBJ}sub_resbio_netcdf.o
