#!/usr/bin/python2.7
from math import *
import string

#####################
## FONCTIONS PARSE ##
#####################

def parsePDB(infile) :
	# lecture du fichier PDB 
	f = open(infile, "r")
	lines = f.readlines()
	f.close()
	# initialisation de variables
	conf= True
	firstline = True
	prevres = None
	dPDB = {}
	dPDB["conformation"] = []
	dPDB["liste_resid"] = []
	# parcoure le fichier PDB
	for line in lines :
		if line[0:5] == "MODEL" :
			conf=(line[10:15]).strip() #attention entier sous forme de chaine de char
			if not conf in dPDB["conformation"] :
				dPDB["conformation"].append(conf)
				dPDB[conf] = {}
				dPDB[conf]["liste_resid"] = []
		elif line[0:4] == "ATOM" :
			curres = "%s"%(line[22:26]).strip()
			if not curres in dPDB[conf]["liste_resid"] :
				dPDB[conf]["liste_resid"].append(curres)
				dPDB[conf][curres] = {}
				dPDB[conf][curres]["nom_resid"] = string.strip(line[17:20])
				dPDB[conf][curres]["liste_atom"] = []
			atomtype = string.strip(line[12:16])
			dPDB[conf][curres]["liste_atom"].append(atomtype)
			dPDB[conf][curres][atomtype] = {}
			dPDB[conf][curres][atomtype]["x"] = float(line[30:38])
			dPDB[conf][curres][atomtype]["y"] = float(line[38:46])
			dPDB[conf][curres][atomtype]["z"] = float(line[46:54])
			dPDB[conf][curres][atomtype]["id"] = line[6:11].strip()
	return dPDB
    
def writePDB(dPDB, filout = "out.pdb") :
	#ecrit les informations globales contenu dans le dictionnaire (a titre indicatif) dans un fichier PDB
	fout = open(filout, "w")
	for conf in dPDB["conformation"]:
		fout.write("MODEL     %4s	RMSD %3.3f	Rgiration=%3.3f\n"%(conf,dPDB[conf]["RMSD"],dPDB[conf]["rayon_giration"]))
		for res in dPDB[conf]["liste_resid"] :
			for atom in dPDB[conf][res]["liste_atom"] :
				fout.write("ATOM  %4s %4s %4s %3s %4s  %8.3f%8.3f%8.3f  %8.3f \n"%(conf, dPDB[conf][res][atom]["id"], atom, res,dPDB[conf][res]["nom_resid"],dPDB[conf][res][atom]["x"], dPDB[conf][res][atom]["y"],dPDB[conf][res][atom]["z"],dPDB[conf][res][atom]["dist_CM"]))
		fout.write("TER\nENDMDL\n")
	fout.close()

def writeglobal(dico, filout = "resultat_global.txt"):
	#ecrit le RMSD et RG de chaque conformation dans un fichier .txt (analyse globale)
	fout = open(filout, "w")
	for conf in dico["conformation"]:
		fout.write("MODEL     %4s	RMSD %3.3f	Rgiration=%3.3f\n"%(conf,dico[conf]["RMSD"],dico[conf]["rayon_giration"]))
	fout.close()

def writelocal (dico, filout = "resultat_local.txt"):
	#ecrit le RMSD et RG de chaque residu, pour chaque conformation (en dehors de la "0" de reference) dans un fichier .txt
	fout = open(filout, "w")
	for conf in dico["conformation"]:
		if conf != "0" :	#attention conf est un string
			fout.write("MODEL     %4s\n"%(conf))
			for res in dico[conf]["liste_resid"]:
				fout.write("\t %3s%4s : RMSD:%3.3f\tdist_res_CM: %3.3f\n"%(dico[conf][res]["nom_resid"],res,dico[conf][res]["RMSDlocal"],dico[conf][res]["dist_res_CM"]))
	fout.close()

###########################
## FONCTIONS PRINCIPALES ##
###########################

def CM_Conf(dPDB,all=True,reslist = False):
	#calcul le centre de masse de la chaque conformation
	if all == True :
		reslist = dPDB["liste_resid"]
	xconf = yconf = zconf = 0.0
	nb_atome_tot=0
	for res in reslist :
		for atom in dPDB[res]["liste_atom"] :
			nb_atome_tot+=1
			xconf +=dPDB[res][atom]["x"]
			yconf +=dPDB[res][atom]["y"]
			zconf +=dPDB[res][atom]["z"]
	Xcm = float(xconf)/nb_atome_tot
	Ycm = float(yconf)/nb_atome_tot
	Zcm = float(zconf)/nb_atome_tot
	dPDB["XCM"] = Xcm
	dPDB["YCM"] = Ycm
	dPDB["ZCM"] = Zcm

def dist_ATOM_CM (dresidu, xconf, yconf, zconf):
	#calcul la distance de chaque atome, par rapport au centre de masse de la conformation
	for atom in dresidu["liste_atom"] :
		dx = dy = dz = 0.0
		dx=xconf-dresidu[atom]["x"]
		dy=yconf-dresidu[atom]["y"]
		dz=zconf-dresidu[atom]["z"]
		dresidu[atom]["dist_CM"]=sqrt(dx**2+dy**2+dz**2)

def CM_Res(dresidu):
	#calcul du centre de masse d'un residu 
		x = y = z = 0.0
		for atom in dresidu["liste_atom"] :
			x +=dresidu[atom]["x"]
			y +=dresidu[atom]["y"]
			z +=dresidu[atom]["z"]
		Xcm = float(x)/len(dresidu["liste_atom"]) 
		Ycm = float(y)/len(dresidu["liste_atom"])
		Zcm = float(z)/len(dresidu["liste_atom"])
		dresidu["XRES"] = Xcm
		dresidu["YRES"] = Ycm
		dresidu["ZRES"] = Zcm

def recherche_point_distant(dconf):
	#retourne le point la distance la plus grande avec le CM pour pour chaque conformation --> Rayon de Giration
	distancemax=0.0
	for res in dconf["liste_resid"]:
		for atom in dconf[res]["liste_atom"]:
			if dconf[res][atom]["dist_CM"]>distancemax:
				distancemax=dconf[res][atom]["dist_CM"]
	return distancemax

def RMSD_conf_global (dico,dicoref,listeresidu):
	#calcul du RMSD (global) pour chaque conformation
	for conf in dico["conformation"]:
		RMSD=0.0
		nb_atom=0 #nb d'atomes dans une conformation
		for res in listeresidu: #on boucle sur l'ensemble des residu de chaque conformation.
			for atom in dico[conf][res]["liste_atom"]:
				nb_atom+=1
				#somme des carres des variations par rapport a l'origine
				RMSD+=((dico[conf][res][atom]["x"] - dicoref["0"][res][atom]["x"])**2+(dico[conf][res][atom]["y"] - dicoref["0"][res][atom]["y"])**2 + (dico[conf][res][atom]["z"] - dicoref["0"][res][atom]["z"])**2)
		dico[conf]["RMSD"]=sqrt(RMSD/nb_atom)
		
def dist_RES_CM (dresidu, xconf, yconf, zconf):
	#calcul de la distance entre le centre de masse d'un residu et le centre de masse de la conformation.
	dx = dy = dz = 0.0
	dx=xconf-dresidu["XRES"]
	dy=yconf-dresidu["YRES"]
	dz=zconf-dresidu["ZRES"]
	dresidu["dist_res_CM"]=sqrt(dx**2+dy**2+dz**2)

def RMSD_res_local(dresid,dicorefresid):
	#calcul du RMSD pour chaque residu
	RMSD=0.0
	nb_atom=0;
	for atom in dresid["liste_atom"]: #pour chaque atome, compteur++
		nb_atom+=1
		#somme des carres des variations par rapport a l'origine
		RMSD+=((dresid[atom]["x"] - dicorefresid[atom]["x"])**2+(dresid[atom]["y"] - dicorefresid[atom]["y"])**2 + (dresid[atom]["z"] - dicorefresid[atom]["z"])**2)
	dresid["RMSDlocal"]=sqrt(RMSD/(nb_atom))
