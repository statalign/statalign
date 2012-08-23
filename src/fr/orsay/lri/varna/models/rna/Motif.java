/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 Universit� Paris-Sud 91405 Orsay Cedex France

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with VARNA version 3.1.
 If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.models.rna;
import java.awt.geom.Point2D;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.lang.Math;





public class Motif implements Serializable{
	StructureTemp _listStrands = new StructureTemp();
	private Double _spaceBetweenBases = 1.0;
	private ArrayList<ModeleBase> _listeBasesM;
	private int _ajustement;
	private double _decalage;
	private RNA _rna;
	
	
	public Motif(RNA rna, ArrayList<ModeleBase> listeBases){		
		this._listeBasesM = listeBases;		
		_rna = rna;
		
	}
	
	
	
	/**
	 * Find all the strands and store them in an arrayList.
	 */
	public void listStrand() {
		int i=0;
		int k=0; //compteur de brin
		int indice_min = 0;
		int indice_max = 0;
		
		getListStrand().clearListStrands();
		
		
		while ( i < getListBasesMotif().size()-1){
			int a =getListBasesMotif().get(i).getBaseNumber();
			int b = getListBasesMotif().get(i+1).getBaseNumber();
			if((b-a)!=1){
				boolean orientation = false;
				indice_max=i;
				ModeleStrand strand = new ModeleStrand();
				for(int j=indice_min; j<=indice_max;j++){
					strand.addBase(getListBasesMotif().get(j));
					getListBasesMotif().get(j).setNumStrand(k);
					
					if(getListStrand().sizeStruct()%2==0 || getListStrand().isEmpty()){
						orientation=true;
					}										
				}
				if(!orientation){
					Collections.sort(strand.getArrayListMB(), Collections.reverseOrder());
				}
				getListStrand().addStrand(strand);
				k++;
				indice_min=i+1;
			}
			i++;
		}
		//On s'occupe du dernier brin. Ou du seul brin de la structure.
		ModeleStrand strand = new ModeleStrand();
		//System.out.print(getListStrand().sizeStruct()+"\n");
		boolean orientation = false;
		for(int j=indice_min;j<getListBasesMotif().size();j++){
			strand.addBase(getListBasesMotif().get(j));
			getListBasesMotif().get(j).setNumStrand(k);
			if(getListStrand().sizeStruct()%2==0 || getListStrand().isEmpty()){
				orientation=true;
			}
		}
		if(!orientation){
			Collections.sort(strand.getArrayListMB(), Collections.reverseOrder());
		}
		getListStrand().addStrand(strand);

	}
	
	public void decalerBase(ModeleStrand strand, double decalage){
		for (int j = 0; j < strand.sizeStrand(); j++ ){
			int indice = strand.getMB(j).getIndex();
			if(!strand.getMB(j).getChecker()){
				getListBasesMotif().get(indice).setCoords(
						new Point2D.Double(getListBasesMotif().get(indice).getCoords().x,
								getListBasesMotif().get(indice).getCoords().y+decalage));
			}
		}
	}
	
	public void initChecker(ModeleStrand strand){
		for (int j = 0; j < strand.sizeStrand(); j++ ){
			strand.getMB(j).setChecker(false);
		}			
	}
	
	public void initCheckerAll(){
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			for (int j = 0; j < getListStrand().getStrand(i).sizeStrand(); j++ ){
				getListStrand().getStrand(i).getMB(j).setChecker(false);
			}		
		}
	}
	
	public void initStrand(){
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			getListStrand().getStrand(i).setStrandLeft(false);
			getListStrand().getStrand(i).setStrandRight(false);
			getListStrand().getStrand(i).setHasBeenPlaced(false);
			getListStrand().getStrand(i).setLevelPosition(1);
		}		
	}
	
	public void initInterv(){
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			for (int j = 0; j < getListStrand().getStrand(i).sizeStrand(); j++ ){
				getListStrand().getStrand(i).getMB(j).setIntervDroite(false);
				getListStrand().getStrand(i).getMB(j).setIntervGauche(false);
			}
		}
	}
	
	/**
	 * Set the coordinates of all bases for all strands in the list of strands
	 */
	public void positionneStrand(){
		int spaceBetweenStrands = 0;
		setAjustement(0);
		boolean droite = true;
		//On boucle sur la liste des brins
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			this.initChecker(getListStrand().getStrand(i));
			if(!this.getListStrand().getStrand(i).getStrandLeft()&& !this.getListStrand().getStrand(i).getStrandRight()){
				if(droite){
					this.getListStrand().getStrand(i).setStrandLeft(true);
					droite=false;
				}
				else{
					this.getListStrand().getStrand(i).setStrandRight(true);
					
					droite=true;
				}
			}
			this.positionneSpecificStrand(i, spaceBetweenStrands);
			spaceBetweenStrands=spaceBetweenStrands+60;
			//S'il ya eu au moins un reajustement intrabrin, il faut decaler davantage le brin suivant
			if(getAjustement()>0){
				spaceBetweenStrands=2*spaceBetweenStrands;
			}
		}
	}
	
	/**
	 * Set the coordinates of all bases for one strand in particular
	 * @param i the index of the strand in the list of strands
	 * @param d the space which separate two strands
	 */
	public void positionneSpecificStrand(int i, double d){
		int dist3=-1;
		int dist5=-1;
		int part3Temp =-1;
		int part5Temp=-1;
			this.initChecker(getListStrand().getStrand(i));
			for (int j = 0; j < getListStrand().getStrand(i).sizeStrand(); j++ ){
				int indice = getListStrand().getStrand(i).getMB(j).getIndex();
				int part3 = getListBasesMotif().get(indice).getStyleBP().getPartner3().getIndex();
				int part5 = getListBasesMotif().get(indice).getStyleBP().getPartner5().getIndex();
				int dist=0;
				int boucle =0;
				int max=0;
				ModeleBase mb3=getListBasesMotif().get(indice).getStyleBP().getPartner3();
				ModeleBase mb5=getListBasesMotif().get(indice).getStyleBP().getPartner5();
				//On positionne la première base du brin
				if(!getListBasesMotif().get(indice).getChecker() && j==0){
					getListBasesMotif().get(indice).setCoords(
							new Point2D.Double(d,0));
				}
				
				//On positionne les autres bases du brin, les unes par rapport aux autres
				if(!getListBasesMotif().get(indice).getChecker() && j!=0){
					dist3=Math.abs(part3Temp-part3);
					dist5=Math.abs(part5Temp-part5);
					
					if(Math.abs(dist3-dist5)>=0 && dist3>=dist5 && part3!=part5 && part3Temp!=-1 && part5Temp!=-1){
						dist=Math.abs(dist3-dist5);
						boucle++;
					}
					getListBasesMotif().get(indice).setCoords(
							new Point2D.Double(d,getListStrand().getStrand(i).getMB(j-1).getCoords().y
									+(dist+1)*_spaceBetweenBases * 50));
				}
				
				//Dans un meme brin, si la paire de base a ete traitee 
				//On reajuste la structure pour eviter d'avoir des angles non droits
				if(getListBasesMotif().get(indice).getChecker() && 
						!getListBasesMotif().get(indice-1).getChecker()){
					int indiceP = getListStrand().getStrand(i).getMB(j-1).getIndex();	
					int indicePP=getListStrand().getStrand(i).getMB(j-2).getIndex();
					
					getListBasesMotif().get(indiceP).setCoords(
							new Point2D.Double(getListBasesMotif().get(indice).getCoords().x,getListBasesMotif().get(indicePP).getCoords().y));
				}
				
				//On verifie si la base est impliquee dans une paire intra-brin
				//On reajuste les coordonnees si c'est le cas :
				//On place le partenaire 3' en face de son partenaire 5'
				if((getListStrand().getStrand(i).existInStrand(part3)&&getListStrand().getStrand(i).existInStrand(part5))&&
						part3!=part5 && getListBasesMotif().get(indice).getStyleBP().isCanonical()
					&& !mb3.getChecker()){
					getListBasesMotif().get(part3).setCoords(new Point2D.Double
							(getListBasesMotif().get(part5).getCoords().x +
									_spaceBetweenBases * 50,getListBasesMotif().get(part5).getCoords().y));
					//positionner les nt qui constituent les boucles.
					if((dist>=0 && boucle>0) || (part3Temp==-1 && part3!=part5)){
						max=part3Temp;
						if(part3Temp==-1){
							max = getListStrand().getStrand(i).sizeStrand();
						}						
						for(int k=part3+1; k<max; k++){
							getListBasesMotif().get(k).setCoords(
									new Point2D.Double(getListBasesMotif().get(k-1).getCoords().x
											,getListBasesMotif().get(k-1).getCoords().y- _spaceBetweenBases * 50));
							getListBasesMotif().get(k).setChecker(true);
						}
					}
					part3Temp=part3;
					part5Temp=part5;
					setAjustement(getAjustement()+1);
					mb3.setChecker(true);
					mb5.setChecker(true);								
				}
			}
			
	}
	
	
	public void ajusteStrand(){		
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			ajusteSpecificStrand(i);
		}		
	}
	
	public void ajusteSpecificStrand(int i){					
		for (int j = 0; j < this.getListStrand().getStrand(i).sizeStrand(); j++ ){
			int indice = this.getListStrand().getStrand(i).getMB(j).getIndex();
			int partner=-1;
			this.initChecker(this.getListStrand().getStrand(i));
			ArrayList<ModeleBase> mb = _rna.getAllPartners(indice);
			// A changer plus tard pour traiter les paires de base multiples...
			
			
			if (mb.size()>0){
				partner = mb.get(0).getIndex();
			}
			this.setAjustement(0);
			this.setDecalage(0);
				
			//On verifie si la base est impliquee dans une paire inter-brin
			//On reajuste les coordonnees si c'est le cas
			if((!getListStrand().getStrand(i).existInStrand(partner)||!getListStrand().getStrand(i).existInStrand(indice))&&
					partner!=-1 && !mb.get(0).getChecker()){
				//On change les coordonnées du nucleotide present sur le brin actuel
					setDecalage(getListBasesMotif().get(partner).getCoords().y - getListBasesMotif().get(indice).getCoords().y);
					
					this.getListBasesMotif().get(indice).setCoords(new Point2D.Double
							(this.getListBasesMotif().get(indice).getCoords().x ,this.getListBasesMotif().get(partner).getCoords().y));
						
				
				 this.setAjustement(getAjustement()+1);
				 this.getListStrand().getStrand(i).getMB(j).setChecker(true);
				 this.getListBasesMotif().get(partner).setChecker(true);								
			}
			
			//decaler les autres bases du brin s'il y a eu reajustement
			if(this.getAjustement()!=0 && this.getDecalage()!=0){
				this.decalerBase(this.getListStrand().getStrand(i), this.getDecalage());
			}				
		}			
	}
	
	public void reajustement(){
		this.initCheckerAll();
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			int partnerLast = -1;
			int partnerFirst = -1;
			int numStrandPartner = -1;
			double decalage = -1;
			int size = getListStrand().getStrand(i).sizeStrand();
			int last = getListStrand().getStrand(i).getMB(size-1).getIndex();
			int first = getListStrand().getStrand(i).getMB(0).getIndex();
			ArrayList<ModeleBase> mbLast = _rna.getAllPartners(last);
			ArrayList<ModeleBase> mbfirst = _rna.getAllPartners(first);
			int strandModified = -1; //brin sur lequel les changements de coord ont t effectu
			if (mbLast.size()>0 && mbfirst.size()>0){
				partnerLast = mbLast.get(0).getIndex();
				numStrandPartner = mbLast.get(0).getNumStrand();
				partnerFirst = mbfirst.get(0).getIndex();
			}

			if(partnerLast !=-1 && partnerFirst !=-1){
				decalage = Math.abs(this.getListBasesMotif().get(partnerFirst).getCoords().y
						-this.getListBasesMotif().get(first).getCoords().y);
				if(decalage == 0){
					decalage = Math.abs(this.getListBasesMotif().get(partnerLast).getCoords().y
							-this.getListBasesMotif().get(last).getCoords().y);
				}
				if(size > getListStrand().getStrand(numStrandPartner).sizeStrand()){
					//On modifiera les coord des bases prsentes sur le brin le plus court
					this.getListBasesMotif().get(partnerLast).setCoords(new Point2D.Double(
							this.getListBasesMotif().get(partnerLast).getCoords().x,
							this.getListBasesMotif().get(last).getCoords().y));

					this.getListBasesMotif().get(partnerFirst).setCoords(new Point2D.Double(
							this.getListBasesMotif().get(partnerFirst).getCoords().x,
							this.getListBasesMotif().get(first).getCoords().y));
					strandModified = numStrandPartner;

				}
				else{
					this.getListBasesMotif().get(last).setCoords(new Point2D.Double(
							this.getListBasesMotif().get(last).getCoords().x,
							this.getListBasesMotif().get(partnerLast).getCoords().y));

					this.getListBasesMotif().get(first).setCoords(new Point2D.Double(
							this.getListBasesMotif().get(first).getCoords().x,
							this.getListBasesMotif().get(partnerFirst).getCoords().y));
					strandModified = i;
				}

				
			}
			
			for (int j = 0; j < this.getListStrand().getStrand(i).sizeStrand(); j++ ){
				int indice = this.getListStrand().getStrand(i).getMB(j).getIndex();
				int partner = -1;
				ArrayList<ModeleBase> amb = _rna.getAllPartners(indice);
				
				if (amb.size()>0){
					partner = amb.get(0).getIndex();
				}
				
				if (partner!=-1 && !this.getListBasesMotif().get(indice).getChecker()){
					double dist = Math.abs(this.getListBasesMotif().get(indice).getCoords().y 
							- this.getListBasesMotif().get(partner).getCoords().y );
					double ccordxIndice = this.getListBasesMotif().get(indice).getCoords().x;
					double coordxPartner = this.getListBasesMotif().get(partner).getCoords().x;
					if(dist <= decalage && ccordxIndice != coordxPartner ){
						if(strandModified == amb.get(0).getNumStrand()){
							this.getListBasesMotif().get(partner).setCoords(new Point2D.Double(
									this.getListBasesMotif().get(partner).getCoords().x,
									this.getListBasesMotif().get(indice).getCoords().y));

						}
						else{
							this.getListBasesMotif().get(indice).setCoords(new Point2D.Double(
									this.getListBasesMotif().get(indice).getCoords().x,
									this.getListBasesMotif().get(partner).getCoords().y));
						}
					}

				}
				
				if(partner!=-1){
					this.getListBasesMotif().get(partner).setChecker(true);
				}
				
				this.getListBasesMotif().get(indice).setChecker(true);
			}
		}
	}
	
	
	
	/**
	 * Find the central strand in a RNA star pattern
	 * @return the index of the strand which is the central strand in a star pattern, or -1 if there is no central strand.
	 */
	public int getCentralStrand(){
		int partnerStrand=-1;
		int partner=-1;
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			int centralStrand =i;
			int k=0;
			for (int j = 0; j < getListStrand().getStrand(i).sizeStrand(); j++ ){
				
				int indice = getListStrand().getStrand(i).getMB(j).getIndex();
				ArrayList<ModeleBase> mb = _rna.getAllPartners(indice);
				// A changer plus tard pour traiter les paires de base multiples...
				
				
				if (mb.size()>0){
					partner = mb.get(0).getIndex();
					partnerStrand=getListBasesMotif().get(partner).getNumStrand();
					if(partner!=-1 && partnerStrand!=i){			
						
						if(centralStrand!=partnerStrand){
							centralStrand=partnerStrand;
							k++;
						}
						if(k>1){
							centralStrand = i;
							return centralStrand;
						}					
					}			
				}			
			}
		}
		return -1;
	
	}
	
	public void orderStrands(int centralStrand){
		int partner=-1;
		boolean droite = false;	
		this.initStrand();
		this.initInterv();
		//this.getListStrand().getStrand(centralStrand).setStrandRight(true);
		Hashtable<Integer,Integer> htableLast = new Hashtable<Integer,Integer>();
		Hashtable<Integer,Integer> htableFirst = new Hashtable<Integer,Integer>();
		
		//Remplissage des tables de hashage 
		//htableLast : les clés correspondent aux brins, les valeurs correspondent 
		//à la dernière base du brin central en interaction avec un autre brin
		//htableFirst : les valeurs correspondent à la permière base du brin central
		//en interaction avec un autre brin
		
		// On boucle sur les bases du brin central
		for(int j = 0; j < this.getListStrand().getStrand(centralStrand).sizeStrand(); j++ ){
			int indice = this.getListStrand().getStrand(centralStrand).getMB(j).getIndex();
			ArrayList<ModeleBase> mb = _rna.getAllPartners(indice);
				
			if (mb.size()>0){
				partner = mb.get(0).getIndex();
			}
			// Si la base est liée à un autre brin
			if (partner!=-1) {
				if (!htableFirst.containsKey(this.getListBasesMotif().get(partner).getNumStrand())) { 
					htableFirst.put(this.getListBasesMotif().get(partner).getNumStrand(), j);
				}		
				htableLast.put(this.getListBasesMotif().get(partner).getNumStrand(), j);
				
			}
		}
		
		
		//int comptDroite=0; //permet de savoir si un brin a été placé à droite.
		//int comptGauche=0;
		int previousStrand =-1;
		//On boucle sur les bases du brin central
		for (int j = 0; j < this.getListStrand().getStrand(centralStrand).sizeStrand(); j++ ){
			int indice = this.getListStrand().getStrand(centralStrand).getMB(j).getIndex();
			ModeleBase baseStrand = this.getListStrand().getStrand(centralStrand).getMB(j);
			//récupération du partenaire
			// A changer plus tard pour traiter les paires de base multiples...
			
			ArrayList<ModeleBase> mb = _rna.getAllPartners(indice); 
			int indexPartenaire = -1;
			ModeleBase partenaire = null;
			if (mb.size()>0)
			{
				indexPartenaire = mb.get(0).getIndex();
				partenaire = mb.get(0);
			}
			double space = _spaceBetweenBases*60;
			if(indexPartenaire!=-1 && !this.getListStrand().getStrand(partenaire.getNumStrand()).hasBeenPlaced()){
				
				//Si la base n'est pas comprise dans un intervalle à droit et à gauche
				if(!baseStrand.getIntervDroite() && !baseStrand.getIntervGauche()){
				
					if(!droite){
						droite=true;
						getListStrand().getStrand(partenaire.getNumStrand()).setStrandRight(true);
						//comptDroite++;
					}
				
					else{
						droite=false;
						space=-space;
						getListStrand().getStrand(partenaire.getNumStrand()).setStrandLeft(true);
						//comptGauche++;
					}
				}
				
				//Si la base est comprise dans un intervalle à droite, on place le brin à gauche
				else if(baseStrand.getIntervDroite() && !baseStrand.getIntervGauche()){
					space=-space;
					getListStrand().getStrand(partenaire.getNumStrand()).setStrandLeft(true);
					
				}
				
				//Si la base est comprise dans un intervalle à gauche, on place le brin à droite
				else if(!baseStrand.getIntervDroite() && baseStrand.getIntervGauche()){
					getListStrand().getStrand(partenaire.getNumStrand()).setStrandRight(true);
					
				}
				
				//Cas ou la base est comprise à la fois dans un intervalle à gauche et à droite
				//On place arbitrairement le brin (a gauche)
				else if (baseStrand.getIntervDroite() && baseStrand.getIntervGauche()){
					space=-space;
					getListStrand().getStrand(partenaire.getNumStrand()).setLevelPosition(2);
					getListStrand().getStrand(partenaire.getNumStrand()).setStrandLeft(true);
				}
				
				space = getListStrand().getStrand(partenaire.getNumStrand()).getLevelPosition()*space;
				
				
				
				for(int i=htableFirst.get(partenaire.getNumStrand()); i < htableLast.get(partenaire.getNumStrand()); i++){
					if(space<0){
						this.getListStrand().getStrand(centralStrand).getMB(i).setIntervGauche(true);
					}
					else{
						this.getListStrand().getStrand(centralStrand).getMB(i).setIntervDroite(true);
					}
				}
				this.positionneSpecificStrand(partenaire.getNumStrand(), space);
				this.ajusteSpecificStrand(partenaire.getNumStrand());
				
				if(previousStrand!=-1 && previousStrand!=partenaire.getNumStrand()){
					//System.out.println("TEST1");
					
					if(getListStrand().getStrand(previousStrand).getStrandLeft()&& 
							getListStrand().getStrand(partenaire.getNumStrand()).getStrandLeft()||
							getListStrand().getStrand(previousStrand).getStrandRight()&& 
							getListStrand().getStrand(partenaire.getNumStrand()).getStrandRight()
							&& getListStrand().getStrand(partenaire.getNumStrand()).getLevelPosition() == 
								getListStrand().getStrand(previousStrand).getLevelPosition()){
						//System.out.println("TEST2");
						
						ArrayList<ModeleBase> mbPartner = _rna.getAllPartners(htableLast.get(previousStrand)); 
						int indexPreviousPartenaire = mbPartner.get(0).getIndex();
						
						//coord Y de la derniere base du brin prcdent en intraction avec le brin central
						double coordYLastBase= getListStrand().getStrand(centralStrand).
							getMB(htableLast.get(previousStrand)).getCoords().y;
						
						//coord Y de la premiere base du brin actuel en interaction avec le brin central
						double coordYFirstBase= getListStrand().getStrand(centralStrand).
							getMB(htableFirst.get(partenaire.getNumStrand())).getCoords().y;
						
						double dist = Math.abs(coordYLastBase-coordYFirstBase);
						int sizeStrand = getListStrand().getStrand(partenaire.getNumStrand()).sizeStrand();
						
						//distance entre la 1ere base du brin partenaire en interaction avec le brin central et 
						double distC = Math.abs(this.getListBasesMotif().get(indexPartenaire).getCoords().y 
								-getListStrand().getStrand(partenaire.getNumStrand()).getMB(sizeStrand-1).getCoords().y);
						
						//distance entre la derniere base du brin prcdent en interaction avec le brin central et 
						//et la 1ere base du brin 
						double distB = Math.abs(this.getListBasesMotif().get(indexPreviousPartenaire).getCoords().y 
								-getListStrand().getStrand(previousStrand).getMB(0).getCoords().y);
						
						//System.out.println("DistC:"+distC+" distB:"+distB+" dist:"+dist);
						
						if((distB+distC+1)> dist){
							//System.out.println("TEST3");
							//on modifie la coordonne y de la base actuelle sur le brin central
							double diff = (distB+distC+1) - dist;
							this.getListBasesMotif().get(indice).setCoords(new Point2D.Double
									(this.getListBasesMotif().get(indice).getCoords().x ,
											this.getListBasesMotif().get(indice).getCoords().y+diff));
							this.ajusteSpecificStrand(partenaire.getNumStrand());
						}
						
					}
				}
				previousStrand = partenaire.getNumStrand();
				this.getListStrand().getStrand(partenaire.getNumStrand()).setHasBeenPlaced(true);
			}
		}
	}
	
	public void setCenterMotif(){
		for (ModeleStrand ms: this.getListStrand().getListStrands())
		{//System.out.println(ms.getStrandLeft()+" "+ ms.getStrandRight());
			for (ModeleBase mb : ms.getArrayListMB())
			{
				int indice = mb.getIndex();
				ArrayList<ModeleBase> amb = _rna.getAllPartners(indice);
				
				boolean useBP = (amb.size()>0); 
				if (useBP)
				{
					ModeleBase partner = amb.get(0);
					if (!mb.getStyleBP().getStyle().isBent())
					{
						mb.setCenter(new Point2D.Double(
								(mb.getCoords().x+partner.getCoords().x)/2.0
								, mb.getCoords().y));
					}
					else
					{useBP = false;}
				}
				
				if (!useBP)
				{
					if(ms.getStrandLeft())
					{
						mb.setCenter(new Point2D.Double(
								mb.getCoords().x+10.0
								, mb.getCoords().y));
					}
					else if((!ms.getStrandRight() && !ms.getStrandLeft())|| ms.getStrandRight())
						
					{// Strandright == true OU Brin central 
						mb.setCenter(new Point2D.Double(
								mb.getCoords().x-10.0
								, mb.getCoords().y));
					}
				}
			}
		}
	}
	
	public void deviationBasePair(){
		ArrayList<ModeleBase> basesPair = listBasePair();
		for (int i=0; i<basesPair.size(); i++){
			int indice = basesPair.get(i).getIndex();
			ArrayList<ModeleBase> mb = _rna.getAllPartners(indice);
			int partner = mb.get(0).getIndex();
			double coordXBP = getListBasesMotif().get(indice).getCoords().x;
			double coordYBP = getListBasesMotif().get(indice).getCoords().y;
			double coordXBPpartner = getListBasesMotif().get(partner).getCoords().x;
			double coordYBPpartner = getListBasesMotif().get(partner).getCoords().y;
			//System.out.println(partner+" "+indice);
			for(int j=0; j<getListBasesMotif().size(); j++){
				
				double coordX = getListBasesMotif().get(j).getCoords().x;
				double coordY = getListBasesMotif().get(j).getCoords().y;
				if(j!=indice){
					if ((coordY==coordYBP && coordY==coordYBPpartner) || 
							(coordX==coordXBP && coordX==coordXBPpartner)){
						if(((coordX<coordXBP && coordX>coordXBPpartner)||(coordX<coordXBPpartner && coordX>coordXBP))
								|| (coordY<coordYBP && coordY>coordYBPpartner)||(coordY<coordYBPpartner && coordY>coordYBP)){
						
							if(!basesPair.get(i).getStyleBP().getStyle().isBent()){
								basesPair.get(i).getStyleBP().getStyle().setBent(1.0);
							}
						}
					}
				}
			}
		}
	}
	
	/**
	 * Create an arrayList which contain all bases which are involved in a base pair. 
	 */
	public ArrayList<ModeleBase> listBasePair(){
		ArrayList<ModeleBase> result = new ArrayList<ModeleBase>();
		for (int i = 0; i < getListStrand().sizeStruct(); i++ ){
			this.initChecker(getListStrand().getStrand(i));		
			for (int j = 0; j < getListStrand().getStrand(i).sizeStrand(); j++ ){
				int indice = this.getListStrand().getStrand(i).getMB(j).getIndex();
				int partner=-1;
				ArrayList<ModeleBase> mb = _rna.getAllPartners(indice);
			
				if (mb.size()>0){
					partner = mb.get(0).getIndex();
				}
				// Si la base est liée à une autre base
				if (partner!=-1 /*&& !getListBasesMotif().get(partner).getChecker()*/ ) {
					result.add(getListBasesMotif().get(indice));
					//getListBasesMotif().get(indice).setChecker(true);
				}
			}
		}
		
		return result;
	}
	
	public ModeleBase getBaseAt(int index)
	    {
	    	return getListBasesMotif().get(index);
	    }
	 
	public int getAjustement(){
		return this._ajustement;
	}
	
	public void setAjustement(int a){
		this._ajustement = a;
	}
	
	public double getDecalage(){
		return this._decalage;
	}
	
	public void setDecalage(double a){
		this._decalage = a;
	}
	
	public StructureTemp getListStrand(){
		return this._listStrands;
	}
	
	public ArrayList<ModeleBase> getListBasesMotif(){
		return this._listeBasesM;
	}
}

