package fr.orsay.lri.varna.models.rna;

public class ExtendedMB extends ModeleBase {

	private boolean intervDroite = false;
	
	private boolean intervGauche = false;
	
	
	public ExtendedMB(){
		
	}
	
	public boolean getIntervGaucheEx(){
		return (this.intervGauche);
	}
	
	public void setIntervGaucheEx(boolean bool){
		this.intervGauche=bool;
	}
	
	public boolean getIntervDroiteEx(){
		return (this.intervDroite);
	}
	
	public void setIntervDroiteEx(boolean bool){
		this.intervDroite=bool;
	}
	
	@Override
	public String getContent() {
		// TODO Auto-generated method stub
		return "";
	}

	@Override
	public int getIndex() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setContent(String s) {
		// TODO Auto-generated method stub
		
	}

}
