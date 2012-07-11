package com.ppfold.main;

public class DataInfo {
	
	String fileName;
	String sequenceName;
	int sequenceID;
	String distFileName;
	int contactDistance;
	
	//Different kinds of data:
	//0: Bars
	//1: General mapping
	//2: Force as hard constraints
	//3: Contact distance (datadistFileName contains number)
	int type;
	
	String iD;

	public String getFileName() {
		return fileName;
	}

	public String getSequenceName() {
		return sequenceName;
	}

	public int getSequenceID() {
		return sequenceID;
	}

	public String getDistFileName() {
		return distFileName;
	}

	public int getContactDistance() {
		return contactDistance;
	}

	public int getType() {
		return type;
	}

	public void setFileName(String datafilename) {
		this.fileName = datafilename;
	}

	public void setSequenceName(String datasequencename) {
		this.sequenceName = datasequencename;
	}

	public void setSequenceID(int datasequenceID) {
		this.sequenceID = datasequenceID;
	}

	public void setDistFileName(String datadistfilename) {
		this.distFileName = datadistfilename;
	}

	public void setContactDistance(int contactDistance) {
		this.contactDistance = contactDistance;
	}

	public void setType(int datatype) {
		this.type = datatype;
	}

	public String getiD() {
		return iD;
	}

	public void setiD(String iD) {
		this.iD = iD;
	} 
	
	
	

}
