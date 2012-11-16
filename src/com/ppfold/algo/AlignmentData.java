package com.ppfold.algo;

import java.io.Serializable;
import java.util.ArrayList;

public class AlignmentData implements Serializable
{
	public AlignmentData(String[] seq){
		for(int i = 0; i<seq.length; ++i){
			sequences.add(seq[i]);
			names.add(Integer.toString(i));
		}
		
		
	}
	public AlignmentData(){
		
	}
	
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -62635487259521461L;
	public ArrayList<String> sequences = new ArrayList<String>();
	public ArrayList<String> names = new ArrayList<String>();
	
	public String toString()
	{
		String ret = "";
		for(int i = 0 ; i < names.size() ; i++)
		{
			ret += ">" + names.get(i)+"\n";
			ret += sequences.get(i)+"\n";
		}
		return ret;
	}
}