package com.ppfold.algo;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

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
}