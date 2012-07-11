package com.ppfold.main;

import java.util.ArrayList;
import java.util.List;

import com.ppfold.algo.MatrixTools;

public class Alignment {

	private List<String> sequences = new ArrayList<String>();
	private List<String> names = new ArrayList<String>();
	
	public Alignment(List<String> seqs, List<String> nams){
		sequences = seqs;
		names = nams; 
	}
	
	public List<String> getSequences() {
		return sequences;
	}
	public List<String> getNames() {
		return names;
	}
	
	public void print(){
		for(int i = 0; i < sequences.size(); i++){
			System.out.println(">" + names.get(i));
			System.out.println(sequences.get(i));
		}
	}

	public int calculateLength(int id) {
		String sequence = sequences.get(id);
		int length = 0;
		for(int i = 0; i<sequence.length(); i++){
			char c = sequence.charAt(i);
			if(!MatrixTools.isGap(c)){
				length++;
			}
		}
		return length;
	}
	

	
}
