package statalign.io;

import java.io.IOException;
import java.util.ArrayList;


/**
 * 
 * @author novak
 *
 */
public class RawSequences {

	/**
	 * Dynamic array of raw sequence data. Sequences can be aligned, in which case the
	 * strings must be of equal length and the '-' character must represent the gaps.
	 */
	public ArrayList<String> sequences;
	
	/**
	 * Dynamic array of sequence names. Any name can be null.
	 */
	public ArrayList<String> seqNames;
	
	/**
	 * Sorted string of characters present in sequences. Must not contain the gap '-'.
	 */
	public String alphabet;
	
	public RawSequences() {
		sequences = new ArrayList<String>();
		seqNames = new ArrayList<String>();
		alphabet = "";
	}
	
	public boolean isAligned() {
		int size;
		if((size = sequences.size()) == 0)
			return true;			// or should it be false?
		int len = sequences.get(0).length();
		for(int i = 1; i < size; i++)
			if(sequences.get(i).length() != len)
				return false;
		return true;
	}
	
	public void add(RawSequences more) {
		sequences.addAll(more.sequences);
		seqNames.addAll(more.seqNames);
		String alpha2 = more.alphabet;
		int len1 = alphabet.length(), len2 = alpha2.length();
		if(len1 == 0) {
			alphabet = alpha2;
			adjustNameLengths();
			return;
		}
		if(len2 == 0)
			return;
		StringBuilder merged = new StringBuilder(len1+len2);
		int i = 0, j = 0;
		char char1 = 0, char2 = 0;
		while(i < len1 || j < len2) {
			if(j == len2 || (i < len1 && (char1=alphabet.charAt(i)) < (char2=alpha2.charAt(j))))
				merged.append(alphabet.charAt(i++));
			else if(i == len1 || char1 > char2)
				merged.append(alpha2.charAt(j++));
			else {
				merged.append(char1);
				i++;
				j++;
			}
		}
		alphabet = merged.toString();
		adjustNameLengths();
	}

	public void adjustNameLengths(){
		int maxLength = 0;
		for(int i = 0; i < seqNames.size(); i++){
			maxLength = Math.max(maxLength, seqNames.get(i).length());
		}
		ArrayList<String> newNames = new ArrayList<String>();
		for(int i = 0; i < seqNames.size(); i++){
			String temp = seqNames.get(i);
			while(temp.length() < maxLength){
				temp += " ";
			}
			newNames.add(temp);
		}
		seqNames = newNames;	
	}
	
	/**
	 * Returns the number of sequences.
	 */
	public int size() {
		return seqNames.size();
	}
	
	/**
	 * Tells whether or not the the sequences are RNA/DNA. The extra letters denote ambiguous nucleotides.
	 * @return True if it is RNA, false otherwise
	 */
	public boolean isRNA() {
		for(int i = 0; i < alphabet.length(); i++) {
			char letter = alphabet.charAt(i);
			if(!(letter == 'A' || letter == 'C' || letter == 'G' 
				|| letter == 'U' || letter == 'T' || letter == 'W' 
				|| letter == 'S' || letter == 'R' || letter == 'Y'
				|| letter == 'K' || letter == 'M' || letter == 'D'
				|| letter == 'D' || letter == 'V' || letter == 'H'
				|| letter == 'B' || letter == 'X' || letter == 'N')) {
				return false;
			}
		}
		return true;
	}

	public static void main(String[] args) throws IOException {
		RawSequences r1 = new RawSequences();
		r1.alphabet = "cefhil";
		RawSequences r2 = new RawSequences();
		r2.alphabet = "aeghijko";
		r1.add(r2);
		System.out.println(r1.alphabet);
	}
}
