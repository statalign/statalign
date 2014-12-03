package statalign.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Dynamic array of raw sequence data. Sequences can be aligned, in which
 * case the strings must be of equal length and the '-' character must
 * represent the gaps.
 * 
 * @author novak, miklos, aszabo
 * 
 */
public class RawSequences implements DataType {

	/**
	 * Dynamic array of raw sequence data.
	 */
	private List<String> sequences = new ArrayList<String>();

	public RawSequences getSeqs() { return this; }
	public void setSeqs(RawSequences rs) { }
	
	/**
	 * Dynamic array of sequence names. Any name can be null.
	 */
	private List<String> seqNames = new ArrayList<String>();
	public List<String> getNames() { return seqNames; }

	/**
	 * Sorted string of characters present in sequences. Does not contain the
	 * gap '-'.
	 */
	public String alphabet;
	
    private int maxNameLength = -1;
	
	public boolean perSequenceData() {
		return false;
	}
	public String getSummaryAssociatedWith(String sequenceName) { 
		return sequenceName;
	}
	public void removeDataAssociatedWith(String sequenceName) { 
		
	}
	
	public RawSequences() {
		sequences = new ArrayList<String>();
		seqNames = new ArrayList<String>();		
	}
	
	public RawSequences(String seq, String name) {
		sequences = new ArrayList<String>();
		sequences.add(seq);
		seqNames = new ArrayList<String>();
		seqNames.add(name);
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

	/**
	 * Returns the common length of the aligned sequences or -1 if sequences are
	 * unaligned (have different length).
	 */
	public int len() {
		int size;
		if ((size = sequences.size()) == 0)
			return 0;
		int len = sequences.get(0).length();
		for (int i = 1; i < size; i++)
			if (sequences.get(i).length() != len)
				return -1;
		return len;
	}

	/**
	 * Add a name and a sequence, that belong together. If another sequence with
	 * the same name is already added then it will be updated.
	 */
	public void add(String name, String sequence) {
//		if (name != null && seqNames.contains(name))
//			throw new Error("Name collision! (" + name
//					+ "), please edit your input files to resolve it.");
		if (seqNames.contains(name)) {
			sequences.set(sequences.indexOf(name),sequence);						
		}
		else {
			seqNames.add(name);
			sequences.add(sequence);
		}
	}

	public void addOrReplace(String name, String sequence) {
		if (seqNames.contains(name)) {
			removeByName(name);
		}
		add(name, sequence);
	}

	public void add(RawSequences more) {
		// so that names are checked
		for (int i = 0; i < more.seqNames.size(); i++) {
			add(more.seqNames.get(i), more.sequences.get(i)); // addOrReplace?
		}

		String alpha1 = getAlphabet();
		String alpha2 = more.getAlphabet();
		int len1 = alpha1.length(), len2 = alpha2.length();
		if (len1 == 0) {
			alphabet = alpha2;
			maxNameLength = more.maxNameLength;
			return;
		}
		if (len2 == 0)
			return;
		StringBuilder merged = new StringBuilder(len1 + len2);
		int i = 0, j = 0;
		char char1 = 0, char2 = 0;
		while (i < len1 || j < len2) {
			if (j == len2 || (i < len1 && (char1 = alpha1.charAt(i)) < (char2 = alpha2.charAt(j))))
				merged.append(alpha1.charAt(i++));
			else if (i == len1 || char1 > char2)
				merged.append(alpha2.charAt(j++));
			else {
				merged.append(char1);
				i++;
				j++;
			}
		}
		alphabet = merged.toString();
		maxNameLength = -1;
	}

	public int getMaxNameLength() {
		if (maxNameLength == -1) {
			maxNameLength = 0;
			for (String name : seqNames) {
				if (name != null && name.length() > maxNameLength) {
					maxNameLength = name.length();
				}
			}
		}
		return maxNameLength;
	}

	public String getAlphabet() {
		if (alphabet == null) {
			boolean present[] = new boolean[256];
			for (String seq : sequences) {
				for (int i = 0; i < seq.length(); i++) {
					char ch = seq.charAt(i);
					if (ch < 256)
						present[Character.toUpperCase(ch)] = true;
				}
			}
			StringBuilder aBuilder = new StringBuilder();
			for (char ch = 'A'; ch <= 'Z'; ch++) {
				if (present[ch])
					aBuilder.append(ch);
			}
			alphabet = aBuilder.toString();
		}
		return alphabet;
	}

//	public void adjustNameLengths(){
//		int maxLength = 0;
//		for(int i = 0; i < seqNames.size(); i++){
//			maxLength = Math.max(maxLength, seqNames.get(i).length());
//		}
//		ArrayList<String> newNames = new ArrayList<String>();
//		for(int i = 0; i < seqNames.size(); i++){
//			String temp = seqNames.get(i);
//			while(temp.length() < maxLength){
//				temp += " ";
//			}
//			newNames.add(temp);
//		}
//		seqNames = newNames;	
//	}

	/**
	 * Returns the number of sequences.
	 */
	public int size() {
		return seqNames.size();
	}

	public String getSeqName(int i) {
		return seqNames.get(i);
	}
	
	public String getSequence(int i) {
		return sequences.get(i);
	}
	
	public List<String> getSeqnames() {
		return Collections.unmodifiableList(seqNames);
	}
	
	public List<String> getSequences() {
		return Collections.unmodifiableList(sequences);
	}
	
	public String getSeqNamePadded(int i) {
		String n = seqNames.get(i);
		int pad = getMaxNameLength()-n.length();
		if(pad > 0) {
			char[] arr = new char[pad];
			Arrays.fill(arr, ' ');
			n += new String(arr);
		}
		return n;
	}

	/**
	 * Removes a sequence and its name.
	 * 
	 * @return true on success
	 */
	public boolean remove(int i) {
		if (size() > i && i >= 0) {
			sequences.remove(i);
			seqNames.remove(i);			
			maxNameLength = -1;
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Removes a sequence and its name (given by name).
	 * 
	 * @return true on success
	 */
	public boolean removeByName(String name) {
		int ind = seqNames.indexOf(name);
		if (ind == -1)
			return false;
		remove(ind);
		return true;
	}

	/**
	 * Removes all sequences (and their names).
	 */
	public void clear() {
		seqNames.clear();
		sequences.clear();
		alphabet = null;
		maxNameLength = -1;
	}

	/**
	 * Removes all gaps from each sequence.
	 */
	public void removeGaps() {
		for (int i = 0; i < sequences.size(); i++) {
			String seq = sequences.get(i);
			StringBuilder builder = new StringBuilder(seq.length());
			for (int j = 0; j < seq.length(); j++) {
				char ch = seq.charAt(j);
				if (ch != '-')
					builder.append(ch);
			}
			sequences.set(i, builder.toString());
		}
	}

	@Override
	public String toString() {
		String s = "";
		for (int i = 0; i < seqNames.size(); i++) {
			s += ">" + seqNames.get(i) + "\n" + sequences.get(i) + "\n";
		}
		return s;
	}

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
		return (alphabet.indexOf("U") != -1);
	}
	
	public static void main(String[] args) throws IOException {
		RawSequences r1 = new RawSequences();
		r1.sequences.add("LceIfhFil");
		RawSequences r2 = new RawSequences();
		r2.sequences.add("aeghijko");
		r1.add(r2);
		System.out.println(r1.alphabet);
	}
}