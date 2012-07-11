package statalign.model.subst.plugins;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import statalign.io.RawSequences;
import statalign.model.score.plugins.DNAScore;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;

/**
 * This is the abstract class for nucleotide models.
 * 
 * @author lyngsoe
 */
public abstract class NucleotideModel extends SubstitutionModel{

	/* Old model parameters */
	protected double oldparams[];
	
	public static String type = "nucleotide";
	
	private boolean printRNA;		// true if should print RNA nucleotides (U for T)

	/**
	 * This constructor reads the alphabet from data/DNAalphabet.dat
	 * @throws IOException
	 */
	protected NucleotideModel() throws IOException{
		attachedScoringScheme = new DNAScore();
		String alphabetFile = "data/DNAalphabet.dat";
		ClassLoader cl = getClass().getClassLoader();
		BufferedReader bf = new BufferedReader(new InputStreamReader(cl.getResource(alphabetFile).openStream()));
		String a = bf.readLine();
		alphabet = new char[a.length()];
		for(int i = 0; i < alphabet.length; i++){
			alphabet[i] = a.charAt(i);
		}
		int size = alphabet.length;
		v = new double[size][size];
		w = new double[size][size];
		d = new double[size];
		e = new double[size];
	}

	/**
	 * Reads rates and equilibriums from given files. Does not set the params array.
	 * 
	 * @param rateFile substitution rates file
	 * @param equFile equilibrium probabilities file
	 */
	protected NucleotideModel(String rateFile, String equFile) throws IOException {
		attachedScoringScheme = new DNAScore();
		String alphabetFile = "data/DNAalphabet.dat";
		ClassLoader cl = getClass().getClassLoader();
		//System.out.println(cl+" "+cl.getResource(alphabetFile)+" "+alphabetFile);
		BufferedReader bf = new BufferedReader(new InputStreamReader(cl.getResource(alphabetFile).openStream()));
//		BufferedReader bf = new BufferedReader(new FileReader(alphabetFile));
		String a = bf.readLine();
		alphabet = new char[a.length()];
		for(int i = 0; i < alphabet.length; i++){
			alphabet[i] = a.charAt(i);
		}
		int size = alphabet.length;
		v = new double[size][size];
		w = new double[size][size];
		d = new double[size];
		e = new double[size];
		bf = new BufferedReader(new InputStreamReader(cl.getResource(rateFile).openStream()));
		//        bf = new BufferedReader(new FileReader(rate));
		String[] temp;
		for(int i = 0; i < size; i++){
			temp = (bf.readLine()).split(" ");
			for(int j = 0; j < size; j++){
				v[i][j] = Double.parseDouble(temp[j]);
			}
		}
		String t = bf.readLine(); //reading an empty line
		for(int i = 0; i < size; i++){
			temp = (bf.readLine()).split(" ");
			for(int j = 0; j < size; j++){
				//    System.out.println(i+" "+j+" "+temp[j]);
				w[i][j] = Double.parseDouble(temp[j]);
			}
		}
		t = bf.readLine(); //reading an empty line
		bf = new BufferedReader(new InputStreamReader(cl.getResource(equFile).openStream()));
		//    bf = new BufferedReader(new FileReader(equilibrium));
		for(int i = 0; i < size; i++){
			t = bf.readLine();
			e[i] = Double.parseDouble(t);
		}
	}

	/**
	 * This function decides if this model can analyze input sequences.
	 * accepted characters are 'acgutryswmkbdhvn'. Our program is
	 * case insensitive.
	 */
	@Override
	public double acceptable(RawSequences r) {
		String acceptableCharacters = "acgturyswmkbdhvn";

		int[] count = new int[5];	// nucleotide counts
		for(int i = 0; i < r.sequences.size(); i++){
			String sequence = r.sequences.get(i);
			for(int j = 0; j < sequence.length(); j++){
				char ch = sequence.charAt(j);
				if(ch != '-' && ch != ' ') {
					ch = Character.toLowerCase(ch);
					int ind = acceptableCharacters.indexOf(ch);
					if(ind == -1){
						throw new RecognitionError(getMenuName()+" cannot accept the sequences because it contains character '"+ch+"'!\n");
					} else if(ind < count.length){
						count[ind]++;
					}

				}
			}
		}
		
		int uind = count.length-1;
		if(count[uind] > count[uind-1])
			printRNA = true;		// automatically choose RNA output if seen more U than T
		else if(count[uind-1] > count[uind])
			printRNA = false;
		
		count[uind-1] += count[uind];	// then treat Us as Ts
		int sum = 0;
		for(int i = 0; i < uind; i++){
			if(count[i] > 0)
				sum++;
		}
		return (double)sum/uind;
	}

	/**
	 * This function assigns colors to characters. 'A' is red, 'C' is blue,
	 * 'G' is orange, 'T' or 'U' is green, the remaining characters (ambiguous characters
	 * and the gap symbol) are grey 
	 */
	@Override
	public Color getColor(char c) {
		c = Character.toUpperCase(c);
		switch (c) {
		case 'A': return Color.RED;
		case 'C': return Color.BLUE;
		case 'G': return Color.ORANGE;
		case 'T':
		case 'U': return Color.GREEN;
		default: return Color.GRAY;
		}
	}

	/**
	 * returns the description of the current parameters of the model
	 */
	@Override
	public String print() {
		return "";
	}

	void setDiagonal(){}
	
	/**
	 * restore the parameters to the old values when a parameter-changing
	 * proposal is not accepted.
	 */
	@Override
	public void restoreParameter() {
		for(int i = 0; i < params.length; i++){
			params[i] = oldparams[i];
		}
		setDiagonal();
	}

	/**
	 * Returns with the most likely character, given a Felsentein
	 * likelihood array.  It can handle ambiguous characters. If the
	 * probabilities in the Felsentein array are all zero, returns '*'.
	 */
	@Override
	public char mostLikely(double[] seq) {
		/* Create binary encoding of characters with probability
		 * higher than 0.5.
		 */
		int code = 0;
		int radix = 1;
		double max = 0.0;
		char character = '*';
		for(int i = 0; i < seq.length; i++){
			if(seq[i] > 0.5){
				/* Character has probability higher than 0.5 */
				if (max > 0.5)
					/* ...but is not the first high probability
					 * character seen.
					 */
					code += radix;
				else{
					/* ...and is the first high probability character seen */
					max = seq[i];
					code = radix;
				}
			} else{
				/* Character has probability at most 0.5 */
				if (seq[i] > max){
					/* ...but is most probable character seen so far */
					max = seq[i];
					code = radix;
				} else if (seq[i] == max)
					/* ...but is among the most probable characters
					 * seen so far.
					 */
					code += radix;
			}
			radix <<= 1;
		}
		if (max == 0.0)
			code = 0;
		/* Look up corresponding character */
		String conversionTable = printRNA ? "*ACMGRSVUWYHKDBN" : "*ACMGRSVTWYHKDBN";
		character = conversionTable.charAt(code);
		return character;
	}
	
	/**
	 * Selects whether to print RNA nucleotides (U for T) or default DNA ones.
	 * Affects the beviour of {@link #mostLikely(double[])}.
	 */
	public void setPrintRNA(boolean printRNA) {
		this.printRNA = printRNA;
	}
}

