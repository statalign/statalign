package statalign.model.score.plugins;

import java.io.FileNotFoundException;
import java.io.IOException;

import statalign.io.input.FileTokenReader;
import statalign.model.score.SubstitutionScore;

/**
 * This is an implementation of a DNA scoring.
 * 
 * @author miklos
 *
 */
public class DNAScore extends SubstitutionScore {
	
	final String dnaCodes = "ACGT";

	/**
	 * This constructor reads the distance matrix from data/dnadist.dat,
	 * and sets the which[][] array.
	 * 
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public DNAScore() throws IOException, FileNotFoundException{
		which = new int[256][4];   /* tells amino acid no. for each valid char */
		dist = new int[4][4];  /* amino acid distances in Gotoh's alg. */

		//	System.out.println("I'm in Blosum62 constructor!");
		//	System.out.println("I'm in Blosum62 constructor! "+which);
		//System.out.println("I'm in Blosum62 constructor! "+which.length);

		FileTokenReader file = new FileTokenReader("data/dnadist.dat");
		for(int i = 0; i < 4; i++){
			for(int j = 0; j <= i; j++){
				dist[i][j] = dist[j][i] = file.readDbl().intValue();
			}
		}
		file.close();

		//for(int i = 0; i < 256; i++)
		//	which[i] = -1;  /* no amino acid associated */
		for(int i = 0; i < dnaCodes.length(); i++){
			which[dnaCodes.charAt(i)][i] = which[Character.toLowerCase(dnaCodes.charAt(i))][i] = 1;
		}
		//IUPAC codes, a, c, g, t
		// puRin: a, g
		which['r'][0] = which['R'][0] = 1;
		which['r'][2] = which['R'][2] = 1;
		// pYrimidin: c, t
		which['y'][1] = which['Y'][1] = 1;
		which['y'][3] = which['Y'][3] = 1;
		// Strong: c,g
		which['s'][1] = which['S'][1] = 1;
		which['s'][2] = which['S'][2] = 1;
		// Weak: a,t
		which['w'][1] = which['W'][1] = 1;
		which['w'][3] = which['W'][3] = 1;
		// Keto: 
		which['k'][2] = which['K'][2] = 1;
		which['k'][3] = which['K'][3] = 1;
		// aMino:  
		which['m'][0] = which['M'][0] = 1;
		which['m'][1] = which['M'][1] = 1;
		// B: not-A, c, g, t
		which['b'][1] = which['B'][1] = 1;
		which['b'][2] = which['B'][2] = 1;
		which['b'][3] = which['B'][3] = 1;
		// D: not-C, a,g,t
		which['d'][0] = which['D'][0] = 1;
		which['d'][2] = which['D'][2] = 1;
		which['d'][3] = which['D'][3] = 1;
		// H: not-G, a,c,t
		which['h'][0] = which['H'][0] = 1;
		which['h'][1] = which['H'][1] = 1;
		which['h'][3] = which['H'][3] = 1;
		// V: not-U, a, c, g
		which['v'][0] = which['V'][0] = 1;
		which['v'][1] = which['V'][1] = 1;
		which['v'][2] = which['V'][2] = 1;
		// aNy: a,c,g,t
		which['n'][0] = which['N'][0] = 1;
		which['n'][1] = which['N'][1] = 1;
		which['n'][2] = which['N'][2] = 1;
		which['n'][3] = which['N'][3] = 1;

		which['u'] = which['U']= which['t'];

	}

}
