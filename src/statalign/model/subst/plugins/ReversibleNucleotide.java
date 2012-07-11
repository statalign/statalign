package statalign.model.subst.plugins;

import java.io.IOException;

import statalign.base.Utils;

/**
 * This class implements the general time reversible parameter model.
 * 
 * @author lyngsoe
 */
public class ReversibleNucleotide extends NucleotideModel{

	public static final String menuName = "General Reversible (Nucl.)";

	static final double span = 0.1; // largest parameter change when resampling

	/**
	 * This constructor uses numerical methods to diagonalise the rate
	 * matrix and to set diagonalising matrices.
	 * @throws IOException
	 */
	public ReversibleNucleotide() throws IOException{
		/* Initialise NucleotideModel part of object */
		super();
		/* Initialise Reversible specific part of object */
		params = new double[8];
		oldparams = new double[8];
		/* Set all frequencies to 1/4 */
		params[0] = params[1] = params[2] = 0.25;
		/* Set all rate ratios to 1 */
		params[3] = params[4] = params[5] = params[6] = params[7] = 1.0;
		/* Set diagonalising matrices and equilibrium frequencies */
		setDiagonal();
	}

	void setDiagonal(){
		/* Set equilibrium frequencies */
		e[alphabet.length - 1] = 1.0;
		for (int i = 0; i < alphabet.length - 1; i++){
			e[i] = params[i];
			e[alphabet.length - 1] -= params[i];
		}
		/* Set diagonalising matrices */
		// Start by constructing rate matrix
		double Q[][] = new double[alphabet.length][alphabet.length];
		for (int i = 0; i < alphabet.length; i++)
			Q[i][i] = 0.0;
		for (int i = 1; i < alphabet.length; i++)
			for (int j = 0; j < i; j++){
				if (j < 2){
					Q[i][j] = e[j] * params[2 + 2 * j + i];
					Q[j][i] = e[i] * params[2 + 2 * j + i];
				}
				else{
					Q[i][j] = e[j];
					Q[j][i] = e[i];
				}
				Q[i][i] -= Q[i][j];
				Q[j][j] -= Q[j][i];
			}
		NucleotideModelUtils.numericalDiagonalisation(this, Q);
	}

	/**
	 * Prints the parameters of the model
	 */
	@Override
	public String print(){
		String s = "";
		for(int i = 3; i < params.length; i++){
			s+="\trho_"+(i-2)+"\t"+params[i];
		}
		return "A\t"+e[0]+"\tC\t"+e[1]+"\tG\t"+e[2]+"\tT\t"+e[3]+s;
	}

	/**
	 * This function implements a proposal for new parameter values.
	 * Returns with the logarithm of the Metropolis-Hastings ratio.
	 */
	@Override
	public double sampleParameter() {
		double x;

		if (Utils.generator.nextInt(8) < 3){
			/* Resample frequency parameter */
			x = NucleotideModelUtils.sampleFrequencies(e, span / 2.0);
			/* Save current parameters as old parameters, store new
			 * frequency parameters, and ensure frequencies sum to one.
			 */
			e[alphabet.length - 1] = 1.0;
			for (int i = 0; i < params.length; i++){
				oldparams[i] = params[i];
				if (i < alphabet.length - 1){
					params[i] = e[i];
					e[alphabet.length - 1] -= e[i];
				}
			}
		}
		else{
			/* Resample rate ratio parameter */
			double r[] = new double[5];
			for (int i = 0; i < 5; i++)
				r[i] = params[3 + i];
			x = NucleotideModelUtils.sampleRates(r, span / 2.0);
			/* save all params (bugfix), we thank William Majoros for reporting */
            for (int i = 0; i < params.length; i++)
            	oldparams[i] = params[i];
            for (int i = 0; i < 5; i++)
				params[3 + i] = r[i];

		}
		setDiagonal();
		return x;
	}
}
