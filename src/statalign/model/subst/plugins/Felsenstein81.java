package statalign.model.subst.plugins;

import java.io.IOException;
import statalign.model.subst.plugins.NucleotideModel;
import statalign.model.subst.plugins.NucleotideModelUtils;

/**
 * This class implements the Felsenstein 1981 parameter model.
 * 
 * @author lyngsoe
 */
public class Felsenstein81 extends NucleotideModel{
	
	public static final String menuName = "Felsenstein 1981";

	static final double span = 0.1; // largest parameter change when resampling

	/**
	 * This constructor uses the analytical solution to diagonalising
	 * the rate matrix to set diagonalising matrices.
	 * 
	 * @throws IOException
	 */
	public Felsenstein81() throws IOException{
		/* Initialise NucleotideModel part of object */
		super();
		/* Initialise Felsenstein81 specific part of object */
		params = new double[3];
		oldparams = new double[3];
		/* Set all frequencies to 1/4 */
		params[0] = params[1] = params[2] = 0.25;
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
		for (int i = 0; i < alphabet.length; i++)
			for (int j = 0; j < alphabet.length; j++){
				if (j == 0)
					v[i][j] = 1;
				else if (i == 0)
					v[i][j] = -e[alphabet.length - j] / e[0];
				else if (i + j == 4)
					v[i][j] = 1;
				else
					v[i][j] = 0;
				if (i == 0)
					w[i][j] = e[j];
				else if (i + j == 4)
					w[i][j] = 1 - e[j];
				else
					w[i][j] = -e[j];
			}
		/* Set diagonal matrix */
		d[0] = 0;
		d[1] = d[2] = d[3] = -1;
	}
	
	/**
	 * Prints the current state of model.
	 */
	public String print(){
		return "A\t"+e[0]+"\tC\t"+e[1]+"\tG\t"+e[2]+"\tT\t"+e[3];
	}

	/**
	 * This function implements a proposal for new parameter values.
	 * Returns with the logarithm of the Metropolis-Hastings ratio.
	 */
	@Override
	public double sampleParameter() {
		double x = NucleotideModelUtils.sampleFrequencies(e, span);
		/* Save current parameters as old parameters, store new
		 * frequency parameters, and ensure frequencies sum to one.
		 */
		e[params.length] = 1.0;
		for (int i = 0; i < params.length; i++){
			oldparams[i] = params[i];
			params[i] = e[i];
			e[params.length] -= e[i];
		}

		setDiagonal();

		return x;
	}
}
