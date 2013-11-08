package statalign.model.subst.plugins;

import java.io.IOException;
import statalign.model.subst.plugins.NucleotideModel;
import statalign.model.subst.plugins.NucleotideModelUtils;
import statalign.base.*;

/**
 * This class implements the Felsenstein 1984 parameter model.
 * 
 * @author lyngsoe
 */
public class Felsenstein84 extends NucleotideModel{

	public static final String menuName = "Felsenstein 1984";
	
	//static final double span = 0.1; // largest parameter change when resampling

	/**
	 * This constructor uses numerical methods to diagonalise the rate
	 * matrix and to set diagonalising matrices.
	 * @throws IOException
	 */
	public Felsenstein84() throws IOException{
		/* Initialise NucleotideModel part of object */
		super();
		/* Initialise Felsenstein 1984 specific part of object */
		params = new double[alphabet.length];
		oldparams = new double[alphabet.length];
		/* Set all frequencies to 1/4 */
		for (int i = 0; i < alphabet.length - 1; i++)
			params[i] = 1.0 / (alphabet.length);
		/* Set all rate ratios to 1 */
		params[alphabet.length - 1] = 1.0;
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
		/* Set diagonal matrix */
		d[0] = 0;
		d[1] = -1;
		d[2] = d[3] = -1 - params[params.length - 1];
		/* Set v matrix */
		for (int i = 0; i < alphabet.length; i++){
			v[i][0] = 1;
			if (i < 2){
				v[i][1] = -(e[2] + e[3]) / (e[0] + e[1]);
				v[i][2] = (i == 0 ? -e[1]/e[0] : 1);
				v[i][3] = 0;
			}
			else{
				v[i][1] = 1;
				v[i][2] = 0;
				v[i][3] = (i == 2 ? -e[3]/e[2] : 1);
			}
		}
		/* Set w matrix */
		for (int i = 0; i < alphabet.length; i++){
			w[0][i] = e[i];
			if (i < 2){
				w[1][i] = -e[i];
				w[2][i] = (i == 0 ? -1 : 1) * e[0] / (e[0] + e[1]);
				w[3][i] = 0;
			}
			else{
				w[1][i] = e[i] * (e[0] + e[1]) / (e[2] + e[3]);
				w[2][i] = 0;
				w[3][i] = (i == 2 ? -1 : 1) * e[2] / (e[2] + e[3]);
			}
		}
	}
	
	/**
	 * Prints the parameters of the model
	 */
	public String print(){
		return "A\t"+e[0]+"\tC\t"+e[1]+"\tG\t"+e[2]+"\tT\t"+e[3]+"\tTs/Tv\t"+params[params.length - 1];
	}

	/**
	 * This function implements a proposal for new parameter values.
	 * Returns with the logarithm of the Metropolis-Hastings ratio.
	 */
	@Override
	public double sampleParameter() {
		double x;

		if (Utils.generator.nextInt(alphabet.length) != 0){
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
			double r[] = new double[1];
			r[0] = params[params.length - 1];
			x = NucleotideModelUtils.sampleRates(r, span / 2.0);
			oldparams[params.length - 1] = params[params.length - 1];
			params[params.length - 1] = r[0];
		}
		setDiagonal();
		return x;
	}
}
