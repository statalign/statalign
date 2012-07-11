package statalign.model.subst.plugins;

import java.io.IOException;
import statalign.model.subst.plugins.NucleotideModel;
import statalign.model.subst.plugins.NucleotideModelUtils;
import statalign.base.*;

/**
 * This class implements the Tamura 1992 model.
 * 
 * @author lyngsoe
 */
public class Tamura92 extends NucleotideModel{

	public static final String menuName = "Tamura 1992";

	static final double span = 0.1; // largest parameter change when resampling

	/**
	 * This constructor uses numerical methods to diagonalise the rate
	 * matrix and to set diagonalising matrices.
	 * 
	 * @throws IOException
	 */
	public Tamura92() throws IOException{
		/* Initialise NucleotideModel part of object */
		super();
		/* Initialise Tamura 1992 specific part of object */
		params = new double[2];
		oldparams = new double[2];
		/* Set all frequencies to 1/4 */
		params[0] = 0.5;
		/* Set all rate ratios to 1 */
		params[1] = 1.0;
		/* Set diagonalising matrices and equilibrium frequencies */
		setDiagonal();
	}

	void setDiagonal(){
		/* Set equilibrium frequencies */
		e[0] = e[3] = (1 - params[0]) / 2.0;
		e[2] = e[1] = params[0] / 2.0;
		/* Set diagonalising matrices */
		// Set diagonal matrix
		d[0] = params[0] * params[1] - params[0] - params[1];
		d[1] = params[0] - params[0] * params[1] - 1;
		d[2] = 0;
		d[3] = (1 - params[0]) / (params[0] - 1);
		// Set v
		for (int i = 0; i < alphabet.length; i++){
			for (int j = 0; j < alphabet.length; j++)
				v[i][j] = 1;
			v[i][1 - i/2] = 0;
		}
		v[0][0] = v[2][1] = -1;
		v[0][3] = v[1][3] = params[0] / (params[0] - 1);
		// Set w
		for (int i = 0; i < alphabet.length; i++){
			w[i / 2][i] = (i % 2 == 0 ? -0.5 : 0.5);
			w[1 - i / 2][i] = 0;
			w[2 + i / 2][i] = (1 - params[0]) / 2.0;
			w[3 - i / 2][i] = (params[0] - 1 + i / 2) / 2.0;
		}
	}
	
	/**
	 * Prints the parameters of the model
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
		double x;

		if (Utils.generator.nextInt(2) == 0){
			/* Resample frequency parameter */
			oldparams[0] = params[0];
			boolean resample = true;
			double r, l;
			do{ /* Frequencies are not allowed to be zero or one */
				l = Math.max(0.0, oldparams[0] - span / 2.0);
				r = Math.min(1.0, oldparams[0] + span / 2.0);
				params[0] = l + Utils.generator.nextDouble() * (r - l);
				if ((params[0] > 0) && (params[0] < 1))
					resample = false;
			} while(resample);
			x = Math.log((Math.min(1.0, params[0] + span / 2.0) - Math.max(0, params[0] - span / 2.0)) / (r - l));
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
