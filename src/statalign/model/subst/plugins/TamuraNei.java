package statalign.model.subst.plugins;

import java.io.IOException;

import statalign.base.Utils;

/**
 * This class implements the Tamura-Nei parameter model.
 * 
 * @author lyngsoe
 */
public class TamuraNei extends NucleotideModel{

	public static final String menuName = "Tamura/Nei";

	//static final double span = 0.1; // largest parameter change when resampling

	/**
	 * This constructor uses numerical methods to diagonalise the rate
	 * matrix and to set diagonalising matrices.
	 * 
	 * @throws IOException
	 */
	public TamuraNei() throws IOException{
		/* Initialise NucleotideModel part of object */
		super();
		/* Initialise Tamura-Nei specific part of object */
		params = new double[alphabet.length + 1];
		oldparams = new double[alphabet.length + 1];
		/* Set all frequencies to 1/4 */
		for (int i = 0; i < alphabet.length - 1; i++)
			params[i] = 1.0 / (alphabet.length);
		/* Set all rate ratios to 1 */
		params[alphabet.length - 1] = params[alphabet.length] = 1.0;
		/* Set diagonalising matrices and equilibrium frequencies */
		setDiagonal();
	}

	@Override
	void setDiagonal(){
		/* Set equilibrium frequencies */
		e[alphabet.length - 1] = 1.0;
		for (int i = 0; i < alphabet.length - 1; i++){
			e[i] = params[i];
			e[alphabet.length - 1] -= params[i];
		}
		/* Set diagonalising matrices */
		// Set d
		d[0] = -params[3] * (e[0] + e[1]) - e[2] - e[3];
		d[1] = -params[4] * (e[2] + e[3]) - e[0] - e[1];
		d[2] = 0;
		d[3] = -1;
		// Set v
		for (int i = 0; i < alphabet.length; i++){
			for (int j = 0; j < alphabet.length; j++)
				v[i][j] = 1;
			v[i][1 - i / 2] = 0;
		}
		v[0][3] = v[1][3] = -(e[2] + e[3]) / (e[0] + e[1]);
		v[0][0] = -e[1] / e[0];
		v[2][1] = -e[3] / e[2];
		// Set w
		for (int j = 0; j < alphabet.length; j++){
			w[j / 2][j] = (j % 2 == 0 ? -1 : 1) * e[(j < 2 ? 0 : 2)] / (e[j < 2 ? 0 : 2] + e[j < 2 ? 1 : 3]);
			w[1 - j / 2][j] = 0;
			w[2][j] = e[j];
			if (j < 2)
				w[3][j] = -e[j];
			else
				w[3][j] = e[j] * (e[0] + e[1]) / (e[2] + e[3]);
		}
	}
	
	/**
	 * Prints the parameters of the model
	 */
	@Override
	public String print(){
		String s = "";
		for(int i = 3; i < params.length; i++){
			s+="\tkappa_"+(i-2)+"\t"+params[i];
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

		if (Utils.generator.nextInt(alphabet.length + 1) > 1){
			/* Resample frequency parameter */
			x = NucleotideModelUtils.sampleFrequencies(e, span / 2.0);
			e[alphabet.length - 1] = 1.0;
			for (int i = 0; i < alphabet.length - 1; i++)
				e[alphabet.length - 1] -= e[i];
			/* Save current parameters as old parameters, store new
			 * frequency parameters, and ensure frequencies sum to one.
			 */
			for (int i = 0; i < params.length; i++){
				oldparams[i] = params[i];
				if (i < alphabet.length - 1)
					params[i] = e[i];
			}
		}
		else{
			/* Resample rate ratio parameter */
			double r[] = new double[2];
			for (int i = 0; i < 2; i++)
				r[i] = params[params.length - 2 + i];
			x = NucleotideModelUtils.sampleRates(r, span / 2.0);
			for (int i = 0; i < params.length; i++)
				oldparams[i] = params[i];
			for (int i = 0; i < 2; i++)
				params[params.length - 2 + i] = r[i];
		}
		setDiagonal();
		return x;
	}
}
