package statalign.model.subst.plugins;

import java.io.IOException;

import statalign.base.Utils;

/**
 * Implements the Kimura 3 parameter model.
 * 
 * @author miklos
 */
public class Kimura3 extends NucleotideModel {

	public static final String menuName = "Kimura 3 parameters";
	
	static final double span = 0.1;

	/**
	 * This constructor reads transition rates from the file data/kimura3_rate.dat,
	 * the alphabet from data/DNAalphabet.dat, and the equilibrium distribution from 
	 * data/kimura3_equilibrium.dat.
	 * 
	 * @throws IOException
	 */
	public Kimura3() throws IOException {
		super("data/kimura3_rate.dat", "data/kimura3_equilibrium.dat");
		double alpha = 0.5;
		double beta  = 0.5;
		double gamma = 0.5;
		params = new double[3];
		oldparams = new double[3];
		params[0] = alpha; params[1] = beta; params[2] = gamma;
		setDiagonal();
	}

	@Override
	void setDiagonal(){
		d[0] = 0;
		d[1]= -2*params[0]-2*params[1];
		d[2] = -2*params[1]-2*params[2];
		d[3] = -2*params[0]-2*params[2]; 

	}

	/**
	 *  Prints the parameters of the model
	 */
	@Override
	public String print(){
		return "\talpha\t"+params[0]+"\tbeta\t"+params[1]+"\tgamma\t"+params[2];
	}

	/**
	 * This function implements a proposal for new parameter values.
	 * Returns with the logarithm of the Metropolis-Hastings ratio.
	 */
	@Override
	public double sampleParameter() {
		for(int i = 0; i < 3; i++){
			oldparams[i] = params[i];
		}
		int k = Utils.generator.nextInt(3);
		double delta = 0.0;
		do{
			delta = Utils.generator.nextDouble()*span - (span/2.0);
		} while(params[k] + delta <= 0.0);
		params[k] += delta;
		setDiagonal();
		return Math.log((params[k] > span/2.0 ? span : span/2.0 + params[k])/
				(oldparams[k] > span/2.0 ? span : span/2.0 + oldparams[k]));
	}

}
