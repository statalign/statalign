package statalign.model.subst.plugins;

import java.io.IOException;

/**
 * Implements the Jukes-Cantor model for nucleic acids.
 *
 * @author miklos, novak
 *
 */
public class JukesCantor extends NucleotideModel {

	public static final String menuName = "Jukes-Cantor";

	static final double span = 0.1;

	/**
	 * This constructor reads transition rates from the file data/jukescantor_rate.dat,
	 * the alphabet from data/DNAalphabet.dat, and the equilibrium distribution from 
	 * data/jukescantor_equilibrium.dat.
	 * 
	 * @throws IOException
	 */
	public JukesCantor() throws IOException {
		super("data/jukescantor_rate.dat", "data/jukescantor_equilibrium.dat");
		params = new double[0];
		d[0] = 0.0;
		d[1] = -4.0/3.0;
		d[2] = d[1];
		d[3] = d[1];
	}

	/**
	 * Empty function (no model parameters)
	 */
	@Override
	public void restoreParameter() {
	}

	/**
	 * It does nothing, and always return with 0, namely, log-probability 1.
	 */
	@Override
	public double sampleParameter() {
		return 0.0;
	}

}
