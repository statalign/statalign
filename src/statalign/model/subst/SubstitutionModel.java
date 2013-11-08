package statalign.model.subst;

import java.awt.*;
import java.util.Arrays;

import statalign.model.score.*;

import statalign.io.*;

/**
 * 
 * Superclass of the substitution models.
 * 
 * @author miklos, novak
 *
 */
public abstract class SubstitutionModel{

	/**
	 * The rate matrix must be diagonalized and represented in a
	 * <b>v d w</b> product, where <b>d</b> is a diagonal matrix
	 */
	public double[][] v;
	/**
	 * The rate matrix must be diagonalized and represented in a
	 * <b>v d w</b> product, where <b>d</b> is a diagonal matrix
	 */
	public double[][] w;
	/**
	 * The rate matrix must be diagonalized and represented in a
	 * <b>v d w</b> product, where <b>d</b> is a diagonal matrix.
	 * Only the diagonal elements must be put into the one dimensional array <b>d</b>.
	 */
	public double[] d;
	
	/**
	 * This array tells the equilibrium distribution
	 */
	public double[] e;

	/**
	 * This array contains the potential characters of the alphabet.
	 */
	public char[] alphabet;

	/**
	 * Model parameters.
	 */
	public double params[];
	
	/**
	 * Tells whether your plugin is a nucleotide or protein model (or something else)
	 * when declared with the same name and type.
	 * Affects grouping of models in the menu.
	 */
	public static String type;
	
	/**
	 * Sets the name of the model plugin that will appear in the menu when declared
	 * with the same name and type.
	 */
	public static String menuName;
	
	/**
	 * This scoring scheme is used for constructing an
	 * initial tree.
	 */
	public SubstitutionScore attachedScoringScheme;
	
	
	/** Determines the span for proposing new parameters. Updated automatically during MCMC. */
	protected double span = 0.1;	
	public void updateSpan(double newSpan) { span = newSpan; }

	/**
	 * 
	 * If a model cannot accept the currently loaded sequences 
	 * (for example, because the sequences contains a character that is not part of the model)
	 * then it returns a negative number. Otherwise it returns a double value between 0 and 1,
	 * the better likes the model the sequences the greater the return number is.
	 * 
	 * @param r raw sequences
	 * @return the value telling how well the substitution model likes the sequences. 
	 */
	public abstract double acceptable(RawSequences r);
	
	/**
	 * Updates externally stored transition matrix using current values of v, w, d, e and current
	 *  edge length. Must always be called in the form m = updateTransitionMatrix(m, len) where m
	 *  can be null as well.
	 * @param transitionMatrix reference to matrix to store results, null to create new
	 * @param edgeLength evolutionary time
	 * @return updated (or newly created) transition matrix
	 */
	public double[][] updateTransitionMatrix(double[][] transitionMatrix, double edgeLength){
		int i, j, k;
		int n = d.length;
		
		if(transitionMatrix == null)
			transitionMatrix = new double[n][n];
		else {
			for(i = 0; i < transitionMatrix.length; i++)
				Arrays.fill(transitionMatrix[i], 0.0);
		}

		double exp_dtk;
		for(k = 0; k < n; k++){
			exp_dtk = Math.exp(d[k]*edgeLength);
			for(i = 0; i < n; i++){
				for(j = 0; j < n; j++){
					transitionMatrix[i][j] += v[i][k] * exp_dtk * w[k][j];
				}
			}
		}
		return transitionMatrix;
	}
	
	/**
	 * Abstract method to sample a new model parameter and update v,w,d,e accordingly.
	 * 
	 * @return log(fwd proposal probability/bwd proposal probability) [i.e. inverse Hastings ratio]
	 */
	public abstract double sampleParameter();
	
	/**
	 * Restores state just before the call of sampleParameter
	 */
	public abstract void restoreParameter();
	
	/**
	 * Represents model state as a String (usually parameter values concatenated)
	 * 
	 * @return the resulting String
	 */
	public abstract String print();
	
	/**
	 * Returns the colour associated with the character. This colour will be
	 * the background colour of the character in the printed alignments.
	 * 
	 * @param c character whose background colour is to be returned
	 */
	public abstract Color getColor(char c);
	
	public abstract char mostLikely(double[] seq);

	/**
	 * Retrieves the static field <code>type</code> of the class of the object the
	 * method is called on.
	 * 
	 * See <code>getType(Class cl)</code> for details.
	 */
	public String getType() {
		return getType(getClass());
	}
	
	/**
	 * Retrieves the static field <code>type</code> of a SubstitutionModel plug-in class.
	 * In the class it must be declared <code>public static String</code>.
	 * Searches the class hierarchy upwards until a class with the field declared is found.
	 * 
	 * @param cl class extending <code>SubstitutionModel</code>
	 * @return value of the static field <code>type</code> or "unknown" if the field 
	 * 		is null or empty
	 */
	@SuppressWarnings("unchecked")
	public static String getType(Class<? extends SubstitutionModel> cl) {
		String type = null;
		while(true) {
			try {
				type = (String)cl.getField("type").get(null);
				break;
			} catch (Exception e) {}
			cl = (Class<? extends SubstitutionModel>)cl.getSuperclass();
		}
		if(type == null || type.equals(""))
			type = "unknown";
		return type;
	}

	/**
	 * Retrieves the static field <code>menuName</code> of the class of the object the
	 * method is called on.
	 * 
	 * See <code>getMenuName(Class cl)</code> for details.
	 */
	public String getMenuName() {
		return getMenuName(getClass());
	}
	
	/**
	 * Retrieves <code>menuName</code> static field of a SubstitutionModel plug-in class.
	 * In the class it must be declared <code>public static String</code>.
	 * 
	 * @param cl class extending <code>SubstitutionModel</code>
	 * @return value of the static field <code>menuName</code> or the simple name of
	 *   <code>cl</code> if the field doesn't exist or is null
	 */
	public static String getMenuName(Class<? extends SubstitutionModel> cl) {
		String name = null;
		try {
			name = (String)cl.getField("menuName").get(null);
		} catch (Exception e) {}
		if(name == null)
			name = cl.getSimpleName();
		return name;
	}
	
	/**
	 * Returns the prior for the substitution parameters. 
	 * By default the prior is uniform and this function returns 1.
	 * @return prior for the parameters of the substitution model.
	 */
	public double getPrior(){
		return 1;
	}
	
	
}