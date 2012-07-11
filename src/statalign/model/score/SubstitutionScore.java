package statalign.model.score;

/**
 * The ancestral class of substitution scores
 * 
 * @author novak
 *
 */
public abstract class SubstitutionScore{

	/**
	 * The distance matrix for subtitutions.
	 * It is a <b>distance</b> and not a similarity matrix!
	 * Distances are used to construct a NJ tree as initial tree.
	 */
  public int dist[][];
  
  /**
   * This array tells the mapping between characters and indices in the distance array.
   * For details, see the Blosum62.java and DNAscore.java, how to fill in this array in
   * implemented subclasses.
   * 
   * It is a two dimensional array to allow ambiguous characters in DNA and RNA sequences.
   */
  public int which[][];
 
  
}