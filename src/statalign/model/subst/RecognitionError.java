/**
 * 
 */
package statalign.model.subst;

/**
 * An error message that is thrown if a model cannot accept a set of sequences.
 * The error is handled, and it is put into an ErrorMessage.
 * 
 * @author miklos
 *
 */
public class RecognitionError extends Error{

	private static final long serialVersionUID = 1L;
	/**
	 * The message of the error
	 */
	public String message;
	
	/**
	 * Trivial constructor, it sets the message 
	 * @param message The message to be set.
	 */
	public RecognitionError(String message){
		this.message = message;
	}

}
