package statalign.exceptions;

public class ExceptionNonFasta extends Exception {

	/**
	 * 
	 * Exception is called if input files/sequences are not in Fasta format.
	 * 
	 * @author Preeti Arunapuram
	 */
	private static final long serialVersionUID = 1L;

	private String errorMessage;
	
	public ExceptionNonFasta(String errorMessage) {
		this.errorMessage = errorMessage;
	}
	
	public String getError() {
		return errorMessage;
	}
	
	public String getMessage() {
		return "Can only input .fasta files!";
	}
	
}
