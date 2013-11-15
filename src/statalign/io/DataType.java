package statalign.io;

/**
 * Ancestor for all input/output data types in StatAlign.
 * 
 * @author novak
 */
public interface DataType {
	
	/**
	 * @return <tt>true</tt> if the DataType contains a List of
	 * items, each of which is associated with one of the input sequences, 
	 * for example protein structure coordinates. 
	 */
	public abstract boolean perSequenceData();
	public RawSequences getSeqs();
	public void setSeqs(RawSequences rs);
	public abstract String getSummaryAssociatedWith(String sequenceName);
	public abstract void removeDataAssociatedWith(String sequenceName);
	
}
