package statalign.io.input;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import statalign.io.RawSequences;

/**
 * 
 * Ancestor for sequence/alignment file reader classes that handle one specific format.
 * 
 * @author novak
 *
 */
public abstract class FileFormatReader {

	/**
	 * Reads a sequence/alignment file and constructs a RawSequences representation.
	 * 
	 * @param file File to read
	 * @return RawSequences object containing the data read
	 * @throws IOException when an I/O error occurs
	 */
	public RawSequences read(File file) throws IOException {
		return read(new FileReader(file));
	}
	
	/**
	 * Equivalent to read(new File(fileName)).
	 * 
	 * @param fileName Name of file to read
	 * @return RawSequences object containing the data read
	 * @throws IOException when an I/O error occurs
	 */
	public RawSequences read(String fileName) throws IOException {
		return read(new FileReader(fileName));
	}

	/**
	 * Reads sequence/alignment from a Reader and constructs a RawSequences representation.
	 * 
	 * @param reader Data source
	 * @return RawSequences object containing the data read
	 * @throws IOException when an I/O error occurs
	 */
	public abstract RawSequences read(Reader reader) throws IOException;
	
}
