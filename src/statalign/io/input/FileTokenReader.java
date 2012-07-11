package statalign.io.input;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * 
 * Class that allows reading the contents of a file token by token, without worrying about
 * line ends and whitespace. E.g. files containing ten integers in a single line or ten
 * lines can be handled in a uniform way. Inspired by the C call fscanf(file, "%d", &value).
 * 
 * Tokens can be integers, doubles or strings. Separators are any combination of whitespace
 * characters as recognized by the regular expression \s+ .
 * 
 * @author novak
 *
 */
public class FileTokenReader {
	private BufferedReader file;
	private String tokens[];
	private int ind;
	boolean EOF;

	/**
	 * Constructs a FileTokenReader reading a given file.
	 * 
	 * @param fname Name of the file to read
	 * @throws IOException if an I/O error occurs
	 */
	 public FileTokenReader(String fname) throws IOException {
		ClassLoader cl = getClass().getClassLoader();
		file = new BufferedReader(new InputStreamReader(cl.getResource(fname).openStream()));
//		file = new BufferedReader(new FileReader(fname));
		tokens = null;
		ind = 0;
		EOF = false;
	}

	protected boolean readTokens() throws IOException {
		String line;
		while(tokens == null || ind == tokens.length || tokens[ind].length() == 0) {
			if((line = file.readLine()) == null) {
				EOF = true;
				return false;
			}
			tokens = line.split("\\s+");
			ind = 0;
		}
		return true;
	}

	/**
	 * Reads the next token and tries to interpret it as a double value.
	 * 
	 * @return Double value of the token or null if end of file has been reached
	 * @throws IOException if an I/O error occurs
	 * @throws NumberFormatException if next token cannot be converted to Double
	 */
	public Double readDbl() throws IOException {
		return readTokens()?Double.valueOf(tokens[ind++]):null;
	}
	
	/**
	 * Reads the next token and tries to interpret it as an integer value.
	 * 
	 * @return Integer value of the token or null if end of file has been reached
	 * @throws IOException if an I/O error occurs
	 * @throws NumberFormatException if next token cannot be converted to Integer
	 */
	public Integer readInt() throws IOException {
		return readTokens()?Integer.valueOf(tokens[ind++]):null;
	}
	
	/**
	 * Reads the next token and returns it as a String.
	 * 
	 * @return Next token as String or null if end of file has been reached
	 * @throws IOException if an I/O error occurs
	 */
	public String readStr() throws IOException {
		return readTokens()?tokens[ind++]:null;
	}
	
	/**
	 * Closes the underlying file.
	 * 
	 * @throws IOException if an I/O error occurs
	 */
	public void close() throws IOException {
		file.close();
	}
	
}
