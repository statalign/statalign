package statalign.io.input.plugins;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import statalign.exceptions.ExceptionNonFasta;
import statalign.io.RawSequences;
import statalign.io.input.FileFormatReader;
import statalign.postprocess.utils.RNAFoldingTools;

/**
 * 
 * Class to read files in Fasta format.
 * 
 * @author novak
 *
 */
public class FastaReader extends FileFormatReader {

	private int errors;

	/**
	 * Reads the contents (aligned/non-aligned sequences) of the given data source in
	 * Fasta format.
	 * 
	 * @param reader Data source
	 * @return RawSequences representation of the contents
	 * @throws IOException if an I/O error occurs
	 * @throws ExceptionNonFasta 
	 * @throws ExceptionNonRNA 
	 */
	@Override
	public RawSequences read(Reader reader) throws IOException, ExceptionNonFasta {
		RawSequences result = new RawSequences();
		BufferedReader br = new BufferedReader(reader);

		String line;
		boolean inSeq = false;
		StringBuilder actSeq = new StringBuilder();
		int errors = 0;
		while(true) {
			line = br.readLine();
			if(line != null && line.length() == 0)
				continue;
			if(line == null || line.charAt(0) == '>') {
				if(inSeq) {
					if(actSeq.length() == 0) {
						errors++;
						result.seqNames.remove(result.seqNames.size()-1);
					} else {
						result.sequences.add(actSeq.toString());
					}
				}
				if(line == null)
					break;
				actSeq.setLength(0);
				inSeq = true;
				int start = 1, index;
				if((index = line.indexOf(' ', 1)) == 1) {
					index = line.indexOf(' ', 2);
					start = 2;
				}
				if(index == -1)
					line = line.substring(start);
				else{
					//line = line.substring(start, index);
					line = line.substring(start);			
				}
				line = line.replaceAll("[ \t]+", "_");
				//line.replaceAll(" ", "_");
				line = line.replaceAll("\\(", "{");
				line = line.replaceAll("\\)", "}");
//				System.out.println("new name: "+line);
				result.seqNames.add(line);
			} else if(inSeq) {
				
				
				int len = line.length();
				char ch;
				for(int i = 0; i < len; i++) {
					if(!Character.isWhitespace(ch = line.charAt(i))) {
						if(Character.isLetter(ch) || ch == '-') {
							actSeq.append(ch);
							if(ch != '-') {
								result.alphabet += ch;
							}
						}
							
						else
							throw new ExceptionNonFasta("Input valid sequences!");
					}
				}
			} else {
				errors++;
			}
		}


//		for(int i = 0; i < result.sequences.size(); i++){
//			System.out.println(">"+result.seqNames.get(i)+"\n"+result.sequences.get(i));
//		}
		if(errors > 0)
			System.out.println("Errors: "+errors);
//		else
//			System.out.println("FastaReader: successfully read "+result.sequences.size()+" sequences.");
		
		return result;
	}
	
	/**
	 * Return the error count of the last read operation.
	 * 
	 * @return number of errors
	 */
	public int getErrors() {
		return errors;
	}

	/**
	 * Only for testing/debugging purposes: reads the given Fasta file
	 * 
	 * @param args First element must be the name of the Fasta file to read
	 * @throws IOException if an I/O error occurs
	 * @throws ExceptionNonFasta 
	 * @throws ExceptionNonRNA 
	 */
	public static void main(String[] args) throws IOException, ExceptionNonFasta {
		FastaReader f = new FastaReader();
		f.read(args[0]);
	}
}