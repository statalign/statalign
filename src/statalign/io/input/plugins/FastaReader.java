package statalign.io.input.plugins;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Arrays;
import java.util.List;

import statalign.io.RawSequences;
import statalign.io.input.DataReader;
import statalign.io.input.IllegalFormatException;

/**
 * 
 * Class to read files in FASTA format.
 * 
 * @author novak
 *
 */
public class FastaReader extends DataReader {
	
	private static boolean[] allowedChars;
	
	static {
		allowedChars = new boolean[255];
		for(int ch = 'a'; ch <= 'z'; ch++)
			allowedChars[ch] = true;
		for(int ch = 'A'; ch <= 'Z'; ch++)
			allowedChars[ch] = true;
		allowedChars['-'] = true;
	}

	@Override
	public List<String> supportedExtensions() {
		return Arrays.asList(new String[] { "fsa","fas","fasta" });
	}
	
	/**
	 * Reads the contents (aligned/non-aligned sequences) of the given data source in
	 * Fasta format.
	 * 
	 * @param reader Data source
	 * @return RawSequences representation of the contents
	 * @throws IOException if an I/O error occurs
	 */
	@Override
	public RawSequences read(Reader reader) throws IOException {
		RawSequences result = new RawSequences();
		BufferedReader br = new BufferedReader(reader);

		String line;
		boolean inSeq = false;
		StringBuilder actSeq = new StringBuilder();
		boolean[] seen = new boolean['Z'-'A'+1];
		while(true) {
			line = br.readLine();
			if(line != null && line.length() == 0)
				continue;
			if(line == null || line.charAt(0) == '>') {
				if(inSeq) {
					if(actSeq.length() == 0) {
						throw new IllegalFormatException("FastaReader: empty sequence "+result.seqNames.get(result.seqNames.size()-1));
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
						if(allowedChars[ch]) {
							actSeq.append(ch);
							if(ch != '-')
								seen[Character.toUpperCase(ch)-'A'] = true;
						} else {
							throw new IllegalFormatException("FastaReader: illegal character "+ch);
						}
					}
				}
			} else {
				throw new IllegalFormatException("FastaReader: data without sequence name");
			}
		}
		if(!inSeq)
			throw new IllegalFormatException("FastaReader: empty file");
		
		StringBuilder alpha = new StringBuilder();
		for(int ch = 0; ch <= 'Z'-'A'; ch++)
			if(seen[ch])
				alpha.append((char)(ch+'A'));
		result.alphabet = alpha.toString();
		
		return result;
	}
	
}