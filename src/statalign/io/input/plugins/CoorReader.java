package statalign.io.input.plugins;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import statalign.io.ProteinSkeletons;
import statalign.io.input.DataReader;

public class CoorReader extends DataReader {

	/** minimum dimension allowed */
	private int allowDimMin = 3;
	/** maximum dimension allowed */
	private int allowDimMax = 3;
	
	public CoorReader() {
	}
	
	public void setAllowedDim(int dim) {
		allowDimMin = allowDimMax = dim;
	}

	@Override
	public List<String> supportedExtensions() {
		return Arrays.asList(new String[] { "coor" });
	}
	
	@Override
	public ProteinSkeletons read(Reader reader) throws IOException {
		ProteinSkeletons data = new ProteinSkeletons();
		BufferedReader br = new BufferedReader(reader);
		String line;
		int cur = -1;
		try {
			while((line = br.readLine()) != null) {
				if(line.charAt(0) == '>') {		// new sequence
					data.names.add(line.substring(1));
					data.coords.add(new ArrayList<double[]>());
					cur++;
				} else {
					if(cur < 0)
						throw new WrongFormatException();
					String[] tokens = line.trim().split("[:,; \t]+");
					if(tokens.length < allowDimMin || tokens.length > allowDimMax)
						throw new WrongFormatException();
					double[] vals = new double[tokens.length];
					for(int i = 0; i < tokens.length; i++)
						vals[i] = Double.parseDouble(tokens[i]);
					data.coords.get(cur).add(vals);
				}
			}
		} catch (NumberFormatException e) {
			throw new WrongFormatException();
		}
		return null;
	}
	
	
}
