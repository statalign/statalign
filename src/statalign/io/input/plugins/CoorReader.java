package statalign.io.input.plugins;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import statalign.io.ProteinSkeletons;
import statalign.io.input.DataReader;
import statalign.io.input.IllegalFormatException;
import statalign.model.ext.plugins.StructAlign;

public class CoorReader extends DataReader {

	/** minimum dimension allowed */
	private int allowDimMin = 3;
	/** maximum dimension allowed */
	private int allowDimMax = 4;
	
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
		String line = null;
		int cur = -1;
		try {
			while((line = br.readLine()) != null) {
				if(line.charAt(0) == '>') {		// new sequence
					if(cur >= 0 && data.coords.get(cur).size() == 0)
						throw new IllegalFormatException("CoorReader: structure "+data.names.get(cur)+" contains no atoms.");
					data.names.add(line.substring(1).trim());
					data.coords.add(new ArrayList<double[]>());
					data.bFactors.add(new ArrayList<Double>());
					cur++;
				} else {
					line = line.trim();
					if(line.isEmpty())
						continue;
					if(cur < 0)
						throw new IllegalFormatException("CoorReader: unnamed structure detected");
					String[] tokens = line.split("[:,; \t]+");
					if(tokens.length < allowDimMin || tokens.length > allowDimMax)
						throw new IllegalFormatException("CoorReader: number of dimensions incorrect");
										
					if (tokens.length >= 3) {
						double[] vals = new double[3];
						for(int i = 0; i < 3; i++) {
							vals[i] = Double.parseDouble(tokens[i]);
						}
						data.coords.get(cur).add(vals);
					}
					if (tokens.length == 4) {
						double b = Double.parseDouble(tokens[3]);
						data.bFactors.get(cur).add(b);						
					}									
				}
			}
		} catch (NumberFormatException e) {
			throw new IllegalFormatException("CoorReader: number formatting error: "+e.getMessage());
		}
		if(cur < 0)
			throw new IllegalFormatException("CoorReader: empty file");
		if(data.coords.get(cur).size() == 0)
			throw new IllegalFormatException("CoorReader: structure "+data.names.get(cur)+" contains no atoms.");
		
		return data;
	}
	
	
}
