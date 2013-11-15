package statalign.io.input.plugins;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import statalign.io.ProteinSkeletons;
import statalign.io.RawSequences;
import statalign.io.input.DataReader;
import statalign.io.input.IllegalFormatException;

public class PDBReader extends DataReader {
		
	public PDBReader() {
	}		

	@Override
	public List<String> supportedExtensions() {
		return Arrays.asList(new String[] { "pdb", "PDB" });
	}
	
	@Override
	public ProteinSkeletons read(Reader reader) throws IOException {
		ProteinSkeletons data = new ProteinSkeletons();
		String sequence = "";
		BufferedReader br = new BufferedReader(reader);
		String line = null;
		int lineNumber = 0;
		try {
			data.names.add(new String());
			data.coords.add(new ArrayList<double[]>());
			data.bFactors.add(new ArrayList<Double>());
			HashMap<Integer,Boolean> seen = new HashMap<Integer,Boolean>();
			while((line = br.readLine()) != null) {				
				if(++lineNumber == 1) {					
					data.names.set(0,line.substring(62,66).toLowerCase());
					continue;
				}
				if(line.startsWith("ATOM")) {
					int resn = Integer.parseInt(line.substring(22, 26).replaceAll("\\s+",""));
					if (seen.containsKey(resn)) continue; // so we don't include alt. conformers					
					if (line.substring(12, 16).equals(" CA ")) {
						sequence += oneLetter(line.substring(17, 20));		
						double x = Double.parseDouble(line.substring(30, 38));
						double y = Double.parseDouble(line.substring(38, 46));
						double z = Double.parseDouble(line.substring(46, 54));
				        
						data.coords.get(0).add(new double[]{x,y,z});
						if (line.length() > 61) {							
							data.bFactors.get(0).add(Double.parseDouble(line.substring(60,66)));
						}
						seen.put(resn,true);
					}
				}
			}						
		} catch (NumberFormatException e) {
			throw new IllegalFormatException("PDBReader: number formatting error: "+e.getMessage());
		}
		if(lineNumber < 1)
			throw new IllegalFormatException("PDBReader: empty file");
		if(data.coords.get(0).size() == 0)
			throw new IllegalFormatException("PDBReader: sequence "+data.names.get(0)+" is without a structure");
		
		data.seqs = new RawSequences(sequence,data.names.get(0));
		for (int i=0; i<data.bFactors.get(0).size(); i++) {				
			System.out.println(i+"\t"+data.seqs.getSequence(0).charAt(i)+"\t"
									+data.coords.get(0).get(i)[0]+"\t"
									+data.coords.get(0).get(i)[1]+"\t"
									+data.coords.get(0).get(i)[2]+"\t"
									+data.bFactors.get(0).get(i));
		}
		return data;
	}
	
	public HashMap<String,String> mp31 = null;
	
	public String oneLetter(String threeLetter) {
        if (mp31 == null) initialiseMp31();
        return mp31.get(threeLetter);
	}
	
	private void initialiseMp31() {
		mp31 = new HashMap<String,String>();
		mp31.put("ALA","A");
        mp31.put("ARG","R");
        mp31.put("ASN","N");
        mp31.put("ASP","D");
        mp31.put("ASX","B");
        mp31.put("CYS","C");
        mp31.put("GLU","E");
        mp31.put("GLN","Q");
        mp31.put("GLX","Z");
        mp31.put("GLY","G");
        mp31.put("HIS","H");
        mp31.put("HSD","H");
        mp31.put("ILE","I");
        mp31.put("LEU","L");
        mp31.put("LYS","K");
        mp31.put("MET","M");
        mp31.put("PHE","F");
        mp31.put("PRO","P");
        mp31.put("SER","S");
        mp31.put("THR","T");
        mp31.put("TRP","W");
        mp31.put("TYR","Y");
        mp31.put("VAL","V");     
	}
}
