package statalign.io.input.plugins;

import java.io.BufferedReader;
import java.io.FileWriter;
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
					if (line.startsWith("HEADER") && line.length() > 62) {
						if (line.length() > 66)
							data.names.set(0,line.substring(62,67).toLowerCase().replaceAll("\\s+",""));
						else 
							data.names.set(0,line.substring(62,66).toLowerCase());
					}
					else { // If HEADER lines have been removed from PDB						
						if (filename.length() > 4) 
							data.names.set(0,filename.toLowerCase().substring(0,5));
						else 
							data.names.set(0,filename.toLowerCase().substring(0,4));
						
						// then use the filename to define the name for this structure
					}
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
//		for (int i=0; i<data.bFactors.get(0).size(); i++) {				
//			System.out.println(i+"\t"+data.seqs.getSequence(0).charAt(i)+"\t"
//									+data.coords.get(0).get(i)[0]+"\t"
//									+data.coords.get(0).get(i)[1]+"\t"
//									+data.coords.get(0).get(i)[2]+"\t"
//									+data.bFactors.get(0).get(i));
//		}
		return data;
	}
	
	public static void writePDB(double[][][] coors, String[] seqs, String[] names, FileWriter fw) throws IOException {
		initialiseMp13();
		String format = "ATOM  %5d  CA  %3s %1c%4d    %8.3f%8.3f%8.3f\n";
		String format2 = "TER   %5d      %3s %1c%4d\n";
		try {
			for (int i=0; i<coors.length; i++) {
				if (coors[i]==null) continue; 
				char chain = (char) ('A' + i);
				fw.write("HEADER "+names[i]+"\n");
				int j=0;
				for (; j<coors[i].length; j++) {						
					fw.write(String.format(format, j+1,mp13.get(seqs[i].substring(j,j+1)),chain,j+1,coors[i][j][0],coors[i][j][1],coors[i][j][2]));
				}			
				fw.write(String.format(format2,j+1,mp13.get(seqs[i].substring(j-1,j)),chain,j));
			}
			fw.write("END\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public HashMap<String,String> mp31 = null;
	
	public String oneLetter(String threeLetter) {
        if (mp31 == null) initialiseMp31();
        String one = mp31.get(threeLetter);
        if (one == null) throw new RuntimeException("Unrecognised amino acid: "+threeLetter);
        return one;
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
	
	public static HashMap<String,String> mp13 = null; 
	private static void initialiseMp13() {
		if (mp13 != null) return;
		mp13 = new HashMap<String,String>();
		mp13.put("A","ALA");
        mp13.put("R","ARG");
        mp13.put("N","ASN");
        mp13.put("D","ASP");
        mp13.put("B","ASX");
        mp13.put("C","CYS");
        mp13.put("E","GLU");
        mp13.put("Q","GLN");
        mp13.put("Z","GLX");
        mp13.put("G","GLY");
        mp13.put("H","HIS");
        mp13.put("H","HSD");
        mp13.put("I","ILE");
        mp13.put("L","LEU");
        mp13.put("K","LYS");
        mp13.put("M","MET");
        mp13.put("F","PHE");
        mp13.put("P","PRO");
        mp13.put("S","SER");
        mp13.put("T","THR");
        mp13.put("W","TRP");
        mp13.put("Y","TYR");
        mp13.put("V","VAL");
        mp13.put("X","XXX");
	}
}
