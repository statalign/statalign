package statalign.io;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents simplified tertiary structures of a set of sequences (primarily proteins),
 * only one 3d-coordinate per character (amino acid) is stored. For amino acids, this is
 * normally the spatial location of the alpha-carbon atom.
 * 
 * @author novak, herman
 */
public class ProteinSkeletons implements DataType {
	
	@Override
	public boolean perSequenceData() {
		return true;
	}
	
	public RawSequences seqs = null;
	
	public RawSequences getSeqs() { return seqs; }
	public void setSeqs(RawSequences rs) { seqs = rs; }
	
	@Override
	public void removeDataAssociatedWith(String sequenceName) {
		for (int i=0; i<names.size(); i++) {
			if (names.get(i).toUpperCase().equals(sequenceName.toUpperCase())) {
				names.remove(i);
				coords.remove(i);
				bFactors.remove(i);
			}
		}
	}
	@Override
	public String getSummaryAssociatedWith(String sequenceName) {
		for (int i=0; i<names.size(); i++) {
			if (names.get(i).toUpperCase().equals(sequenceName.toUpperCase())) {
				return (names.get(i)+".pdb ("+coords.get(i).size()+" residues)");
			}
		}
		return "";
	}
	/** names of sequences */
	public List<String> names = new ArrayList<String>();		
	
	/** coordinates */
	public List<List<double[]>> coords = new ArrayList<List<double[]>>();

	public List<List<Double> > bFactors = new ArrayList<List<Double> >();
}
