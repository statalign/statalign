package statalign.io;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents simplified tertiary structures of a set of sequences (primarily proteins),
 * only one 3d-coordinate per character (amino acid) is stored. For amino acids, this is
 * normally the spatial location of the alpha-carbon atom.
 * 
 * @author novak
 */
public class ProteinSkeletons {
	
	/** names of sequences */
	public List<String> names = new ArrayList<String>();
	
	/** coordinates */
	public List<List<double[]>> coords = new ArrayList<List<double[]>>();
}
