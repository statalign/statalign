package statalign.model.ext.plugins;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JToggleButton;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

import statalign.base.InputData;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.io.DataType;
import statalign.io.ProteinSkeletons;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;
import statalign.postprocess.PluginParameters;

public class StructAlign extends ModelExtension implements ActionListener {
	
	/** the command line identifier of this plugin */
	private static final String CMD_LINE_PLUGIN_ID = "structal";
	
	JToggleButton myButton;

	/** alpha-C atomic coordinate for each sequence and each residue */
	double[][][] coords;
	
	/** alpha-C atomic coordinates under the current set of rotations/translations */
	double[][][] rotCoords;
	
	/** axis of rotation for each sequence */
	double[][] axes;
	/** rotation angle for each protein along the rotation axis */
	double[] rots;
	/** translation vector for each protein */
	double[][] xlats;
	
	/** covariance matrix implied by current tree topology */
	double[][] fullCovar;
	
	/** parameters of structural drift */
	double theta = .1; // CONSTANT VALUE FOR NOW
	double sigma2 = 1; // SHOULD BE UPDATED WITH MCMC
	
	@Override
	public List<JComponent> getToolBarItems() {
		myButton = new JToggleButton(new ImageIcon(ClassLoader.getSystemResource("icons/protein.png")));
    	myButton.setToolTipText("Structural alignment mode (for proteins only)");
    	myButton.addActionListener(this);
    	myButton.setEnabled(true);
    	myButton.setSelected(false);
    	return Arrays.asList((JComponent)myButton);

	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		setActive(myButton.isSelected());
	}

	@Override
	public void setActive(boolean active) {
		super.setActive(active);
		System.out.println("StructAlign plugin is now "+(active?"enabled":"disabled"));
	}
	
	@Override
	public void init(ModelExtManager manager, PluginParameters params) {
		if(params != null && params.getParameter(CMD_LINE_PLUGIN_ID) != null) {
			// TODO parse plugin parameters
			setActive(true);
		}
	}
	
	@Override
	public void initRun(InputData inputData) throws IllegalArgumentException {
		HashMap<String, Integer> seqMap = new HashMap<String, Integer>();
		int i = 0;
		for(String name : inputData.seqs.seqNames)
			seqMap.put(name, i++);
		coords = new double[inputData.seqs.seqNames.size()][][];
		for(DataType data : inputData.auxData) {
			if(!(data instanceof ProteinSkeletons))
				continue;
			ProteinSkeletons ps = (ProteinSkeletons) data;
			for(i = 0; i < ps.names.size(); i++) {
				String name = ps.names.get(i);
				if(!seqMap.containsKey(name))
					throw new IllegalArgumentException("structalign: missing sequence or duplicate structure for "+name);
				int ind = seqMap.get(name);
				int len = inputData.seqs.sequences.get(ind).replaceAll("-", "").length();
				List<double[]> cl = ps.coords.get(i);
				if(len != cl.size())
					throw new IllegalArgumentException("structalign: sequence length mismatch with structure file for seq "+name);
				coords[ind] = new double[len][];
				for(int j = 0; j < len; j++)
					coords[ind][j] = Utils.copyOf(cl.get(j));
				seqMap.remove(name);
			}
		}
		if(seqMap.size() > 0)
			throw new IllegalArgumentException("structalign: missing structure for sequence "+seqMap.keySet().iterator().next());
	}
	
	@Override
	public double logLikeFactor(Tree tree) {
		// simplest (and slowest) approach: get alignment of all leaves and compute likelihood from there
		String[] align = tree.getState().getLeafAlign();
		fullCovar = calcFullCovar(tree);
		double logli = 0;
		int[] inds = new int[align.length];		// current char indices
		int[] col = new int[align.length];  // SOMETHING IS WRONG WITH THE WAY INDICES ARE HANDLED HERE
		for(int i = 0; i < align[0].length(); i++) {
			for(int j = 0; j < align.length; j++)
				col[j] = align[j].charAt(i) == '-' ? -1 : inds[j]++;
			logli += Math.log(columnContrib(col));
		}
		return logli;
	}
	
	/**
	 * Calculates the structural likelihood contribution of a single alignment column
	 * @param col the column, id of the residue for each sequence (or -1 if gapped in column)
	 * @return the likelihood contribution
	 */
	public double columnContrib(int[] col) {
		// count the number of ungapped positions in the column
		int numMatch = 0;
		for(int i = 0; i < col.length; i++){
			if(col[i] != -1)
				numMatch++;
		}
		if(numMatch == 0)  // CHRIS: this shouldn't happen, but some columns are all gaps
			return 0;
		// collect indices of ungapped positions
		int[] notgap = new int[numMatch];
		int j = 0;
		for(int i = 0; i < col.length; i++)
			if(col[i] != -1)
				notgap[j++] = i;
		
		// extract covariance corresponding to ungapped positions
		double[][] subCovar = getSubMatrix(fullCovar, notgap, notgap);
		// create normal distribution with mean 0 and covariance subCovar
		MultivariateNormalDistribution multiNorm = new MultivariateNormalDistribution(new double[numMatch], subCovar);
		
		double li = 1;
		double[] vals = new double[numMatch];
		// CHRIS: FOR TESTING ONLY, REMOVE WHEN rotCoords PROPERLY FILLED
		rotCoords = coords;
		// loop over all 3 coordinates
		for(j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++)
				vals[i] = rotCoords[notgap[i]][col[notgap[i]]][j];
			li *= multiNorm.density(vals);
		}
		return li;
	}

	/**
	 * extracts the specified rows and columns of a 2d array
	 * @param matrix, 2d array from which to extract; rows, rows to extract; cols, columns to extract
	 * @return submatrix
	 */
	public double[][] getSubMatrix(double[][] matrix, int[] rows, int[] cols) {
		double[][] submat = new double[rows.length][cols.length];
		for(int i = 0; i < rows.length; i++)
			for(int j = 0; j < cols.length; j++)
				submat[i][j] = matrix[rows[i]][cols[j]];
		return submat;
	}
	
	/**
	 * return the full covariance matrix for the tree topology and branch lengths
	 */	
	public double[][] calcFullCovar(Tree tree) {
		// I'm assuming that tree.names.length is equal to the number of vertices here
		double[][] distMat = new double[tree.names.length][tree.names.length];
		calcDistanceMatrix(tree.root, distMat);
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i; j < tree.names.length; j++)
				distMat[i][j] = sigma2 / (2 * theta ) * Math.exp(-theta * distMat[i][j]);
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i + 1; j < tree.names.length; j++)
				distMat[j][i] = distMat[i][j];
		return distMat;
	}
	
	/**
	 * recursive algorithm to traverse tree and calculate distance matrix between leaves 
	 */		
	public int[] calcDistanceMatrix(Vertex vertex, double[][] distMat){
		int[] subTree = new int[distMat.length + 1];
		// either both left and right are null or neither is
		if(vertex.left != null){
			int[] subLeft  = calcDistanceMatrix(vertex.left, distMat);
			int[] subRight = calcDistanceMatrix(vertex.right, distMat);
			int i = 0;
			while(subLeft[i] > -1){
				subTree[i] = subLeft[i];
				i++;
			}
			for(int j = 0; i+j < subTree.length; j++)
				subTree[i+j] = subRight[j];
		}
		else{
			subTree[0] = vertex.index;
			for(int j = 1; j < subTree.length; j++)
				subTree[j] = -1;
		}
		addEdgeLength(distMat, subTree, vertex.edgeLength);
		return subTree;
	}
		
	// adds the length of the current edge to the distance between all leaves
	// of a subtree to all other leaves
	// 'rows' contains the indices of vertices in the subtree
	public void addEdgeLength(double[][] distMat, int[] subTree, double edgeLength){
		
		int i = 0;
		while(subTree[i] > -1){
			for(int j = 0; j < distMat.length; j++)  
				distMat[subTree[i]][j] = distMat[subTree[i]][j] + edgeLength;
			i++;		
		}
			
		// edge length should not be added to distance between vertices in the subtree
		// subtract the value from these entries of the distance matrix
		i = 0;
		while(subTree[i] > -1){
			int j = 0;
			while(subTree[j] > -1){
				distMat[subTree[i]][subTree[j]] = 
						distMat[subTree[i]][subTree[j]] - edgeLength;
				j++;
			}
			i++;
		}
	}


	@Override
	public int getParamChangeWeight() {
		// TODO test converge and tune value
		return 25;
	}
	
	@Override
	public double proposeParamChange(Tree tree) {
		return 0;
	}

	@Override
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {
		if(!accepted)
			return;
		// TODO recalculate things
	}
	
	@Override
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {
		if(!accepted)
			return;
		// TODO recalculate
	}

}
