package statalign.model.ext.plugins.structalign;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.model.ext.plugins.structalign.Transformation;

public class Funcs{

	/** Calculates the cross product of 2 vectors
	 * 
	 * @param a first vector
	 * @param b second vector
	 * @return cross product
	 */
	/*static RealVector crossProduct(RealVector a, RealVector b){
		RealMatrix skew = new Array2DRowRealMatrix(
				new double[][] {{0, -a.getEntry(2), a.getEntry(1)}, {a.getEntry(2), 0, -a.getEntry(0)},
						{-a.getEntry(1), a.getEntry(0), 0}});
		return skew.operate(b);
	}		*/

	/** 
	 * Calculates the rotation matrix from an axis and angle
	 * with axis x and angle r, the rotation matrix is
	 * cos(r)I + sin(r)cross(x) + (1-cos(r))outer(x)
	 * where cross(x) is the crossproduct matrix of x and outer(x) is the outer (tensor) product
	 * 
	 * @return Transformation object with rotation matrix calculated from axis and angle
	 */
	/*static RealMatrix calcRotationMatrix(RealVector axis, double rot){
		RealMatrix outer = axis.outerProduct(axis);
		// use transpose of cross product matrix (as given on wikipedia) because
		// we post-multiply by rotation matrix (other components of final step are symmetric)
		RealMatrix crossTranspose = new Array2DRowRealMatrix(new double[][] { 
				{0, axis.getEntry(2), -axis.getEntry(1)},
				{-axis.getEntry(2), 0, axis.getEntry(0)},
				{axis.getEntry(1), -axis.getEntry(0), 0} 
		});
		RealMatrix Icos = new Array2DRowRealMatrix(new double[3][3]);
		for(int i = 0; i < 3; i++)
			Icos.setEntry(i, i, Math.cos(rot));
		crossTranspose = crossTranspose.scalarMultiply(Math.sin(rot));
		outer = outer.scalarMultiply(1 - Math.cos(rot));
		return Icos.add(crossTranspose).add(outer);
	} */

	/**
	 * For an n X 3 coordinate matrix, calculate the 1 X 3 mean vector
	 * @param A - coordinate matrix
	 * @return mean vector
	 */
	
	static RealVector meanVector(RealMatrix A){
		RealVector mean = new ArrayRealVector(new double[3]);
		for(int i = 0; i < 3; i ++){
			for(int j = 0; j < A.getColumn(0).length; j++)
				mean.addToEntry(i, A.getEntry(j, i));
			mean.setEntry(i, mean.getEntry(i) / A.getColumn(0).length);
		}
		return mean;
	}


	public Vertex sampleVertex(Tree tree){
		/* Vertex v;
		int n = tree.vertex.length;		// number of vertices
		int l = coords.length;			// number of leaves
		
		
		List<Integer> inds = findRefSubtrees(tree, 0);	// returns indices of all ancestor vertices of reference protein
		if(inds.size() + l < n){
			int prop = Utils.generator.nextInt(n-l) + l;		// don't choose a leaf vertex			
			while(inds.contains(new Integer(prop)))
				prop = Utils.generator.nextInt(n-l) + l;		// don't choose a subtree containing reference protein
			v = tree.vertex[prop];
			return v;
		} else
			return tree.root; */
		return tree.vertex[Utils.generator.nextInt(tree.vertex.length - 1)]; // -1 excludes root
	}
	
	public ArrayList<Integer> findRefSubtrees(Tree tree, int refInd){
		ArrayList<Integer> inds = new ArrayList<Integer>(0);
		moveUp(tree.vertex[refInd].parent, inds);
		return inds;
	}
	
	public void moveUp(Vertex v, List<Integer> inds){
		inds.add(v.index);
		if(v.parent != null)
			moveUp(v.parent, inds);
	}
	
	public ArrayList<Integer> collectLeaves(Vertex v){
		ArrayList<Integer> inds = new ArrayList<Integer>(0);
		moveDown(v, inds);
		return inds;
	}
	
	public void moveDown(Vertex v, List<Integer> inds){
		if(v.left != null){
			moveDown(v.left, inds);
			moveDown(v.right, inds);
		}
		else
			inds.add(v.index);
	}
	
	public static void writeRotationFiles(Tree tree, double[][][] coords, 
			double[][] xlats, double[][] axes, double[] angles){
		try{
			FileWriter fstream = new FileWriter("names.txt");
			BufferedWriter out = new BufferedWriter(fstream);
			for(int i = 0; i < coords.length; i++)
				out.write(tree.vertex[i].name + "\n");
			out.close();
	
			FileWriter fstream2 = new FileWriter("trans.txt");
			BufferedWriter out2 = new BufferedWriter(fstream2);
			for(int i = 0; i < coords.length; i++){
				for(int j = 0; j < 3; j++)
					out2.write(xlats[i][j] + "\t");
				out2.write("\n");
			}
			out2.close();
			
			FileWriter fstream3 = new FileWriter("rots.txt");
			BufferedWriter out3 = new BufferedWriter(fstream3);
			for(int i = 0; i < coords.length; i++){
				Transformation trans = new Transformation(axes[i], angles[i], xlats[i]);
				RealMatrix temp = trans.getRealRotation();
				for(int j = 0; j < 3; j++){
					for(int k = 0; k < 3; k++)
						out3.write(temp.getEntry(j, k) + "\t");
					out3.write("\n");
				}
			}
			out3.close();
			
		} catch (Exception e){
			System.err.println("File writing error: " + e.getMessage());
		}
		
	}
	
	public void initLSRotations(Tree tree, double[][][] coords, 
			double[][] xlats, double[][] axes, double[] angles){
		String[] align = tree.getState().getLeafAlign();
		String ref = align[0];
		for(int i = 1; i < align.length; i++){
			ArrayList<Integer> refInds = new ArrayList<Integer>(0);
			ArrayList<Integer> otherInds = new ArrayList<Integer>(0);
			String other = align[i];
			
			int r = 0, o = 0;
			for(int j = 0; j < align[0].length(); j++){
				if(ref.charAt(j) != '-' & other.charAt(j) != '-'){
					refInds.add(r);
					otherInds.add(o);
				}
				r += ref.charAt(j) != '-' ? 1 : 0;
				o += other.charAt(j) != '-' ? 1 : 0;
			}
			
			double[][] refSub = getRowSub(coords[0], refInds);
			double[][] otherSub = getRowSub(coords[i], otherInds);
			
			Transformation trans = new Transformation();
			
			trans.lsTrans(new Array2DRowRealMatrix(refSub), new Array2DRowRealMatrix(otherSub));
			axes[i] = trans.rotMatrix.getAxis().toArray();
			xlats[i] = trans.xlat.toArray();
			angles[i]  = trans.rotMatrix.getAngle();
		}
	}
	
	public double[][] getRowSub(double[][] coord, ArrayList<Integer> rows){
		double[][] sub = new double[rows.size()][coord[0].length];
		
		for(int i = 0; i < sub.length; i++)
			for(int j = 0; j < sub[0].length; j++)
				sub[i][j] = coord[rows.get(i)][j];
		return sub;
	}

}


