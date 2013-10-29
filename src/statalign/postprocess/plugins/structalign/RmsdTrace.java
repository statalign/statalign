package statalign.postprocess.plugins.structalign;

import java.awt.BorderLayout;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import statalign.base.InputData;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.mcmc.McmcMove;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.plugins.StructAlign;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.StructAlignTraceGUI;


public class RmsdTrace extends Postprocess {

	
	public StructAlign structAlign;
	
	public double[][] distanceMatrix;

	public RmsdTrace() {
		screenable = false;
		outputable = true;
		postprocessable = true;
		postprocessWrite = false;
		selected = false;
		active = false;
	}

	
	@Override
	public String getFileExtension() {
		return "rmsd";
	}

	@Override
	public void setSampling(boolean enabled) {
	}
	
	@Override
	public void init(ModelExtManager modelExtMan) {
		for(ModelExtension modExt : modelExtMan.getPluginList()) {
			if(modExt instanceof StructAlign) {
				structAlign = (StructAlign) modExt;
				structAlign.connectRmsdTrace(this);
			}
		}
		postprocessWrite = active;
	}
	
	@Override
	public void beforeFirstSample(InputData inputData) {
//		for(ModelExtension modExt : getModExtPlugins()) {
//			if(modExt instanceof StructAlign) {
//				structAlign = (StructAlign) modExt;
//			}
//		}				
		if(!active)	return;									
		
		//double[] rad = calcGyration();
		int leaves = structAlign.coords.length;
		try {
			for(int i = 0; i < leaves-1; i++)
				for(int j = i+1; j < leaves; j++)
					outputFile.write("msd" + i + "_" + j + "\t");
			for(int i = 0; i < leaves-1; i++)
				for(int j = i+1; j < leaves; j++)
					outputFile.write("t" + i + "_" + j + "\t");
			for(int i = 0; i < leaves-1; i++)
				for(int j = i+1; j < leaves; j++)
					outputFile.write("seqID" + i + "_" + j + "\t");
			
			
//			for(int i = 0; i < rad.length; i++)
//				outputFile.write(mcmc.tree.names[i] + "\t");
//			outputFile.write("\n");
//			for(int i = 0; i < rad.length; i++)
//				outputFile.write(rad[i] + "\t"); 

			outputFile.write("\n");
		} catch (IOException e){}		
	}
	
	@Override
	public void newSample(State state, int no, int total) {
		if(!active)
			return;
		if(postprocessWrite) {
			double[][] msd = calcMSD();
			double[][] seqID = calcSeqID();			
			
			try {
				for(int i = 0; i < msd.length-1; i++)
					for(int j = i+1; j < msd.length; j++)
						outputFile.write(msd[i][j] + "\t");
				for(int i = 0; i < msd.length-1; i++)
					for(int j = i+1; j < msd.length; j++) 
						outputFile.write(distanceMatrix[i][j] + "\t");
				for(int i = 0; i < msd.length-1; i++)
					for(int j = i+1; j < msd.length; j++)
						outputFile.write(seqID[i][j] + "\t");
				outputFile.write("\n");
			} catch (IOException e){
				e.printStackTrace();
			}		
		}
	}
	
	@Override
	public void afterLastSample() {
		if(!active)
			return;		
		try {
			outputFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	public static void printMatrix(double[][] m) {
		for(int i = 0; i < m.length; i++)
			System.out.println(Arrays.toString(m[i]));
		System.out.println();
	}
	
	public double[][] calcMSD(){
		double[][][] coor = structAlign.rotCoords;
		String[] align = structAlign.curAlign;
		int leaves = coor.length;
		boolean igap, jgap;
		double[][] msd = new double[leaves][leaves];
		for(int i = 0; i < leaves-1; i++){
			for(int j = i+1; j < leaves; j++){
				int ii = 0, jj = 0, n = 0;
				for(int k = 0; k < align[0].length(); k++){
					igap = align[i].charAt(k) == '-';
					jgap = align[j].charAt(k) == '-';
					if(!igap & !jgap){
						msd[i][j] += sqDistance(coor[i][ii], coor[j][jj]);
						n++;
					}
					ii += igap ? 0 : 1;
					jj += jgap ? 0 : 1;
				}
				msd[i][j] /= n;
			}
		}
		return msd;
	}
	
	public double sqDistance(double[] x, double[] y){
		double d = 0;
		for(int i = 0; i < x.length; i++)
			d += Math.pow(x[i] - y[i], 2.0);
		return d;
	}
	
	public double[] calcGyration(){
		double[][][] coor = structAlign.coords;
		int leaves = coor.length;
		double[] radii = new double[leaves];
		for(int i = 0; i < leaves; i++){
			radii[i] = 0;
			// coordinates are centered in StructAlign.initRun()
			for(int j = 0; j < coor[i].length; j++)
				for(int k = 0; k < coor[i][0].length; k++)
					radii[i] += Math.pow(coor[i][j][k], 2.0);
			radii[i] /= coor[0].length;
			radii[i] = Math.pow(radii[i], 0.5);
		}
		return radii;
	}
	
	public double[][] calcSeqID(){
		String[] align = structAlign.curAlign;
		int leaves = align.length;
		boolean igap, jgap;
		double[][] seqID = new double[leaves][leaves];
		for(int i = 0; i < leaves-1; i++){
			for(int j = i+1; j < leaves; j++){
				double match = 0, id = 0;
				for(int k = 0; k < align[0].length(); k++){
					igap = align[i].charAt(k) == '-';
					jgap = align[j].charAt(k) == '-';
					if(!igap & !jgap){
						id += align[i].charAt(k) == align[j].charAt(k) ? 1 : 0;
						match++;
					}
				}
				seqID[i][j] = id / match;
			}
		}
		return seqID;
	}


	@Override
	public String getTabName() {
		// TODO Auto-generated method stub
		return "RMSD trace";
	}


	@Override
	public Icon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public JPanel getJPanel() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public String getTip() {
		// TODO Auto-generated method stub
		return "";
	}
}
