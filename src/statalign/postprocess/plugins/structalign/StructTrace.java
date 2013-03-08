package statalign.postprocess.plugins.structalign;

import java.awt.BorderLayout;
import java.io.IOException;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import statalign.base.InputData;
import statalign.base.Mcmc;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.base.Utils;
import statalign.mcmc.McmcMove;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.plugins.StructAlign;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.StructAlignTraceGUI;


public class StructTrace extends Postprocess {

	FileWriter rmsdOut;
	FileWriter radiiOut;
	
	public StructAlign structAlign;
	List<StructAlignTraceParameters> parameterHistory;
	
	public int burninLength; 
	public final int MAX_HISTORY_SIZE = 1000;
	public int refreshRate;
	
	public List<StructAlignTraceParameters> getParameterHistory() {
		return parameterHistory;
	}
	
	JPanel pan;
	int current;
	private int count;
	public int getCount() {
		return count;
	}
	private StructAlignTraceGUI gui;

	public StructTrace() {
		screenable = true;
		outputable = true;
		postprocessable = true;
		postprocessWrite = true;
		selected = true;
		active = true;
	}

	@Override
	public String getTabName() {
		return "Structural parameters";
	}

	@Override
	public double getTabOrder() {
		return 5.0d;
	}
	
	@Override
	public Icon getIcon() {
		return new ImageIcon(ClassLoader.getSystemResource("icons/loglikelihood1.gif"));
	}

	@Override
	public JPanel getJPanel() {
		pan = new JPanel(new BorderLayout());
		return pan;
	}

	@Override
	public String getTip() {
		return "StructAlign parameter values";
	}
	
	@Override
	public String getFileExtension() {
		return "struct";
	}

	@Override
	public void setSampling(boolean enabled) {
	}
	
	@Override
	public void init(ModelExtManager modelExtMan) {
		for(ModelExtension modExt : modelExtMan.getPluginList()) {
			if(modExt instanceof StructAlign) {
				structAlign = (StructAlign) modExt;
			}
		}
	}
	
	@Override
	public void beforeFirstSample(InputData inputData) {
//		for(ModelExtension modExt : getModExtPlugins()) {
//			if(modExt instanceof StructAlign) {
//				structAlign = (StructAlign) modExt;
//			}
//		}
		active = structAlign.isActive();
		if(!active)
			return;
		try {
			for (McmcMove mcmcMove : structAlign.getMcmcMoves()) {
				if (mcmcMove.getParam() != null) {
					outputFile.write(mcmcMove.name+"\t");
					outputFile.write(mcmcMove.name+" (Proposed)\t");
				}
			}				
			outputFile.write("\n");
		} catch (IOException e) {
		}
		
		if(show) {
			pan.removeAll();
			gui = new StructAlignTraceGUI(pan, this);
			pan.add(gui);
			pan.getParent().getParent().getParent().validate();
		}
		
		parameterHistory = new ArrayList<StructAlignTraceParameters>();
		burninLength = inputData.pars.burnIn;
		current = 0;
		//refreshRate = inputData.pars.burnIn / (2*MAX_HISTORY_SIZE);
		refreshRate = inputData.pars.burnIn / (MAX_HISTORY_SIZE);
		// Means we will have the whole burnin in one window, but then it
		// will start to shift.
		count = 0;
		
		try{
			rmsdOut = new FileWriter("rmsd.txt");
		} catch (IOException e){}
		try{
			radiiOut = new FileWriter("radii.txt");
		} catch (IOException e){}
		
		double[] rad = calcGyration();
		for(int i = 0; i < rad.length; i++){
			try {radiiOut.write(rad[i] + "\n");
			} catch (IOException e){}
		}
	}
	
	@Override
	public void newSample(State state, int no, int total) {
		if(!active)
			return;
		if(postprocessWrite) {
			try {
				for (McmcMove mcmcMove : structAlign.getMcmcMoves()) {
					if (mcmcMove.getParam() != null) {
						outputFile.write(mcmcMove.getParam().get()+"\t");
						if (mcmcMove.moveProposed) {
							outputFile.write(mcmcMove.getParam().get()+"\t");
						}
						else {
							outputFile.write(-1+"\t");
						}
					}
				}				
				outputFile.write("\n");
			} catch (IOException e) {
				e.printStackTrace(); 
			}
			//structAlign.setAllMovesNotProposed();
			
			double[][] msd = calcMSD();
			
			try {
				for(int i = 0; i < msd.length-1; i++)
					for(int j = i+1; j < msd.length; j++)
						rmsdOut.write(msd[i][j] + "\t" + structAlign.distanceMatrix[i][j] + "\t");
				rmsdOut.write("\n");
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
		if(Utils.DEBUG) {
			System.out.println("final rotation matrices:");
			for(int i = 1; i < structAlign.xlats.length; i++) {
				Rotation rot = new Rotation(new Vector3D(structAlign.axes[i]), structAlign.angles[i]);
				printMatrix(rot.getMatrix());
			}
			System.out.println("final translations:");
			for(int i = 0; i < structAlign.xlats.length; i++) {
				System.out.println(Arrays.toString(structAlign.xlats[i]));
			}
			System.out.println();				
			System.out.println("Acceptance rates:");
			for (McmcMove mcmcMove : structAlign.getMcmcMoves()) {
				System.out.println(mcmcMove.name+"\t"+mcmcMove.acceptanceRate());
			}
		}
	}
	
	public static void printMatrix(double[][] m) {
		for(int i = 0; i < m.length; i++)
			System.out.println(Arrays.toString(m[i]));
		System.out.println();
	}
	
	@Override
	public void newStep(McmcStep mcmcStep) {
		if(!active)
			return;
		if (screenable && (count % refreshRate == 0)) {
			StructAlignTraceParameters currentParameters = 
				new StructAlignTraceParameters(this,mcmcStep.burnIn);
			
			currentParameters.globalSigma = structAlign.globalSigma;		
			if (count > 0) {
				currentParameters.setProposalFlags(parameterHistory.get(parameterHistory.size()-1));
			}
			if(parameterHistory.size() <= MAX_HISTORY_SIZE){
				parameterHistory.add(currentParameters);
			} else {
				parameterHistory.remove(0);
				parameterHistory.add(currentParameters);				
			}
			if(show) {
				gui.repaint();
			}
		}
		++count;
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
}
