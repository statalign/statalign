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
import statalign.base.Utils;
import statalign.io.input.plugins.PDBReader;
import statalign.mcmc.McmcMove;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.plugins.StructAlign;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.StructAlignTraceGUI;


public class StructTrace extends Postprocess {

	List<StructAlignTraceParameters> parameterHistory;

	public StructAlign structAlign;
	
	public int burninLength; 
	public int MAX_HISTORY_SIZE = 1000;
	public int refreshRate;
	
	double maxLikelihood = Double.NEGATIVE_INFINITY;
	int sampleNumberMLE;
	double[][][] coorMLE;
	String[] seqs, names;
	
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
		postprocessWrite = false;
		selected = false;
		active = false;
	}

	@Override
	public void init(ModelExtManager modelExtMan) {
		for(ModelExtension modExt : modelExtMan.getPluginList()) {
			if(modExt instanceof StructAlign) {
				structAlign = (StructAlign) modExt;
				structAlign.connectStructTrace(this);
			}
		}
		active = structAlign.isActive();
		postprocessWrite = active;
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
		return "struct.params";
	}

	@Override
	public void setSampling(boolean enabled) {
	}
	

	
	@Override
	public void beforeFirstSample(InputData inputData) {
//		for(ModelExtension modExt : getModExtPlugins()) {
//			if(modExt instanceof StructAlign) {
//				structAlign = (StructAlign) modExt;
//			}
//		}
		
		if(!active)
			return;
		try {
			for (McmcMove mcmcMove : structAlign.getMcmcMoves()) {
				if (mcmcMove.getParam() != null) {
					outputFile.write(mcmcMove.name);
					// print a tab separator after all but the last column
					if (mcmcMove != structAlign.getMcmcMoves().get(structAlign.getMcmcMoves().size()-1)) {
						outputFile.write("\t");
					}
					//outputFile.write(mcmcMove.name+" (Proposed)\t");
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
		MAX_HISTORY_SIZE = Math.min(burninLength, MAX_HISTORY_SIZE);
		current = 0;
		//refreshRate = inputData.pars.burnIn / (2*MAX_HISTORY_SIZE);
		refreshRate = burninLength / (MAX_HISTORY_SIZE);
		// Means we will have the whole burnin in one window, but then it
		// will start to shift.
		count = 0;
		
	}
	private void doUpdate(State state, int sampleNumber) {				
		if (!state.isBurnin && state.logLike > maxLikelihood) {
			maxLikelihood = state.logLike;			
			sampleNumberMLE = sampleNumber;
			coorMLE = structAlign.rotCoords.clone();
			if (seqs == null) { 
				seqs = state.seq; 
				names = state.name;
			}
		}
	}
	@Override
	public void newPeek(State state) {
		if (!active) return;
		//if (show) {
			doUpdate(state,0);
		//}
	}
	@Override
	public void newSample(State state, int no, int total) {
		if(!active)
			return;
		if(postprocessWrite) {
			try {
				for (McmcMove mcmcMove : structAlign.getMcmcMoves()) {
					if (mcmcMove.getParam() != null) {
						outputFile.write(mcmcMove.getParam().get()+"");
						// print a tab separator after all but the last column
						if (mcmcMove != structAlign.getMcmcMoves().get(structAlign.getMcmcMoves().size()-1)) {
							outputFile.write("\t");
						}
//						if (mcmcMove.moveProposed) {
//							outputFile.write(mcmcMove.getParam().get()+"\t");
//						}
//						else {
//							outputFile.write(-1+"\t");
//						}
					}
				}				
				outputFile.write("\n");
			} catch (IOException e) {
				e.printStackTrace(); 
			}
			//structAlign.setAllMovesNotProposed();
		}				
		doUpdate(state,no);
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
		if (coorMLE != null) {
			try {
				FileWriter mle = new FileWriter(getBaseFileName()+"mle.super.pdb");
				PDBReader.writePDB(coorMLE, seqs, names, mle);				
				mle.close();
			}catch (IOException e) {
				e.printStackTrace();
			}
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
}
