package statalign.postprocess.plugins.structalign;

import java.awt.BorderLayout;
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
import statalign.model.ext.ModelExtension;
import statalign.model.ext.McmcMove;
import statalign.model.ext.plugins.StructAlign;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.StructAlignTraceGUI;
import statalign.postprocess.utils.StructAlignTraceParameters;


public class StructTrace extends Postprocess {
	
	StructAlign structAlign;
	List<StructAlignTraceParameters> parameterHistory;
	
	public int burninLength; 
	public final int MAX_HISTORY_SIZE = 1000;
	public int refreshRate;
	
	public List<StructAlignTraceParameters> getParameterHistory() {
		return parameterHistory;
	}
	
	JPanel pan;
	int current;
	int count;
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
	public void beforeFirstSample(InputData inputData) {
		for(ModelExtension modExt : getModExtPlugins()) {
			if(modExt instanceof StructAlign) {
				structAlign = (StructAlign) modExt;
			}
		}
		active = structAlign.isActive();
		if(!active)
			return;
		try {
			for (McmcMove mcmcMove : structAlign.getMcmcMoves()) {
				if (mcmcMove.getParam() != null) {
					outputFile.write(mcmcMove.name+"\t");
				}
			}				
			outputFile.write("\n");
		} catch (IOException e) {
		}
		lastSigmaProp = null;
		lastTauProp = 0;
		lastEpsilonProp = 0;
		lastSigma2HProp = 0;
		lastNuProp = 0;
		
		sigma2Proposed = null;
		tauProposed = 0;
		epsilonProposed = 0;
		sigma2HProposed = 0;
		nuProposed = 0;
		
		if(show) {
			pan.removeAll();
			gui = new StructAlignTraceGUI(pan, this);
			pan.add(gui);
			pan.getParent().getParent().getParent().validate();
		}
		
		parameterHistory = new ArrayList<StructAlignTraceParameters>();
		burninLength = inputData.pars.burnIn;
		current = 0;
		refreshRate = inputData.pars.burnIn / (2*MAX_HISTORY_SIZE);
		count = 0;
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
						if (mcmcMove.lastMoveAccepted) {
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
		}
		System.out.println("Acceptance rates:");
		for (McmcMove mcmcMove : structAlign.getMcmcMoves()) {
			System.out.println(mcmcMove.name+"\t"+mcmcMove.acceptanceRate());
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
		++count;
		if (count > burninLength / 2) {
			structAlign.MIN_EPSILON = 0.1;
		}
		if (count % refreshRate == 0) {
			StructAlignTraceParameters currentParameters = new StructAlignTraceParameters(mcmcStep.burnIn);
			/*
			Parameters currentParameters = structAlign.getParameters().clone(); 
			*/
			currentParameters.globalSigma = structAlign.globalSigma;

			currentParameters.tau = structAlign.tau;
			if (structAlign.proposalCounts[structAlign.tauInd] != tauProposed) {
				currentParameters.tauProposed = true;
				tauProposed = structAlign.proposalCounts[structAlign.tauInd];
			}
			else {
				currentParameters.tauProposed = false;
			}
			currentParameters.epsilon = structAlign.epsilon;
			if (structAlign.proposalCounts[structAlign.epsilonInd] != epsilonProposed) {
				currentParameters.epsilonProposed = true;
				epsilonProposed = structAlign.proposalCounts[structAlign.epsilonInd];
			}
			else {
				currentParameters.epsilonProposed = false;
			}
			if (!currentParameters.globalSigma) {
				currentParameters.nu = structAlign.nu;
				if (structAlign.proposalCounts[structAlign.nuInd] != nuProposed) {
					currentParameters.nuProposed = true;
					nuProposed = structAlign.proposalCounts[structAlign.nuInd];
				}
				else {
					currentParameters.nuProposed = false;
				}
				currentParameters.sigma2Hier = structAlign.sigma2Hier;
				if (structAlign.proposalCounts[structAlign.sigma2HInd] != sigma2HProposed) {
					currentParameters.sigma2HProposed = true;
					sigma2HProposed = structAlign.proposalCounts[structAlign.sigma2HInd];
				}
				else {
					currentParameters.sigma2HProposed = false;
				}
			}
			currentParameters.sigma2 = structAlign.sigma2.clone();
			int sigLen = currentParameters.sigma2.length;
			currentParameters.sigma2Proposed = new boolean[sigLen];
			if(sigma2Proposed == null || sigma2Proposed.length != sigLen) {
				sigma2Proposed = new int[sigLen];
				for (int i=0; i<sigLen; i++) {
					sigma2Proposed[i] = 0;
				}
			}	
			for (int i=0; i<sigLen; i++) {
				if (structAlign.proposalCounts[i] != sigma2Proposed[i]) {
					currentParameters.sigma2Proposed[i] = true;
					sigma2Proposed[i] = structAlign.proposalCounts[i];
				}
				else {
					currentParameters.sigma2Proposed[i] = false;
				}
			}
			
								
			if(parameterHistory.size() < MAX_HISTORY_SIZE){
				parameterHistory.add(currentParameters);
			} else {
				parameterHistory.remove(0);
				parameterHistory.add(currentParameters);				
			}
			if(show) {
				gui.repaint();
			}
		}
	}
	
}
