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
	
	public double[] getAcceptanceRates() {
		int n = structAlign.proposalCounts.length;
		double[] acceptanceRates = new double[n];
		for (int i=0; i<n; i++) {
			acceptanceRates[i] = (double) structAlign.acceptanceCounts[i] / (double) structAlign.proposalCounts[i]; 
		}
		return acceptanceRates;
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
		return "Structural alignment";
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
			int sigLen = structAlign.globalSigma ? 1 : structAlign.sigma2.length;
			outputFile.write("sigma2_1");
			for(int i = 1; i < sigLen; i++)
				outputFile.write("\tsigma2_"+(i+1));
			outputFile.write("\ttau");
			outputFile.write("\tepsilon");
			if (!structAlign.globalSigma) {
				outputFile.write("\tsigma2H");
				outputFile.write("\tnu");
			}
			for(int i = 0; i < sigLen; i++)
				outputFile.write("\tsigma2_"+(i+1)+"_proposed");
			outputFile.write("\ttau_proposed");
			outputFile.write("\tepsilon_proposed");
			if (!structAlign.globalSigma) {
				outputFile.write("\tsigma2H_proposed");
				outputFile.write("\tnu_proposed");
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
	
	int[] lastSigmaProp;
	int lastTauProp;
	int lastEpsilonProp;
	int lastSigma2HProp;
	int lastNuProp;
	
	int[] sigma2Proposed;
	int tauProposed;
	int epsilonProposed;
	int sigma2HProposed;
	int nuProposed;
	
	@Override
	public void newSample(State state, int no, int total) {
		if(!active)
			return;
		if(postprocessWrite) {
			try {
				int sigLen = (structAlign.globalSigma ? 1 : structAlign.sigma2.length);
				for(int i = 0; i < sigLen; i++)
					outputFile.write(structAlign.sigma2[i]+"\t");
				outputFile.write(structAlign.tau+"\t");
				outputFile.write(structAlign.epsilon+"\t");
				if (!structAlign.globalSigma) {
					outputFile.write(structAlign.sigma2Hier+"\t");
					outputFile.write(structAlign.nu+"\t");
				}
				
				if(lastSigmaProp == null || lastSigmaProp.length != sigLen)
					lastSigmaProp = new int[sigLen];
				int i=0;
				for(i = 0; i < sigLen; i++) {
					outputFile.write(lastSigmaProp[i] != structAlign.proposalCounts[i] ? structAlign.sigma2[i]+"\t" : -1+"\t");
					lastSigmaProp[i] = structAlign.proposalCounts[i];
				}
				outputFile.write(lastTauProp != structAlign.proposalCounts[i] ? structAlign.tau+"\t" : -1+"\t");
				lastTauProp = structAlign.proposalCounts[i];
				++i;
				outputFile.write(lastEpsilonProp != structAlign.proposalCounts[i] ? structAlign.epsilon+"\t" : -1+"\t");
				lastEpsilonProp = structAlign.proposalCounts[i];
				if (!structAlign.globalSigma) {
					++i;
					outputFile.write(lastSigma2HProp != structAlign.proposalCounts[i] ? structAlign.sigma2Hier+"\t" : -1+"\t");
					lastSigma2HProp = structAlign.proposalCounts[i];
					++i;
					outputFile.write(lastNuProp != structAlign.proposalCounts[i] ? structAlign.nu+"\t" : -1+"\t");
					lastNuProp = structAlign.proposalCounts[i];
				}

				outputFile.write("\n");
			} catch (IOException e) {
				e.printStackTrace(); 
			}
		}
		
//		if(sampling) {
//			try {
//				file.write("Sample "+no+"\tStructure:\t");
//			} catch (IOException e) {
//			}
//		}
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

			System.out.print("Sigma2 (prop): ");
			for (int i=0; i<structAlign.proposalCounts.length; i++) {
				System.out.print(structAlign.proposalCounts[i]+" ");
			}
			System.out.println("");
			System.out.print("Sigma2 (acce): ");
			for (int i=0; i<structAlign.acceptanceCounts.length; i++) {
				System.out.print(structAlign.acceptanceCounts[i]+" ");
			}
			System.out.println("Rotation: " + structAlign.rotProposed + " " + structAlign.rotAccept);
			System.out.println("Xlat: " + structAlign.xlatProposed + " " + structAlign.xlatAccept);
			System.out.println("Library: " + structAlign.libProposed + " " + structAlign.libAccept);
		}
		
		System.out.println("final translations:");
		for(int i = 0; i < structAlign.xlats.length; i++) {
			System.out.println(Arrays.toString(structAlign.xlats[i]));
		}
		
		System.out.println();
		System.out.println("Acceptance rates:");
		//System.out.println("Sigma2: " + structAlign.sigProposed + " " + structAlign.sigAccept);
		System.out.println("Rotation: " + structAlign.rotProposed + " " + structAlign.rotAccept);
		System.out.println("Xlat: " + structAlign.xlatProposed + " " + structAlign.xlatAccept);
		System.out.println("Library: " + structAlign.libProposed + " " + structAlign.libAccept);
		//System.out.println("Tau: " + structAlign.tau + "  Epsilon: " + structAlign.epsilon);
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
		if (count % refreshRate == 0) {
			StructAlignTraceParameters currentParameters = new StructAlignTraceParameters(mcmcStep.burnIn);

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
