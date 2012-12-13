package statalign.postprocess.plugins.structalign;

import java.io.IOException;
import java.util.Arrays;

import javax.swing.Icon;
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

public class StructTrace extends Postprocess {
	
	StructAlign structAlign;

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
	public Icon getIcon() {
		return null;
	}

	@Override
	public JPanel getJPanel() {
		return null;
	}

	@Override
	public String getTip() {
		return "Structural protein alignment";
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
		try {
			outputFile.write("tau");
			int sigLen = structAlign.globalSigma ? 1 : 2*inputData.seqs.sequences.size()-1;
			for(int i = 0; i < sigLen; i++)
				outputFile.write("\tsigma2_"+(i+1));
			outputFile.write("\ttau_proposed");
			for(int i = 0; i < sigLen; i++)
				outputFile.write("\tsigma2_"+(i+1)+"_proposed");
			outputFile.write("\n");
		} catch (IOException e) {
		}
		lastSigmaProp = null;
		lastTauProp = 0;
	}
	
	int[] lastSigmaProp;
	int lastTauProp;
	
	@Override
	public void newSample(State state, int no, int total) {
		if(postprocessWrite) {
			try {
				int sigLen = structAlign.globalSigma ? 1 : structAlign.sigma2.length;
				outputFile.write(structAlign.tau+"\t");
				for(int i = 0; i < sigLen; i++)
					outputFile.write(structAlign.sigma2[i]+"\t");
				
				if(lastSigmaProp == null || lastSigmaProp.length != sigLen)
					lastSigmaProp = new int[sigLen];
				outputFile.write(lastTauProp != structAlign.tauProposed ? ""+structAlign.tau : "");
				lastTauProp = structAlign.tauProposed;
				for(int i = 0; i < sigLen; i++) {
					outputFile.write(lastSigmaProp[i] != structAlign.proposalCounts[i] ? "\t"+structAlign.sigma2[i] : "\t");
					lastSigmaProp[i] = structAlign.proposalCounts[i];
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
	}
	
	public static void printMatrix(double[][] m) {
		for(int i = 0; i < m.length; i++)
			System.out.println(Arrays.toString(m[i]));
		System.out.println();
	}
	
	@Override
	public void newStep(McmcStep mcmcStep) {
	}
	
}
