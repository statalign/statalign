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
			outputFile.write("sigma2\ttheta\tsigma2_proposed\ttheta_proposed\n");
		} catch (IOException e) {
		}
	}
	
	int lastSigmaProp;
	int lastThetaProp;
	
	@Override
	public void newSample(State state, int no, int total) {
		if(postprocessWrite) {
			// TODO decide if this is still needed for new parameterization
			/** try {
				outputFile.write(structAlign.sigma2+"\t"+structAlign.theta+"\t");
				int newSigmaProp = structAlign.sigProposed;
				outputFile.write(lastSigmaProp != newSigmaProp? ""+structAlign.sigma2 : "");
				outputFile.write("\t");
				lastSigmaProp = newSigmaProp;
				int newThetaProp = structAlign.thetaProposed;
				outputFile.write(lastThetaProp != newThetaProp ? ""+structAlign.theta : "");
				outputFile.write("\n");
				lastThetaProp = newThetaProp; 
			} catch (IOException e) {
				e.printStackTrace(); 
			} **/
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
