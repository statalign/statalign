package statalign.postprocess.plugins.structalign;

import java.io.FileWriter;
import java.io.IOException;

import javax.swing.Icon;
import javax.swing.JPanel;

import statalign.base.InputData;
import statalign.base.McmcStep;
import statalign.base.State;
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
		return ".struct";
	}

	@Override
	public void setSampling(boolean enabled) {
	}
	
	FileWriter sigmaWrite;
	FileWriter thetaWrite;
	
	@Override
	public void beforeFirstSample(InputData inputData) {
		for(ModelExtension modExt : getModExtPlugins()) {
			if(modExt instanceof StructAlign) {
				structAlign = (StructAlign) modExt;
			}
		}
		try {
			sigmaWrite = new FileWriter("sigma2");
			thetaWrite = new FileWriter("theta");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	@Override
	public void newSample(State state, int no, int total) {
		try {
			sigmaWrite.write(structAlign.sigma2+"\n");
			thetaWrite.write(structAlign.theta+"\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	@Override
	public void afterLastSample() {
		try {
			outputFile.close();
			sigmaWrite.close();
			thetaWrite.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	@Override
	public void newStep(McmcStep mcmcStep) {
	}
	
}
