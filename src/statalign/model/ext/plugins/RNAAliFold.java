package statalign.model.ext.plugins;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import statalign.base.Tree;
import statalign.model.ext.ModelExtension;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.utils.RNAalifold;

public class RNAAliFold extends ModelExtension {

	private final String pluginID = "rnaalifold";	
	
	String params = "";
		
	@Override
	public String getPluginID() {
		return pluginID;
	}
	
	@Override
	public double logLikeFactor(Tree tree) {
		return 0;
	}

	@Override
	public void setActive(boolean active) {
		super.setActive(active);
		PostprocessManager.pluginParameters.setParameter("rnaalifold", "");
		PostprocessManager.rnaMode = true;
		System.out.println("RNAalifold plugin is now "+(active?"enabled":"disabled"));		
	}
	@Override
	public void setParam(String paramName, boolean paramValue) {
		if (paramName.equals("experimental")) {
			PostprocessManager.pluginParameters.setParameter("experimental", "");
		}
		else if (paramName.equals("circ")) {
			params += " -circ";
		}
	}	
	
	public void setParam(String paramName, String paramValue) {
		
		
		if (paramName.equals("exe")) {
			RNAalifold.executable = paramValue;
			System.out.println("Setting rnaalifold executable path to "+paramValue);
		}		
		else if (paramName.equals("T")) {
			params += " -T "+paramValue;
		}
		else if (paramName.equals("cov")) {
			params += " --cfactor "+paramValue;
		}
		else if (paramName.equals("n")) {
			params += " --nfactor "+paramValue;
		}
	}
	
	@Override
	public void init() {
		params = RNAalifold.executable + params;
		PostprocessManager.pluginParameters.setParameter("rnaalifold",params);
	}
	
	@Override
	public String getUsageInfo() {
		StringBuilder usage = new StringBuilder();
		usage.append("________________________\n\n");
		usage.append("  RNAalifold plugin \n\n");
		usage.append("^^^^^^^^^^^^^^^^^^^^^^^^\n\n");
		usage.append("java -jar statalign.jar -plugin:rnaalifold[OPTION1,OPTION2,...]\n");
		usage.append("OPTIONS: \n");		
		usage.append("\tcirc\t\t(Activates circular mode.)\n");
		usage.append("\texperimental\t\t(Activates experimental features.)\n");
		usage.append("\tT=TEMP\t\t(Sets the temperature.)\n");
		usage.append("\tcov=COV\t\t(Sets the covariance factor.)\n");
		usage.append("\tn=NFACTOR\t\t(Sets the non-compatibility factor.)\n");
		usage.append("\texe=EXE\t\t(Sets the path to the RNAalifold executable.)\n");
		usage.append("\t\t\tThis should generally be set to \"/path/to/statalign/lib/RNAalifold.exe\"\n");
		usage.append("\nSee online tutorial for further details.");
		usage.append("\nNote that the above syntax is designed to work in bash shells. " +
				"Other shells such as csh may require square brackets to be preceded by a backslash.");
												
		return usage.toString();
	}

}
