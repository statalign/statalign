package statalign.model.ext.plugins;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import statalign.base.Tree;
import statalign.model.ext.ModelExtension;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.utils.RNAalifold;

public class PPFold extends ModelExtension {

	private final String pluginID = "ppfold";	
	
	@Override
	public String getPluginID() {
		return pluginID;
	}
	
	@Override
	public double logLikeFactor(Tree tree) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setActive(boolean active) {
		super.setActive(active);
		PostprocessManager.pluginParameters.setParameter("ppfold", "");
		PostprocessManager.rnaMode = true;
		System.out.println("PPFold plugin is now "+(active?"enabled":"disabled"));		
	}
	@Override
	public void setParam(String paramName, boolean paramValue) {
		if (paramName.equals("experimental")) {
			PostprocessManager.pluginParameters.setParameter("experimental", "");
		}
	}	
	
	@Override
	public String getUsageInfo() {
		StringBuilder usage = new StringBuilder();
		usage.append("____________________\n\n");
		usage.append("  PPFold plugin \n\n");
		usage.append("^^^^^^^^^^^^^^^^^^^^\n\n");
		usage.append("java -jar statalign.jar -plugin:ppfold[[OPTION1,OPTION2,...]]\n");
		usage.append("OPTIONS: \n");		
		usage.append("\texperimental\t\t(Activates experimental features. See online tutorial for details.)\n");
		usage.append("\nNote that the above syntax is designed to work in bash shells. " +
				"Other shells such as csh may require square brackets to be preceded by a backslash.");
												
		return usage.toString();
	}
}
