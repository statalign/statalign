package statalign.model.ext;

import java.util.List;

/**
 * Manager for {@link ModelExtension} plugins. Provides utility methods for the plugins.
 * 
 * @author novak
 */
public class ModelExtManager {
	
	public ModelExtInterface parent;
	
	public ModelExtManager(ModelExtInterface parent) {
		this.parent = parent;
	}
	
	/**
	 * Returns the (unmodifiable) list of recognised {@link ModelExtension} plugins.
	 */
	public List<ModelExtension> getPluginList() {
		return parent.getPluginList();
	}
}
