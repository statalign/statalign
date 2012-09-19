package statalign.postprocess;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * A class which for storing plugin parameters.
 * 
 * Any parameters starting with "plugin:" will be passed from the CommandLine class to 
 * a static PluginParameters object which is visible to all Postprocess classes. It is
 * the responsibility of the individual Postprocess classes to retrieve the parameters specific to that
 * they require.
 */
public class PluginParameters {
	
	Hashtable<String, String> parameters = new Hashtable<String, String>();
	
	public PluginParameters()
	{
		
	}
		
	/**
	 * An alternative way to initialise this class.
	 * @param args a list of parameters and values of the form "parameter=value".
	 */
	public PluginParameters(ArrayList<String> args)
	{
		for(int i = 0 ; i < args.size() ; i++)
		{
			String [] split = args.get(i).split("=", 2);
			String param = split[0];
			if(split.length == 1)
			{
				parameters.put(param, "");
			}
			else
			{
				parameters.put(param, split[1]);
			}
		}
	}
	
	/**
	 * Returns the corresponding parameter value or null if the parameter does not exist.
	 * @param name the name of the parameter to retrieve, excluding the "plugin:" suffix.
	 * @return the corresponding value.
	 */
	public String getParameter(String name)
	{
		if(parameters.containsKey(name))
		{
			return parameters.get(name);
		}
		return null;
	}
	
	/**
	 * Given a parameter and a value, sets the corresponding parameter value.
	 * @param name the name of the parameter to set, excluding the "plugin:" suffix.
	 * @param value the value of the parameter.
	 */
	public void setParameter(String name, String value)
	{
		parameters.put(name, value);
	}
	
	public void removeParameter(String name)
	{
		if(parameters.containsKey(name))
		{
			parameters.remove(name);
		}
	}
	
	public void print()
	{
		Enumeration<String> keys = parameters.keys();
		String key = null;
		while(keys.hasMoreElements())
		{
			key = keys.nextElement();
			System.out.println(key+":"+parameters.get(key));
		}
	}

}
