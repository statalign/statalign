package statalign.postprocess;

import java.util.ArrayList;
import java.util.Hashtable;

public class PluginParameters {
	
	Hashtable<String, String> parameters = new Hashtable<String, String>();
	
	public PluginParameters()
	{
		
	}
		
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
	
	public String getParameter(String name)
	{
		if(parameters.containsKey(name))
		{
			return parameters.get(name);
		}
		return null;
	}
	
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

}
