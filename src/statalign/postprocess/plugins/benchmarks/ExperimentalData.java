package statalign.postprocess.plugins.benchmarks;
import java.util.ArrayList;

import statalign.postprocess.utils.RNAFoldingTools;


public class ExperimentalData
{
	int [] pairedSites;
	ArrayList<String> sequences;
	ArrayList<String> sequenceNames;
	
	public String toString ()
	{
		String ret = "";
		ret += RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites) + "\n";
		for(int i = 0 ; i < Math.min(sequences.size(), sequenceNames.size()) ; i++)
		{
			ret += ">" + sequenceNames.get(i) + "\n";
			ret += sequences.get(i) + "\n";
			
		}
		
		return ret;
	}
}
