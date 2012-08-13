package statalign.postprocess.plugins;

import java.io.Serializable;

import statalign.postprocess.utils.RNAFoldingTools;


public class RNAalifoldResult implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 957877165963191430L;
	
	public double [][] matrix;
	public int [] pairedSites;
	public double reliablityScore;
	public double pairsOnlyReliablityScore;
	
	public RNAalifoldResult getSmallResult()
	{
		RNAalifoldResult result = new RNAalifoldResult();
		result.pairedSites = this.pairedSites;
		result.reliablityScore = RNAFoldingTools.calculatePPfoldReliabilityScore(pairedSites, matrix);
		result.pairsOnlyReliablityScore = RNAFoldingTools.calculatePairsOnlyReliabilityScore(pairedSites, matrix);
		return result;
	}
}
