package statalign.postprocess.plugins.benchmarks;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;

import statalign.postprocess.plugins.RNAalifoldResult;
import statalign.postprocess.utils.RNAFoldingTools;

import com.ppfold.algo.AlignmentData;
import com.ppfold.algo.FuzzyAlignment;
import com.ppfold.algo.ResultBundle;

public class Dataset implements Serializable {
	
	public String title;
	public long randomSeed = 1;
	public int burnIn = 0;
	public int mcmcSteps = 0;
	public int samplingRate = 0;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3205713573185024480L;
	
	public File fastaFile;
	
	public AlignmentData inputAlignment = new AlignmentData();
	public AlignmentData mpdAlignment = new AlignmentData();
	
	/*
	public double reliabilityOnSamplingAndAveraging = -1;
	public double reliabilityOnReference = -1;
	public double reliabilityOnMPD = -1;
	public double reliabilityOnEntropyExp = -1;
	public double reliabilityOnEntropyObs = -1;
	*/
	
	public ArrayList<Double> posteriors = new ArrayList<Double>();
	public double posteriorsAverage = -1;
	public double mpdVsInputSim = -1;
	
	public ArrayList<AlignmentData> sampledAlignments = new ArrayList<AlignmentData>();
	public ArrayList<Double> logLikelihoods = new ArrayList<Double>();
	public ArrayList<ResultBundle> sampledStructures = new ArrayList<ResultBundle>();	
	public ArrayList<int[]> pairedSitesProjectedSamples = new ArrayList<int[]>();
	public ArrayList<int[]> projectPairedSitesOnCumulativeMatrix = new ArrayList<int[]>();
	
	public ArrayList<FuzzyAlignment> cumulativeFuzzyAlignment = new ArrayList<FuzzyAlignment>();
	public ArrayList<ResultBundle> cumulativeFuzzyExpResults = new ArrayList<ResultBundle>(); 
	public ArrayList<ResultBundle> cumulativeFuzzyObsResults = new ArrayList<ResultBundle>();
	//public ArrayList<Double> cumulativeEntropyExp = new ArrayList<Double>();
	//public ArrayList<Double> cumulativeEntropyObs = new ArrayList<Double>();
	//public double entropyExpSamplignAndAveraging;
	//public double entropyObsSamplignAndAveraging;
	public ResultBundle resultBundleReference;
	public ResultBundle resultBundleMPD;
	public ResultBundle resultBundleEntropyExp;
	public ResultBundle resultBundleEntropyObs;
	public ResultBundle resultBundlePPfold;
	
	public double fuzzyAlignmentObsEntropyVal= 0;
	public double fuzzyAlignmentObsEntropyPerc= 0;
	public double fuzzyAlignmentObsEntropyMax= 0;
	public double ppfoldReliabilityEntropyObs;
	public double pairsOnlyReliabilityEntropyObs;
	public double pairsOnlyReliabilityEntropyObsPosteriorWeighted = -1;
	
	public double ppfoldReliabilityScoreSamplingAndAveraging = -1;
	public double ppfoldReliabilityScoreSamplingAndAveragingWeighted = -1;
	public double pairsOnlyReliabilityScoreSamplingAndAveragingPosteriorWeighted = -1;
	public double ppfoldReliabilityMPD = -1;
	
	public double pairsOnlyReliabilityScoreSamplingAndAveraging = -1;
	public double pairsOnlyReliabilityScoreSamplingAndAveragingWeighted = -1;
	public double pairsOnlyReliabilityMPD = -1;
	public double pairsOnlyMPDPosteriorWeighted = -1;
	
	
	public int [] pairedSites;
	public int [] pairedSitesWeighted;	
	public int [] pairedSitesMPD; // also projected
	public int [] pairedSitesPPfoldProjected;
	
	
	public String pairedSitesRefSeq;
	public String pairedSitesWeightedRefSeq;
	

	public int [] pairedSitesEntropyExp;
	public int [] pairedSitesEntropyObs;
	
	//public ArrayList<RNAalifoldResult> rnaAlifoldSamples = new ArrayList<RNAalifoldResult>();
	public ArrayList<int[]> pairedSitesProjectedRnaAlifoldSamples = new ArrayList<int[]>();
	//public RNAalifoldResult rnaAlifoldFinal;
	public int [] pairedSitesRNAalifold;
	public double ppfoldReliabilityScoreRNAalifold;
	public double pairsOnlyReliabilityScoreRNAalifold;
	
	public RNAalifoldResult rnaAlifoldMPD;
	//public int[] pairedSitesMPDprojected;
	
	public RNAalifoldResult rnaAlifoldRef;
	public int[] pairedSitesRNAalifoldMPDProjected;
	public int[] pairedSitesRNAalifoldRefProjected;
	//public int [] pairedSitesMPD;
	
	/*public int [] samplingAndAveragingPairedSites;
	public int [] refencedPairedSites;
	public int [] mpdPairedSites;*/
	public int [] pairedSitesCombined;
	public double ppfoldReliabilityScoreCombined;
	public double pairsOnlyReliabilityScoreCombined;
	
	public ArrayList<ResultBundle> matrixFolds = new ArrayList<ResultBundle>();
	
	
	public void saveDatasetResult(File outFile)
	{
		try
		{
			ObjectOutput out = new ObjectOutputStream(new FileOutputStream(outFile));
		    out.writeObject(this);
		    out.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	public static Dataset loadDatasetResult(File inFile)
	{
		try
		{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
			Dataset ret = (Dataset) in.readObject();
			in.close();
			return ret;
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
}
