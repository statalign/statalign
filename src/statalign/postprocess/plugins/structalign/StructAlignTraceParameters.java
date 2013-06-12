package statalign.postprocess.plugins.structalign;

import java.util.List;
import java.util.ArrayList;

import statalign.mcmc.McmcMove;
import statalign.model.ext.plugins.structalign.ContinuousPositiveStructAlignMove;
import statalign.model.ext.plugins.StructAlign;

/**
 * A container class that stores a parameter values from the StructAlign plugin
 * and a logical value indicating whether the sample comes from the burn-in period or not.
 * 
 * @author herman
 *
 */
public class StructAlignTraceParameters {

	public boolean burnin; // True if this sample of parameters came from the burnin
	public boolean globalSigma; 
	
	public class PlottableParameter {
		public String name;
		public Number value;
		public double acceptanceRate;
		public int proposalCount;
		public int plotSide;
		public boolean wasProposed = true;
		public boolean fixedToParent;
		PlottableParameter (String n,Number v, double a, int p, int s, boolean f) {
			name = n;
			value = v;			
			acceptanceRate = a;
			proposalCount = p;
			plotSide = s;
			fixedToParent = f;
		}
	}
	public List<PlottableParameter> plottableParameters = new ArrayList<PlottableParameter>(); 

	public StructAlignTraceParameters(StructTrace s, boolean isBurnin){
		burnin = isBurnin;
		globalSigma = s.structAlign.globalSigma;
		for (McmcMove mcmcMove : s.structAlign.getMcmcMoves()) {
			if (mcmcMove instanceof ContinuousPositiveStructAlignMove) {
				if (((ContinuousPositiveStructAlignMove) mcmcMove).moveParams.isPlottable()) {
					plottableParameters.add(
							new PlottableParameter(
									mcmcMove.name,
									mcmcMove.getParam().get(),
									mcmcMove.acceptanceRate(),
									mcmcMove.proposalCount,
									((ContinuousPositiveStructAlignMove) mcmcMove).moveParams.plotSide(),
									((ContinuousPositiveStructAlignMove) mcmcMove).moveParams.fixedToParent()));							
									
				}
			}
		}
	}
	public void setProposalFlags(StructAlignTraceParameters previous) {
		if (previous.plottableParameters.size() == plottableParameters.size()) {
			for (int i=0; i<plottableParameters.size(); i++) {
				plottableParameters.get(i).wasProposed = 
					(plottableParameters.get(i).proposalCount != 
						previous.plottableParameters.get(i).proposalCount);
			}
			// We assume that the two lists being the same length implies
			// that the same set of parameters was recorded at the previous step.
			// This might not be true if somebody adds a method for altering
			// whether a parameter is plottable during the simulation.
		}
		else {
			throw new IllegalArgumentException("The number of plottable parameters has changed during the simulation.");
		}
	}
}
