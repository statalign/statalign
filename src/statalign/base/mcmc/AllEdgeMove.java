package statalign.base.mcmc;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.ContinuousPositiveParameterMove;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;

public class AllEdgeMove extends ContinuousPositiveParameterMove {

	double multiplier = 1.0;
	class AllEdgeInterface implements ParameterInterface {
		public double get() {
			return 1;
		}
		public void set(double x) {
			multiplier = x;
			for (int i=0; i<tree.vertex.length-1; i++) {
				tree.vertex[i].edgeLength *= x;
			}
		}
	}
	
	public AllEdgeMove (McmcModule m, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,null,pr,prop,n);
		param = new AllEdgeInterface();
		minValue = 0.01;
	}
	@Override
	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("AllEdgeMove.move must take an argument of type Tree.");
		}
		((CoreMcmcModule) owner).getModelExtMan().beforeEdgeLenChange(tree,tree.vertex[0]);
		// Should pass the index of a tip Vertex?
		super.copyState(externalState);
	}
	@Override
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterEdgeLenChange(tree,tree.vertex[0],lastMoveAccepted);
		// Should pass the index of a tip Vertex?
	}
	@Override
	public double proposal(Object externalState) {
		double logProposalRatio = super.proposal(externalState);
		for (int i=0; i<tree.vertex.length-1; i++) {
			tree.vertex[i].edgeChangeUpdate();
		}
//		tree.root.edgeChangeUpdateRecursively();
		tree.root.calcFelsRecursively();
		return logProposalRatio;
	}
	@Override
	public double logPriorDensity(Object externalState) {
		double logDensity = 0;
		for (int i=0; i<tree.vertex.length-1; i++) {
			logDensity +=  prior.logDensityUnnormalised(tree.vertex[i].edgeLength);
		}
		return logDensity;
	}
	public void updateLikelihood(Object externalState) {
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeEdgeLenChange(tree, tree.vertex[0]));
	    // Should pass the index of a tip Vertex?
	}
	@Override
	public void restoreState(Object externalState) {		
		owner.setLogLike(oldll);
		for (int i=0; i<tree.vertex.length-1; i++) {
			tree.vertex[i].edgeLength /= multiplier;
			tree.vertex[i].edgeChangeUpdate();
		}
//		tree.root.edgeChangeUpdateRecursively();
		tree.root.calcFelsRecursively();
	}
}
