package statalign.base.mcmc;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Tree;
import statalign.base.Vertex;
import statalign.mcmc.ContinuousPositiveParameterMove;
import statalign.mcmc.McmcModule;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class EdgeMove extends ContinuousPositiveParameterMove {

	private int index;

	class EdgeInterface implements ParameterInterface {
		public double get() {
			return tree.vertex[index].edgeLength;
		}
		public void set(double x) {
			tree.vertex[index].edgeLength = x;
		}
	}
	
	public EdgeMove (McmcModule m, 
			int edgeIndex,
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,null,pr,prop,n);
		param = new EdgeInterface();
		setMinValue(0.01);
	}

	@Override
	public void move(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("IndelMove.move must take an argument of type Tree.");
		}
		((CoreMcmcModule) owner).getModelExtMan().beforeEdgeLenChange(tree,tree.vertex[index]);
		super.move(externalState);
		((CoreMcmcModule) owner).getModelExtMan().afterEdgeLenChange(tree,tree.vertex[index],lastMoveAccepted);
	}
	
	public void updateLikelihood(Object externalState) {
		if (param.get() > minValue && param.get() < maxValue) {
			tree.vertex[index].edgeChangeUpdate();
			tree.vertex[index].calcAllUp();
			owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeEdgeLenChange(tree, tree.vertex[index]));
		}
	}
	@Override
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		updateLikelihood(externalState);
		// Need to do the second step because updateLikelihood 
		// has repercussions, i.e. updates various stored values
	}
}
