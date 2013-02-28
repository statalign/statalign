package statalign.base.mcmc;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Tree;
import statalign.base.Utils;
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
		int ind;
		EdgeInterface (int i) {
			ind = i;
		}
		public double get() {
			return tree.vertex[ind].edgeLength;
		}
		public void set(double x) {
			tree.vertex[ind].edgeLength = x;
		}
	}
	
	public EdgeMove (McmcModule m, 
			int edgeIndex,
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,null,pr,prop,n);
		index = edgeIndex;
		param = new EdgeInterface(index);
		minValue = 0.01;
	}

	@Override
	public void move(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("EdgeMove.move must take an argument of type Tree.");
		}
//		if (Utils.DEBUG) {
//			System.out.println("EdgeMove"+index);
//		}
		((CoreMcmcModule) owner).getModelExtMan().beforeEdgeLenChange(tree,tree.vertex[index]);
		super.move(externalState);
		((CoreMcmcModule) owner).getModelExtMan().afterEdgeLenChange(tree,tree.vertex[index],lastMoveAccepted);
	}
	
	public void updateLikelihood(Object externalState) {
		if (param.get() >= minValue && param.get() < maxValue) {
			tree.vertex[index].edgeChangeUpdate();
			tree.vertex[index].calcAllUp();
			owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeEdgeLenChange(tree, tree.vertex[index]));
		}
	}
	@Override
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		tree.vertex[index].edgeChangeUpdate();
		tree.vertex[index].calcAllUp();
	}
}
