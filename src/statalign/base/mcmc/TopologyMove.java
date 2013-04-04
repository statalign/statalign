package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;

public class TopologyMove extends McmcMove {

	Tree tree = null;
	int vnum;
	Vertex nephew;
	Vertex uncle;
	
	public TopologyMove (McmcModule m, String n) {
		owner = m;
		name = n;
		autoTune = false;
	}

	@Override
	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("TopologyMove.copyState must take an argument of type Tree.");
		}
		vnum = tree.vertex.length;
		if (vnum <= 3) { // then we can't do the nephew-uncle swap
			return;
		}
		
		int rnd = Utils.generator.nextInt(vnum);
		int vertId = tree.getTopVertexId(rnd);
		// Need to choose a valid nephew, i.e. not the root or either of its children
		while (vertId >= 0) { // while nephew is invalid
			rnd = Utils.generator.nextInt(vnum); // select a new nephew
			vertId = tree.getTopVertexId(rnd);
		}
		

		nephew = tree.vertex[rnd];
		uncle = nephew.parent.brother();
		((CoreMcmcModule) owner).getModelExtMan().beforeTreeChange(tree, nephew);		
		// Should also do a beforeAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		// The rest of the state copying is handled inside the Vertex
	}
	@Override
	public double proposal(Object externalState) {
		double logProposalRatio = nephew.fastSwapWithUncle();
		// Below is another version, slow and slightly better mixing
		// double logProposalRatio = nephew.swapWithUncleAlignToParent();
		return logProposalRatio;
	}
	@Override
	public double logPriorDensity(Object externalState) {
		return 0.0;
	}
	@Override
	public void updateLikelihood(Object externalState) {
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeTreeChange(tree, nephew));
	}
	@Override
	public void restoreState(Object externalState) {
		uncle.fastSwapBackUncle();
		// If using the alternative move:
        // uncle.swapBackUncleAlignToParent();
	}
	
	@Override
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterTreeChange(tree,lastMoveAccepted ? uncle : nephew,lastMoveAccepted);
		// Should also do an afterAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		
		if (Utils.DEBUG) {
			tree.checkPointers();
		}
	}
	
	 
	
}
