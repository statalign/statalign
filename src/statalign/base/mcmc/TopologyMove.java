package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;

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

	public void copyState(Object externalState) {
		// This is handled inside the Vertex
	}
	public double proposal(Object externalState) {
		int vertId, rnd = Utils.generator.nextInt(vnum - 3);
		vertId = tree.getTopVertexId(rnd);
		if (vertId != -1) {
			int lastId[] = new int[3], num = 0, newId = vertId;
	
			for (int i = vnum - 3; i < vnum; i++) {
				int id = tree.getTopVertexId(i);
				if (id == -1)
					lastId[num++] = i;
				else if (id < vertId)
					newId--;
			}
			rnd = lastId[newId];
		}
		nephew = tree.vertex[rnd];
		uncle = nephew.parent.brother();
		((CoreMcmcModule) owner).getModelExtMan().beforeTreeChange(tree, nephew);		
		// Should also do a beforeAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		double logProposalRatio = nephew.fastSwapWithUncle();
		return logProposalRatio;
	}
	public double logPriorDensity(Object externalState) {
		return 0.0;
	}
	public void updateLikelihood(Object externalState) {
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeTreeChange(tree, nephew));
	}
	public void restoreState(Object externalState) {
		uncle.fastSwapBackUncle();
	}
	
	public void move(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("AlignmentMove.copyState must take an argument of type Tree.");
		}
		vnum = tree.vertex.length;
		if (vnum <= 3) {
			return;
		}
		
		super.move(externalState);
		((CoreMcmcModule) owner).getModelExtMan().afterTreeChange(tree,lastMoveAccepted ? uncle : nephew,lastMoveAccepted);
	}
	
	 
	
}
