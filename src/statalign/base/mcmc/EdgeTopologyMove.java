package statalign.base.mcmc;

import statalign.base.AlignColumn;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.GammaPrior;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import statalign.mcmc.MultiplicativeProposal;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.mcmc.UniformProposal;

public class EdgeTopologyMove extends McmcMove {

	Tree tree = null;
	int vnum;
	Vertex nephew;
	Vertex uncle;
	EdgeMove nephewEdgeMove;
	EdgeMove uncleEdgeMove;
	
	public EdgeTopologyMove (McmcModule m, String n) {
		owner = m;
		name = n;
		autoTune = false; 
		// Because in this case making the variance higher is
		// just as likely to increase acceptance.
		nephewEdgeMove = new EdgeMove(owner,0,new GammaPrior(1,1),
				//new MultiplicativeProposal(),"NephewEdge");
				new UniformProposal(),"NephewEdge");
		uncleEdgeMove = new EdgeMove(owner,0,new GammaPrior(1,1),
				//new MultiplicativeProposal(),"UncleEdge");
				new UniformProposal(),"UncleEdge");
	}

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
		if (vnum <= 3) {
			return;
		}
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
		nephewEdgeMove.setEdgeIndex(nephew.index);
		uncleEdgeMove.setEdgeIndex(uncle.index);
		nephewEdgeMove.copyState(externalState); // calls beforeEdgeLenChange
		uncleEdgeMove.copyState(externalState); // calls beforeEdgeLenChange
		//((CoreMcmcModule) owner).getModelExtMan().beforeTreeChange(tree, nephew);		
		// Should also do a beforeAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		// The rest of the state copying is handled inside the Vertex
	}
	public double proposal(Object externalState) {
		nephewEdgeMove.proposalWidthControlVariable = proposalWidthControlVariable;
		uncleEdgeMove.proposalWidthControlVariable = proposalWidthControlVariable;
//		double meanLength = (nephewEdgeMove.getParam().get()+ uncleEdgeMove.getParam().get())/2;
//		nephewEdgeMove.updateProposal(proposalWidthControlVariable,meanLength);
//		uncleEdgeMove.updateProposal(proposalWidthControlVariable,meanLength);
//		nephewEdgeMove.proposalWidthControlVariable = 
//			proposalWidthControlVariable;// * edgeDiff;
//		uncleEdgeMove.proposalWidthControlVariable = 
//			proposalWidthControlVariable;// * edgeDiff;
		
//		System.out.println("Before proposal: "+nephewEdgeMove.getEdgeIndex()+"="+nephewEdgeMove.getParam().get()
//				+"\t"+uncleEdgeMove.getEdgeIndex()+"="+uncleEdgeMove.getParam().get());
		double logProposalRatio = nephewEdgeMove.proposal(externalState);
		logProposalRatio += uncleEdgeMove.proposal(externalState);
//		System.out.println("After proposal: "+nephewEdgeMove.getParam().get()+"\t"+uncleEdgeMove.getParam().get());
		logProposalRatio += nephew.fastSwapWithUncle();
		// Below is another version, slow and slightly better mixing
		// double logProposalRatio = nephew.swapWithUncleAlignToParent();
		return logProposalRatio;
	}
	public double logPriorDensity(Object externalState) {
		return nephewEdgeMove.logPriorDensity(externalState) + uncleEdgeMove.logPriorDensity(externalState);
	}
	public void updateLikelihood(Object externalState) {
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeTreeChange(tree, nephew));
	}
	public void restoreState(Object externalState) {
		uncle.fastSwapBackUncle();
		// If using the alternative move:
        // uncle.swapBackUncleAlignToParent();
		//System.out.println("Before restore: "+nephewEdgeMove.getParam().get()+"\t"+uncleEdgeMove.getParam().get());
		nephewEdgeMove.restoreState(externalState);
		uncleEdgeMove.restoreState(externalState);
		//System.out.println("After restore: "+nephewEdgeMove.getParam().get()+"\t"+uncleEdgeMove.getParam().get());
	}
	
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterTreeChange(tree,lastMoveAccepted ? uncle : nephew,lastMoveAccepted);
		// Should also do an afterAlignChange here, but not obvious what to pass
		// as the selectedRoot argument.
		
		if (Utils.DEBUG) {
			for (int i = 0; i < tree.vertex.length; i++) {
				if (tree.vertex[i].left != null && tree.vertex[i].right != null) {
					tree.vertex[i].checkPointers();
					AlignColumn p;
					// checking pointer integrity
					for (AlignColumn c = tree.vertex[i].left.first; c != null; c = c.next) {
						p = tree.vertex[i].first;
						while (c.parent != p && p != null)
							p = p.next;
						if (p == null)
							throw new Error(
									"children does not have a parent!!!"
											+ tree.vertex[i] + " "
											+ tree.vertex[i].print());
					}
					for (AlignColumn c = tree.vertex[i].right.first; c != null; c = c.next) {
						p = tree.vertex[i].first;
						while (c.parent != p && p != null)
							p = p.next;
						if (p == null)
							throw new Error(
									"children does not have a parent!!!"
											+ tree.vertex[i] + " "
											+ tree.vertex[i].print());
					}
	
				}
			}
		}
	}
	
	 
	
}
