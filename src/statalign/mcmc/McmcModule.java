package statalign.mcmc;

import java.util.ArrayList;
import java.util.List;

import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.base.Utils;

public abstract class McmcModule {
	
	protected Mcmc mcmc;
	public void setMcmc(Mcmc m) {
		mcmc = m;
	}
	
	/** Current log-likelihood contribution */
	public double curLogLike = 0;
	
	protected List<McmcMove> mcmcMoves = new ArrayList<McmcMove>();
	protected List<Integer> mcmcMoveWeights = new ArrayList<Integer>();
	public int getParamChangeWeight() {
		int w = 0;
		for (int i=0; i<mcmcMoveWeights.size(); i++) {
			w += mcmcMoveWeights.get(i);
		}
		return w;
	}
	public void addMcmcMove(McmcMove m, int weight) {
		mcmcMoves.add(m);
		mcmcMoveWeights.add(weight);
	}
	public List<McmcMove> getMcmcMoves() {
		return mcmcMoves;
	}
	public void setAllMovesNotProposed() {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.moveProposed = false;
		}
	}
	public McmcMove getMcmcMove(String name) {
		for (McmcMove mcmcMove : mcmcMoves) {
			if (mcmcMove.name.equals(name)) {
				return mcmcMove;
			}
		}
		throw new RuntimeException("McmcMove "+name+" not found.");
	}
	
	/**
	 * Called before the start of MCMC sampling, but after the initial tree, alignment etc. have been
	 * generated. Override to initialise data structures etc.
	 * @param tree the starting tree
	 */
	public void beforeSampling(Tree tree) {}
	
	public void afterSampling() {}
	
	/**
	 * This should return the log of the model's contribution to the likelihood, it will be added on to
	 * the log-likelihood of the current point in the MCMC state space. Normally it will be called once at the
	 * initialisation of the MCMC process and from then on once in each MCMC step, when proposing any change.
	 * In debug mode, will be called more often (including after proposed changes) to ensure consistency.
	 * @param tree current tree
	 * @return log of model extension likelihood, conditional on current tree, alignment and params
	 */
	public abstract double logLikeFactor(Tree tree);
	
	public double getLogLike() {
		return curLogLike;
	}
	public void setLogLike(double ll) {
		curLogLike = ll;
	}

	/**
	 * This should return the log of the total prior calculated for the model parameters. It is only used
	 * in parallel mode when proposing swaps between chains. By default returns 0.
	 */
	public double logPrior(Tree tree) {
		return 0;
	}
	
	public boolean proposeParamChange(Tree tree) {
		int selectedMoveIndex = Utils.weightedChoose(mcmcMoveWeights);
		McmcMove selectedMove = mcmcMoves.get(selectedMoveIndex); 
		selectedMove.move(tree);
		return selectedMove.lastMoveAccepted;
	}
	
	public void modifyProposalWidths() {
		for (McmcMove m : mcmcMoves) {
			if (!m.autoTune) { continue; }
			if (m.proposalCount > Utils.MIN_SAMPLES_FOR_ACC_ESTIMATE) {
				if (m.acceptanceRate() < m.minAcceptance) {
					m.proposalWidthControlVariable *= m.spanMultiplier;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
				}
				else if (m.acceptanceRate() > m.maxAcceptance) {
					m.proposalWidthControlVariable /= m.spanMultiplier;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
				}
			}
		}
	}
	public boolean isParamChangeAccepted(double logProposalRatio) {
		return mcmc.isParamChangeAccepted(logProposalRatio);
	}
}
