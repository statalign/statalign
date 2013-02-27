package statalign.mcmc;

import java.util.ArrayList;
import java.util.List;

import statalign.base.Tree;
import statalign.base.Utils;

public abstract class McmcModule {
	
	/** Current log-likelihood contribution */
	public double curLogLike = 0;
	
	protected List<McmcMove> mcmcMoves = new ArrayList<McmcMove>();
	protected List<Integer> mcmcMoveWeights = new ArrayList<Integer>();
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
	public double logPrior() {
		return 0;
	}
	
	public void proposeParamChange(Tree tree) {
		int selectedMove = Utils.weightedChoose(mcmcMoveWeights);
		mcmcMoves.get(selectedMove).move(tree);
	}
	
	public void modifyProposalWidths() {
		for (McmcMove m : mcmcMoves) {
			if (!m.autoTune) { continue; }
			if (m.proposalCount > Utils.MIN_SAMPLES_FOR_ACC_ESTIMATE) {
				if (m.acceptanceRate() < Utils.MIN_ACCEPTANCE) {
					m.proposalWidthControlVariable *= Utils.SPAN_MULTIPLIER;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
				}
				else if (m.acceptanceRate() > Utils.MAX_ACCEPTANCE) {
					m.proposalWidthControlVariable /= Utils.SPAN_MULTIPLIER;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
				}
			}
		}
	}
	public abstract boolean isParamChangeAccepted(double logProposalRatio);
}
