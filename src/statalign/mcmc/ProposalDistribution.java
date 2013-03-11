package statalign.mcmc;

public abstract class ProposalDistribution <T> {

	public abstract double logDensity(T x);
	public abstract T sample();
	public abstract void updateProposal(double proposalWidthControlVariable, 
			T currentParam);
}
