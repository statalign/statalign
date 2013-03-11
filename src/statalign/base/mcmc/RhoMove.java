package statalign.base.mcmc;

import statalign.mcmc.McmcModule;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;

public class RhoMove extends IndelMove {
	
	public RhoMove (McmcModule m, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,pr,prop,n);	
		param = new RhoInterface();
	}
}
