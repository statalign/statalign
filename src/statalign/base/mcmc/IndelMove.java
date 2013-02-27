package statalign.base.mcmc;

import java.util.List;
import java.util.ArrayList;

import statalign.mcmc.ContinuousPositiveParameterMove;
import statalign.mcmc.McmcModule;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public abstract class IndelMove extends ContinuousPositiveParameterMove {
	
	class RInterface implements ParameterInterface {
		public double get() {
			return tree.hmm2.params[0];
		}
		public void set(double x) {
			tree.hmm2.params[0] = x;
		}
	}
	class RhoInterface implements ParameterInterface {
		// rho = lambda + mu
		public double get() {
			double lambda = tree.hmm2.params[1];
			double mu = tree.hmm2.params[2];
			return(lambda + mu);
		}
		public void set(double x) {
			double lambda = tree.hmm2.params[1];
			double mu = tree.hmm2.params[2];
			double theta = lambda / (lambda + mu);
			// Now change rho, keeping theta fixed
			tree.hmm2.params[1] = x * theta;
			tree.hmm2.params[2] = x - tree.hmm2.params[1];
		}
	}
	class ThetaInterface implements ParameterInterface {
		// theta = lambda / (lambda + mu)
		public double get() {
			double lambda = tree.hmm2.params[1];
			double mu = tree.hmm2.params[2];
			return(lambda/(lambda + mu));
		}
		public void set(double x) {
			double lambda = tree.hmm2.params[1];
			double mu = tree.hmm2.params[2];
			double rho = lambda + mu;
			// Now change theta, keeping rho fixed
			tree.hmm2.params[1] = x * rho;
			tree.hmm2.params[2] = rho - tree.hmm2.params[1];
		}
	}
	class PhiInterface implements ParameterInterface {
		// phi = lambda / mu
		public double get() {
			double lambda = tree.hmm2.params[1];
			double mu = tree.hmm2.params[2];
			return(lambda/mu);
		}
		public void set(double x) {
			// Now change phi, keeping lambda fixed.
			// i.e. set mu' = lambda / phi'
			tree.hmm2.params[2] = tree.hmm2.params[1] / x;
		}
	}
	class LambdaInterface implements ParameterInterface {
		public double get() {
			return tree.hmm2.params[1];
		}
		public void set(double x) {
			tree.hmm2.params[1] = x;
		}
	}
	class MuInterface implements ParameterInterface {
		public double get() {
			return tree.hmm2.params[2];
		}
		public void set(double x) {
			tree.hmm2.params[2] = x;
		}
	}
	
	public IndelMove (McmcModule m, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,null,pr,prop,n);		
	}

	public void move(Object externalState) {
		owner.getModelExtMan().beforeIndelParamChange(tree,tree.hmm2,this);
		super.move(externalState);
		owner.getModelExtMan().afterIndelParamChange(tree,tree.hmm2,this);
	}
	public void updateLikelihood(Object externalState) {
		if (param.get() > minValue && param.get() < maxValue) {
			for (int i = 0; i < tree.vertex.length; i++) {
				tree.vertex[i].updateHmmMatrices();
			}
			tree.root.calcIndelLikeRecursively();
			owner.setLogLike(owner.getModelExtMan().logLikeIndelParamChange(tree, tree.hmm2, this));

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
