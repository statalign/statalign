package statalign.model.ext.plugins.structalign;

import statalign.mcmc.ParameterInterface;
import statalign.model.ext.plugins.StructAlign;

// I wish Java allowed pointers and proper functors...

public class StructAlignParameterInterface {
	private StructAlign structAlign;
	public StructAlignParameterInterface(StructAlign s) {
		structAlign = s;
	}
	public class Sigma2Interface implements ParameterInterface {
		private int index;
		public int getIndex() {
			return index;
		}
		public Sigma2Interface(int i) {
			index = i;
		}
		public double get() {
			return structAlign.sigma2[index];
		}
		public void set(double x) {
			structAlign.sigma2[index] = x;
		}
	}
	public class Sigma2HInterface implements ParameterInterface {
		public double get() {
			return structAlign.sigma2Hier;
		}
		public void set(double x) {
			structAlign.sigma2Hier = x;
		}
	}
	public class NuInterface implements ParameterInterface {
		public double get() {
			return structAlign.nu;
		}
		public void set(double x) {
			structAlign.nu = x;
		}
	}
	public class TauInterface implements ParameterInterface {
		public double get() {
			return structAlign.tau;
		}
		public void set(double x) {
			structAlign.tau = x;
		}
	}
	public class EpsilonInterface implements ParameterInterface {
		public double get() {
			return structAlign.epsilon;
		}
		public void set(double x) {
			structAlign.epsilon = x;
		}
	}
}
