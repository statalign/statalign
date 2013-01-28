package statalign.model.ext.plugins.structalign;

import statalign.model.ext.plugins.StructAlign;

// I wish Java allowed pointers and proper functors...

public class StructAlignParameterInterface {
	public interface ParameterInterface {
		abstract double get();
		abstract void set(double x);
	}
	public class Sigma2Interface implements ParameterInterface {
		private StructAlign structAlign;
		private int index;
		public Sigma2Interface(StructAlign s,int i) {
			structAlign = s;
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
		private StructAlign structAlign;
		public Sigma2HInterface(StructAlign s) {
			structAlign = s;
		}
		public double get() {
			return structAlign.sigma2Hier;
		}
		public void set(double x) {
			structAlign.sigma2Hier = x;
		}
	}
	public class NuInterface implements ParameterInterface {
		private StructAlign structAlign;
		public NuInterface(StructAlign s) {
			structAlign = s;
		}
		public double get() {
			return structAlign.nu;
		}
		public void set(double x) {
			structAlign.nu = x;
		}
	}
	public class TauInterface implements ParameterInterface {
		StructAlign structAlign;
		public TauInterface(StructAlign s) {
			structAlign = s;
		}
		public double get() {
			return structAlign.tau;
		}
		public void set(double x) {
			structAlign.tau = x;
		}
	}
	public class EpsilonInterface implements ParameterInterface {
		StructAlign structAlign;
		public EpsilonInterface(StructAlign s) {
			structAlign = s;
		}
		public double get() {
			return structAlign.epsilon;
		}
		public void set(double x) {
			structAlign.epsilon = x;
		}
	}
}
