package statalign.postprocess.utils;

import statalign.postprocess.utils.StructAlignTraceParameters;

public class StructAlignTraceParameterGetters {
	public interface ParameterGetter {
		abstract double getParameter(StructAlignTraceParameters s);
		abstract boolean wasProposed(StructAlignTraceParameters s);
	}
	public class Sigma2Getter implements ParameterGetter {
		private int index;
		public Sigma2Getter(int i) {
			index = i;
		}
		public double getParameter(StructAlignTraceParameters s) {
			return s.sigma2[index];
		}
		public boolean wasProposed(StructAlignTraceParameters s) {
			return s.sigma2Proposed[index];
		}
	}
	public class Sigma2HGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.sigma2Hier;
		}
		public boolean wasProposed(StructAlignTraceParameters s) {
			return s.sigma2HProposed;
		}
	}
	public class NuGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.nu;
		}
		public boolean wasProposed(StructAlignTraceParameters s) {
			return s.nuProposed;
		}
	}
	public class TauGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.tau;
		}
		public boolean wasProposed(StructAlignTraceParameters s) {
			return s.tauProposed;
		}
	}
	public class EpsilonGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.epsilon;
		}
		public boolean wasProposed(StructAlignTraceParameters s) {
			return s.epsilonProposed;
		}
	}
}
