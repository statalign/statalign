package statalign.postprocess.utils;

import statalign.postprocess.utils.StructAlignTraceParameters;

public class StructAlignTraceParameterGetters {
	public interface ParameterGetter {
		abstract double getParameter(StructAlignTraceParameters s);
	}
	public class Sigma2Getter implements ParameterGetter {
		private int index;
		public Sigma2Getter(int i) {
			index = i;
		}
		public double getParameter(StructAlignTraceParameters s) {
			return s.sigma2[index];
		}
	}
	public class Sigma2HGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.sigma2Hier;
		}
	}
	public class NuGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.nu;
		}
	}
	public class TauGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.tau;
		}
	}
	public class EpsilonGetter implements ParameterGetter {
		public double getParameter(StructAlignTraceParameters s) {
			return s.epsilon;
		}
	}
}
