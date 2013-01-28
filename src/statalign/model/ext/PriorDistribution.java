package statalign.model.ext;

public interface PriorDistribution<T> {
	public abstract double logDensity(T x);
}
