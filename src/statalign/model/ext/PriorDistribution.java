package statalign.model.ext;

public interface PriorDistribution<T> {
	public abstract double logDensity(Double x);
	public abstract double logDensityUnnormalised(Double x);
	//TODO implement this for arbitrary arguments
}
