package statalign.model.ext;

public interface PriorDistribution<T> {
	public abstract double logDensity(Double x);
	//TODO implement this for arbitrary arguments
}
