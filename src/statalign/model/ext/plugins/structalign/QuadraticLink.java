package statalign.model.ext.plugins.structalign;

import statalign.utils.LinkFunction;

public class QuadraticLink implements LinkFunction<Double> {
	public Double f(Double t) {
		return Math.pow(t,2);
	}
}