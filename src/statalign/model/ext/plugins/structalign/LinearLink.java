package statalign.model.ext.plugins.structalign;

import statalign.utils.LinkFunction;

public class LinearLink implements LinkFunction<Double> {
	public Double f(Double t) {
		return t;
	}
}